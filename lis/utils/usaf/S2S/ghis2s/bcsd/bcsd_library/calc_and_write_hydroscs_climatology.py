#!/usr/bin/env python
"""
# Author: Shrad Shukla
# Usage: This is a part of the BCSD code. It creates sorted time series of observed data
# This has been created seperately so during the bias-correction process we don't need to calculate
# This module bias corrects a forecasts following probability mapping approach as described in
# Wood et al. 2002
# Date: August 06, 2015
"""

import os
from glob import glob
import sys
import gc
import numpy as np
import xarray as xr
import yaml
from ghis2s.shared.utils import get_domain_info, load_ncdata, get_chunk_sizes, write_ncfile
from ghis2s.shared.logging_utils import TaskLogger

# This function takes in a time series as input and provides sorted times series of values and
# quantiles in return
## Using Plotting position formula to get quantiles

VAR_REPLACE = {'LWdown': 'LWGAB',
               'Rainf': 'PRECTOT',
               'Psurf': 'PS',
               'Qair': 'QV2M',
               'SWdown': 'SWGDN',
               'Tair': 'T2M',
               'Wind': 'U10M'}

def magnitude(_a, _b):
    ''' computes wind magnitude u^2 + v^2'''
    func = lambda x, y: np.sqrt(x**2 + y**2)
    return xr.apply_ufunc(func, _a, _b)

def create_sorted_ts(data_ts):
    """create sorted ts"""
    clim_sort = np.sort(data_ts) ## np.sort a function to sort data into ascending order

    ## Now storing plotting position in an array
    ## formula for plotting position is (1/(n+1))
    quant_ts = np.arange(1, len(data_ts)+1)/(len(data_ts)+1)
    ## np.arange creates an array where values range from 1 to the length of time series and then
    ## simply dividing every value by length of the time series (n) + 1
    ## Now passing on sorted observed and forecast value and quantile time series
    return clim_sort, quant_ts

# Directory and file addresses
CMDARGS = str(sys.argv)
# variable name:
#  1) if 'LWdown', 'Rainf', 'Psurf', 'Qair', 'SWdown', 'Tair', or 'Wind' creates observational climatology
#  2) 1 or 2: writes monthly files for the 6-month segment of the year str(sys.argv[4])
VAR_NAME = str(sys.argv[1]) #
CONFIGFILE = str(sys.argv[2])
OUTDIR =  str(sys.argv[3])
with open(CONFIGFILE, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)

INDIR = config['SETUP']['AF10KM']
CLIM_SYR = config['BCSD']['clim_start_year']
CLIM_EYR = config['BCSD']['clim_end_year']

LATS, LONS = get_domain_info(CONFIGFILE, coord=True)
lat1, lat2, lon1, lon2 = get_domain_info(CONFIGFILE, extent=True)
lat_chunk, lon_chunk = get_chunk_sizes(None, dim_in=[len(LATS), len(LONS)])
COMP = {'zlib':True, 'complevel':6, 'shuffle':True, 'missing_value': -9999., '_FillValue': -9999.}

def write_monthly_files(segment):
    ''' writes monthly file '''
    main_task = str(sys.argv[4]) + '-' + segment
    task_name = os.environ.get('SCRIPT_NAME')
    logger = TaskLogger(task_name,
                        os.getcwd(),
                        f'bcsd/bcsd_library/calc_and_write_hydroscs_climatology.py: {main_task}')
    infile_template = '{}/{:04d}{:02d}/HydroSCS_{:04d}{:02d}*.d01.nc'
    outfile_template = '{}/HydroSCS_{:04d}{:02d}.nc'
    
    year = int(sys.argv[4])
    if segment == '1':
        month_range = [0, 3]
    if segment == '2':
        month_range = [3, 6]
    if segment == '3':
        month_range = [6, 9]
    if segment == '4':
        month_range = [9, 12]
        
    for mon in range(month_range[0], month_range[1]):
        subtask = main_task + f' month: {mon+1:02d}'
        infile = infile_template.format(INDIR, year, mon+1, year, mon+1)
        outfile = outfile_template.format(OUTDIR, year, mon+1)
        logger.info(f"Reading Observed Data {infile}",  subtask=subtask)
        file_list = sorted(glob(infile))
        datasets = []
        for i, file in enumerate(file_list):
            logger.info(f"Reading Observed Data {file}",  subtask=subtask)
            ds = xr.open_dataset(file, engine='netcdf4', chunks={'lat': 'auto', 'lon': 'auto'})
            if 'time' in ds.dims:
                ds = ds.drop_dims('time')
            if 'time' in ds.variables:
                ds = ds.drop_vars('time')
            ds = ds.expand_dims('time')
            datasets.append(ds)
        print(xr.concat(datasets,  dim='time'))
        ds = xr.concat(datasets,  dim='time').mean(dim='time')
        encoding = {var: COMP for var in ds.data_vars}
        logger.info(f"Writing merged output to {outfile}", )
        write_ncfile(ds, outfile, encoding, [logger, subtask])
        ds.close()
        del ds
        for ds in datasets:
            ds.close()
            del ds
        gc.collect()
        
if VAR_NAME.isdigit():
    write_monthly_files(VAR_NAME)
    sys.exit()

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
SUBTASK = VAR_NAME
task_name = os.environ.get('SCRIPT_NAME')
logger = TaskLogger(task_name,
                    os.getcwd(),
                    f'bcsd/bcsd_library/calc_and_write_hydroscs_climatology.py: {VAR_NAME}')
INFILE_TEMPLATE = '{}/HydroSCS_{:04d}{:02d}.nc'
OUTFILE_TEMPLATE = '{}/{}_obs_clim.nc'

# Read in LDT landmask - impose on precip field prior to BC:
ldt_xr = xr.open_dataset(config['SETUP']['supplementarydir'] + '/lis_darun/' + \
        config['SETUP']['ldtinputfile'])
mask = np.array(ldt_xr['LANDMASK'].values)

## Defining array to store observed data
OBS_DATA_COARSE = np.empty((((CLIM_EYR-CLIM_SYR)+1)*12, len(LATS), len(LONS)))

MON_COUNTER = 0
for YEAR in range(CLIM_SYR, CLIM_EYR+1):
    for MON in range(0, 12):
        INFILE = INFILE_TEMPLATE.format(INDIR, YEAR, MON+1)
        logger.info(f"Reading Observed Data {INFILE}",  subtask=SUBTASK)
        ds = xr.open_dataset(INFILE, engine='netcdf4', chunks={'lat': "auto", 'lon': "auto"})
        if VAR_NAME == 'Wind':
            OBS_DATA_COARSE[MON_COUNTER, ] = \
                magnitude(ds['Wind_E'].load(), ds['Wind_N'].load()).values
        else:
            OBS_DATA_COARSE[MON_COUNTER, ] = \
                ds[VAR_NAME].values
#       Impose mask on precip values:
        if VAR_NAME == 'Rainf':
            OBS_DATA_COARSE[MON_COUNTER,mask == 0] = -9999.
        MON_COUNTER+=1
        ds.close()
        del ds
        gc.collect()

## Looping through each month and creating time series of quantiles and observed climatology
## (i.e. sorted observed values) for each grid cells with the mask
CLIM_ARRAY = np.full((13, (CLIM_EYR-CLIM_SYR)+1, len(LATS), len(LONS)), -9999., dtype=np.float32)
for lat_num in range(0, len(LATS)):
    for lon_num in range(0, len(LONS)):
        ## Only work with grid cells that are within the given mask
        if ((lat1<=LATS[lat_num]) and (LATS[lat_num]<=lat2) and
            (lon1<=LONS[lon_num]) and (LONS[lon_num]<=lon2)):
            # Check if this is a land pixel
            if mask[lat_num, lon_num] > 0:
                ## 1st column is for sorted quantile array
                ## The other twelve columns have sorted climatology values for all years
                for MON_NUM in range(0, 12): ## Looping from month 1 to 12
                    ## First storing time series for the given month, latitude and longitude
                    OBS_TS = OBS_DATA_COARSE[MON_NUM::12, lat_num, lon_num]
                    
                    # Filter out fill values before creating sorted time series
                    valid_data = OBS_TS[OBS_TS != -9999.]
                    
                    if len(valid_data) > 0:
                        ## Now creating sorted time series of observed data and a time series of
                        ## quantile values (ranging from 1/(n+1) to n/(n+1)
                        clim_sort, quant_ts = create_sorted_ts(valid_data)
                        
                        # For precipitation, replace any zeros with small value
                        if VAR_NAME == 'Rainf':
                            clim_sort = np.where(clim_sort == 0, 1.e-7, clim_sort)
                        
                        # Pad arrays if needed to match expected length
                        expected_len = (CLIM_EYR-CLIM_SYR)+1
                        if len(clim_sort) < expected_len:
                            # Pad with the last valid value
                            clim_sort = np.pad(clim_sort, (0, expected_len - len(clim_sort)), 
                                             mode='edge')
                            quant_ts = np.pad(quant_ts, (0, expected_len - len(quant_ts)), 
                                            mode='edge')
                        
                        CLIM_ARRAY[MON_NUM+1, :, lat_num, lon_num] = clim_sort
                        CLIM_ARRAY[0, :, lat_num, lon_num] = quant_ts

## finished storing climatology for all months
OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, VAR_REPLACE.get(VAR_NAME))
## opening outfile to write
print (CLIM_ARRAY.min(), CLIM_ARRAY.max())
# Now writing output file
CLIM_XR = xr.Dataset()
CLIM_XR['clim'] = (('DIST', 'time', 'latitude', 'longitude'), CLIM_ARRAY)
CLIM_XR.coords['DIST'] = (('DIST'), np.arange(13))
CLIM_XR.coords['time'] = (('time'), np.arange((CLIM_EYR-CLIM_SYR)+1))
CLIM_XR.coords['latitude'] = (('latitude'), LATS)
CLIM_XR.coords['longitude'] = (('longitude'), LONS)

#Now selecting only the data over the given lat lon box
SLICED_CLIM_XR = CLIM_XR.sel(longitude=slice(lon1, lon2), latitude=slice(lat1, lat2))
edict = {}
for var in SLICED_CLIM_XR.data_vars:
    edict[var] = {"zlib":True, "complevel":6, "shuffle":True, "missing_value": np.nan, "_FillValue": np.nan}
SLICED_CLIM_XR.to_netcdf(OUTFILE, format="NETCDF4", encoding=edict)
