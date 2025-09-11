#!/usr/bin/env python
"""
# Author: Shrad Shukla
# coding: utf-8
#Author: Shrad Shukla
#Usage: This is a module for the BCSD code.
#This module bias corrects a forecasts following
#probability mapping approach as described in Wood et al. 2002
#Date: August 06, 2015
"""



import calendar
import sys
from datetime import datetime
import os
from dateutil.relativedelta import relativedelta
import xarray as xr
import numpy as np
import yaml
# pylint: disable=import-error
from ghis2s.bcsd.bcsd_library.shrad_modules import read_nc_files
from ghis2s.bcsd.bcsd_library.bcsd_stats_functions import write_4d_netcdf
from ghis2s.shared.utils import get_domain_info
from ghis2s.bcsd.bcsd_library import bcsd_function
from ghis2s.bcsd.bcsd_library.bcsd_function import VarLimits as lim
from ghis2s.shared.logging_utils import TaskLogger
# pylint: enable=import-error

limits = lim()
CF2VAR = {
    'PRECTOT': 'PRECTOT',
    'LWS': 'LWS',
    'SLRSF': 'SLRSF',
    'PS': 'PS',
    'Q2M':'Q2M',
    'T2M': 'T2M',
    'WIND10M': 'WIND',
    }

## Usage: <Name of variable in observed climatology>
## <Name of variable in reforecast climatology
## (same as the name in target forecast> <forecast model number>
CMDARGS = str(sys.argv)
OBS_VAR = str(sys.argv[1])
FCST_VAR = str(sys.argv[2])
BC_VAR = str(sys.argv[3])
## This is used to figure out if the variable is a precipitation variable or not
UNIT = str(sys.argv[4])
INIT_FCST_MON = int(sys.argv[5])

task_name = os.environ.get('SCRIPT_NAME')
subtask = f'{FCST_VAR}'
logger = TaskLogger(task_name,
                    os.getcwd(),
                    f'bcsd/bcsd_library/bias_correction_modulefast.py processing {FCST_VAR} for month {INIT_FCST_MON:02d}')

# Forecast model and ensemble input arguments:
LEAD_FINAL = int(sys.argv[6])
ENS_NUM = int(sys.argv[7])

FCST_SYR = int(sys.argv[8])
TARGET_FCST_SYR = int(sys.argv[8])
TARGET_FCST_EYR = int(sys.argv[9])
CLIM_SYR = int(sys.argv[10])
CLIM_EYR = int(sys.argv[11])

# Directory and file addresses
OBS_INDIR = str(sys.argv[12])
FCST_INDIR = str(sys.argv[13])

# Observation climatology filename templates:
OBS_CLIM_FILE_TEMPLATE = '{}/raw/Climatology/{}_obs_clim.nc'
FCST_CLIM_FILE_TEMPLATE = '{}/raw/Climatology/{}/{}_fcst_clim.nc'
MONTH_NAME_TEMPLATE = '{}01'
# CFSV2 monthly filename template:
FCST_INFILE_TEMPLATE = '{}/raw/Monthly/{}/{:04d}/ens{:01d}/{}.cfsv2.{:04d}{:02d}.nc'

CONFIG_FILE = str(sys.argv[14])
LAT1, LAT2, LON1, LON2 = get_domain_info(CONFIG_FILE, extent=True)
LATS, LONS = get_domain_info(CONFIG_FILE, coord=True)
with open(CONFIG_FILE, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)
ldt_xr = xr.open_dataset(config['SETUP']['supplementarydir'] + '/lis_darun/' + \
        config['SETUP']['ldtinputfile'])
land_mask = np.array(ldt_xr['LANDMASK'].values)

### Output directory
OUTFILE_TEMPLATE = '{}/{}.CFSv2.{}_{:04d}_{:04d}.nc'
OUTDIR = str(sys.argv[15])

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR, exist_ok=True)

NUM_YRS = (CLIM_EYR-CLIM_SYR)+1
TINY = ((1/(NUM_YRS))/ENS_NUM)/2   # Adjust quantile, if it is out of bounds
             #          This value represents 1/NYRS/NENS/2, so about
             #          half the prob. interval beyond the lowest value
             #          (arbitrary choice) */
	     ## This is probably used for real-time forecasts when a
             ## forecasted value happened to be an outlier of the reforecast
             ## climatology

EPS = 1.0e-5

##### Starting bias-correction from here

# First read observed climatology for the given variable
OBS_CLIM_FILE = OBS_CLIM_FILE_TEMPLATE.format(OBS_INDIR, OBS_VAR)
OBS_CLIM_ARRAY = xr.open_dataset(OBS_CLIM_FILE)

def get_index(ref_array, my_value):
    """
      Function for extracting the index of a Numpy array (ref_array)
      which value is closest to a given number.

      Input parameters:
        - ref_array: reference Numpy array
        - my_value:  floating point number

      Returned value:
        - An integer corresponding to the index
    """
    return np.abs(ref_array - my_value).argmin()

def monthly_calculations(mon):
    """Monthly Bias Correction"""
    month_name = MONTH_NAME_TEMPLATE.format((calendar.month_abbr[mon]).lower())
    ## This provides abbrevated version of the name of a month:
    ## (e.g. for January (i.e. Month number = 1) it will return "Jan").
    ##The abbrevated name is used in the forecasts file name
    print(f"Forecast Initialization month is {month_name}")
    #First read forecast climatology for the given variable and forecast
    #initialzation month
    fcst_clim_infile = FCST_CLIM_FILE_TEMPLATE.format(FCST_INDIR,\
                                                      month_name, FCST_VAR)
    logger.info(f"Reading forecast climatology {fcst_clim_infile}", subtask=subtask)
    fcst_clim_array = xr.open_dataset(fcst_clim_infile)
    #First read raw forecasts
    fcst_coarse = np.empty(((TARGET_FCST_EYR-TARGET_FCST_SYR)+1,
                            LEAD_FINAL, ENS_NUM, len(LATS), len(LONS)))
    for lead_num in range(0, LEAD_FINAL): ## Loop from lead =0 to Final Lead
        for ens in range(ENS_NUM):
            for init_fcst_year in range(TARGET_FCST_SYR, TARGET_FCST_EYR+1):
                ## Reading forecast file
                fcst_date = datetime(init_fcst_year, INIT_FCST_MON, 1) + \
                            relativedelta(months=lead_num)
                fcst_year, fcst_month = fcst_date.year, fcst_date.month
                infile = FCST_INFILE_TEMPLATE.format(FCST_INDIR, month_name, \
                init_fcst_year, ens+1, month_name, fcst_year, fcst_month)
                logger.info(f"Reading forecast file {infile}", subtask=subtask)
                fcst_coarse[init_fcst_year-TARGET_FCST_SYR, lead_num, ens, ] = \
                read_nc_files(infile, FCST_VAR)[:]
    # Defining array to store bias-corrected monthly forecasts
    #correct_fcst_coarse = np.ones(((TARGET_FCST_EYR-TARGET_FCST_SYR)+1, \
    #LEAD_FINAL, ENS_NUM, len(LATS), len(LONS)))*-9999.
    #print("shape of fcst_coarse: ", fcst_coarse.shape)

    # Get the lat/lon indexes for the ranges
    ilat_min = get_index(LATS, LAT1)
    ilat_max = get_index(LATS, LAT2)
    ilon_min = get_index(LONS, LON1)
    ilon_max = get_index(LONS, LON2)
    nlats = len(LATS)
    nlons = len(LONS)

    # Get the values (Numpy array) for the lat/lon ranges
    np_obs_clim_array = OBS_CLIM_ARRAY.clim.sel(longitude=slice(LON1, LON2), \
                        latitude=slice(LAT1, LAT2)).values
    np_fcst_clim_array = fcst_clim_array.clim.sel(longitude=slice(LON1, LON2), \
                         latitude=slice(LAT1, LAT2)).values

    #print("Latitude:  ", nlats, ilat_min, ilat_max)
    #print("Longitude: ", nlons, ilon_min, ilon_max)
    #print("np_obs_clim_array:", np_obs_clim_array.shape, type(np_obs_clim_array))
    #print("np_fcst_clim_array:", np_fcst_clim_array.shape, type(np_fcst_clim_array))
    logger.info(f"Calling bcsd_function.latlon_calculations", subtask=subtask)
    correct_fcst_coarse = \
        bcsd_function.latlon_calculations(ilat_min, ilat_max, ilon_min, ilon_max,
                                          nlats, nlons, np_obs_clim_array, np_fcst_clim_array,
                                          LEAD_FINAL, TARGET_FCST_SYR, TARGET_FCST_EYR, FCST_SYR,
                                          ENS_NUM, mon, BC_VAR, TINY, fcst_coarse, land_mask)

    correct_fcst_coarse = np.ma.masked_array(correct_fcst_coarse, \
                                             mask=correct_fcst_coarse == -9999.)

    # clip limits - monthly BC NMME precip:
    correct_fcst_coarse = limits.clip_array(correct_fcst_coarse, \
            var_name=CF2VAR.get(FCST_VAR))

    outfile = OUTFILE_TEMPLATE.format(OUTDIR, FCST_VAR, month_name, \
              TARGET_FCST_SYR, TARGET_FCST_EYR)
    logger.info(f"Writing {outfile}", subtask = subtask)
    sdate = datetime(TARGET_FCST_SYR, mon, 1)
    dates = [sdate+relativedelta(years=n) for n in range(correct_fcst_coarse.shape[0])]
    write_4d_netcdf(outfile, correct_fcst_coarse, FCST_VAR, 'CFSv2', \
    'Bias corrected', UNIT, 5, LONS, LATS, ENS_NUM, LEAD_FINAL, sdate, dates)

# Then for forecast files:
for MON in [INIT_FCST_MON]:
    monthly_calculations(MON)
