#!/usr/bin/env python
"""
# Author: Shrad Shukla
# coding: utf-8
#Author: Shrad Shukla
#Usage: This is a module for the BCSD code.
#This module bias corrects a forecasts following probability
#mapping approach as described in Wood et al. 2002
#Date: August 06, 2015
"""
import os
import sys
from datetime import datetime
import calendar
from time import ctime as t_ctime
from time import time as t_time
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
from concurrent.futures import ProcessPoolExecutor
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
# pylint: enable=no-name-in-module
# pylint: disable=import-error
from ghis2s.shared.utils import get_domain_info
from ghis2s.bcsd.bcsd_library.bcsd_function import VarLimits as lim
from ghis2s.shared.logging_utils import TaskLogger
# pylint: enable=import-error

limits = lim()

CF2VAR = {
    'PRECTOT': 'PRECTOT',
    'LWGAB': 'LWS',
    'SWGDN': 'SLRSF',
    'PS': 'PS',
    'QV2M':'Q2M',
    'T2M': 'T2M',
    'U10M': 'WIND',
    'WIND10M': 'WIND',
    }

def scale_forcings (mon_bc_value, mon_raw_value, input_raw_data, bc_var = None):
    ''' perform scaling '''
    output_bc_data = np.ones(len(input_raw_data))*-9999.

    if bc_var == 'PRCP':
        if mon_raw_value == 0.:
            output_bc_data[:] = mon_bc_value
        else:
            correction_factor = mon_bc_value/mon_raw_value
            output_bc_data[:] = input_raw_data[:]*correction_factor
    else:
        correction_factor = mon_bc_value - mon_raw_value
        output_bc_data[:] = input_raw_data[:] + correction_factor

    return output_bc_data

def write_bc_netcdf(outfile, var, varname, description, source, var_units, \
var_standard_name, lons, lats, sdate, dates, sig_digit, north_east_corner_lat, \
north_east_corner_lon, south_west_corner_lat, south_west_corner_lon, \
resolution_x, resolution_y, time_increment):
    """write netcdf"""
    rootgrp = nc4_dataset(outfile, 'w', format='NETCDF4_CLASSIC')
    time = rootgrp.createDimension('time', None)
    longitude = rootgrp.createDimension('lon', len(lons))
    latitude = rootgrp.createDimension('lat', len(lats))

    longitudes = rootgrp.createVariable('lon', 'f4', ('lon',))
    latitudes = rootgrp.createVariable('lat', 'f4', ('lat',))
    times = rootgrp.createVariable('time', 'f4', ('time', ))

    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname, 'f4', ('time', 'lat', \
    'lon',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    rootgrp.missing_value = -9999
    rootgrp.description = description
    rootgrp.zenith_interp = "true,false,"
    rootgrp.MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL"
    rootgrp.conventions = "CF-1.6"
    rootgrp.south_west_corner_lat = float(south_west_corner_lat)
    rootgrp.south_west_corner_lon = float(south_west_corner_lon)
    rootgrp.north_east_corner_lat = float(north_east_corner_lat)
    rootgrp.north_east_corner_lon = float(north_east_corner_lon)
    rootgrp.DX = resolution_x
    rootgrp.DY = resolution_y
    rootgrp.history = 'Created ' + t_ctime(t_time())
    rootgrp.source = source
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    varname.units = var_units
    varname.standard_name = var_standard_name
    string_date = datetime.strftime(sdate, "%Y-%m-%d %H:%M:%S")
    times.units = 'minutes since ' + string_date
    times.time_increment = time_increment
    times.begin_date = datetime.strftime(sdate, "%Y%m%d")
    times.begin_time = '000000'
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    varname[:, :, :] = var
    times[:] = nc4_date2num(dates, units=times.units, calendar=times.calendar)
    rootgrp.close()

task_name = os.environ.get('SCRIPT_NAME')
subtask = f'{sys.argv[1]}'
logger = TaskLogger(task_name,
                    os.getcwd(),
                    f'bcsd/bcsd_library/temporal_disaggregation_6hourly_module.py processing {sys.argv[1]} for month {int(sys.argv[4]):02d}')

## Usage: <Name of variable in observed climatology> <Name of variable in
## reforecast climatology (same as the name in target forecast>
## <forecast model number>

def process_ensemble(MON, ens):
    task_label = subtask + f'-ens{ens:02d}'
    # All file formats
    MONTHLY_BC_INFILE_TEMPLATE = '{}/{}.{}.{}_{:04d}_{:04d}.nc'
    MONTHLY_RAW_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.cfsv2.{:04d}{:02d}.nc'
    SUBDAILY_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.cfsv2.{:04d}{:02d}.nc'
    SUBDAILY_OUTFILE_TEMPLATE = '{}/{}.{:04d}{:02d}.nc4'
    MONTHLY_NMME_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.nmme.monthly.{:04d}{:02d}.nc'
    
    OBS_VAR = str(sys.argv[1]) ##
    FCST_VAR = str(sys.argv[2]) ##
    INIT_FCST_YEAR = int(sys.argv[3])
    ## initial forecast year for which to downscale the data
    INIT_FCST_MON = int(sys.argv[4])
    ## initial forecast month for which to downscale the data
    BC_VAR = str(sys.argv[5])
    ## This is used to figure out if the variable a precipitation variable or not
    UNIT = str(sys.argv[6])
    MODEL_NAME = str(sys.argv[7])
    ENS_NUM = int(sys.argv[8])
    LEAD_FINAL = int(sys.argv[9])
    MONTH_NAME_TEMPLATE = '{}01'
    MONTH_NAME = MONTH_NAME_TEMPLATE.format(calendar.month_abbr[INIT_FCST_MON])
    
    BC_FCST_SYR, BC_FCST_EYR = int(sys.argv[10]), int(sys.argv[11])
    CONFIG_FILE = str(sys.argv[12])
    LAT1, LAT2, LON1, LON2 = get_domain_info(CONFIG_FILE, extent=True)
    LAT_LDT, LON_LDT = get_domain_info(CONFIG_FILE, coord=True)
    
    MONTHLY_BC_FCST_DIR = str(sys.argv[13])
    MONTHLY_RAW_FCST_DIR = str(sys.argv[14])
    SUBDAILY_RAW_FCST_DIR = str(sys.argv[15])
    BASE_OUTDIR = str(sys.argv[16])
    OUTDIR_TEMPLATE = '{}/{:04d}/ens{:01d}'
    
    DOMAIN = str(sys.argv[17])

    MONTH_NAME = MONTH_NAME_TEMPLATE.format((calendar.month_abbr[MON]).lower())
    ## This provides abbrevated version of the name of a month: (e.g. for
    ## January (i.e. Month number = 1) it will return "Jan"). The abbrevated
    ## name is used in the forecasts file name
    
    logger.info(f"Forecast Initialization month is {MONTH_NAME}", subtask=task_label)
    ### First read bias corrected monthly forecast data
    BC_INFILE = MONTHLY_BC_INFILE_TEMPLATE.format(MONTHLY_BC_FCST_DIR,\
                                                  FCST_VAR, MODEL_NAME, MONTH_NAME, BC_FCST_SYR, BC_FCST_EYR)
    
    logger.info(f"Reading bias corrected monthly forecasts {ens} {BC_INFILE}", subtask=task_label)

    MON_BC_DATAG = xr.open_dataset(BC_INFILE)

    LONS = MON_BC_DATAG['longitude'].values
    LATS = MON_BC_DATAG['latitude'].values
    II1 = np.min(np.where (LONS >= LON1))
    II2 = np.max(np.where (LONS <= LON2))
    JJ1 = np.min(np.where (LATS >= LAT1))
    JJ2 = np.max(np.where (LATS <= LAT2))

    MON_BC_DATAG = MON_BC_DATAG.rename({"longitude": "lon", "latitude" : "lat"})
    MON_BC_DATA = MON_BC_DATAG.sel(lon=slice(LON1,LON2),lat=slice(LAT1,LAT2))

    OUTDIR = OUTDIR_TEMPLATE.format(BASE_OUTDIR, INIT_FCST_YEAR, ens+1)
    if os.path.isdir(OUTDIR):
        pass
    else:
        os.makedirs(OUTDIR, exist_ok=True)

    for LEAD_NUM in range(0, LEAD_FINAL): ## Loop from lead =0 to Final Lead
        FCST_DATE = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1, 6) + \
        relativedelta(months=LEAD_NUM)
        FCST_YEAR, FCST_MONTH = FCST_DATE.year, FCST_DATE.month

        # Number of subdaily time steps in the target forecast month
        NUM_TIMESTEPS = 4*calendar.monthrange(FCST_YEAR, FCST_MONTH)[1]

        # Using number of days above to read input daily forecasts
        # and define array to store output file
        OUTFILE = SUBDAILY_OUTFILE_TEMPLATE.format(OUTDIR, OBS_VAR, \
        FCST_YEAR, FCST_MONTH)
        OUTPUT_BC_DATA = np.ones((NUM_TIMESTEPS, len(LATS), len(LONS)))*-9999.
        # Monthly raw data
        if FCST_VAR != 'PRECTOT':
            MONTHLY_INFILE = MONTHLY_RAW_INFILE_TEMPLATE.format(\
            MONTHLY_RAW_FCST_DIR, INIT_FCST_YEAR, ens+1, MONTH_NAME, \
            FCST_YEAR, FCST_MONTH)
        else:
            MONTHLY_INFILE = MONTHLY_NMME_INFILE_TEMPLATE.format(\
            MONTHLY_RAW_FCST_DIR, INIT_FCST_YEAR, ens+1, MONTH_NAME, \
            FCST_YEAR, FCST_MONTH)
        logger.info(f"Reading raw monthly forecast {MONTHLY_INFILE}", subtask=task_label)
        MONTHLY_INPUT_RAW_DATAG = xr.open_dataset(MONTHLY_INFILE)
        MONTHLY_INPUT_RAW_DATA = MONTHLY_INPUT_RAW_DATAG.sel(lon=slice(LON1,LON2),
                                                             lat=slice(LAT1,LAT2))
        MONTHLY_INPUT_RAW_DATA = MONTHLY_INPUT_RAW_DATA[FCST_VAR][:,:]

        # Sub-Daily raw data
        SUBDAILY_INFILE = SUBDAILY_INFILE_TEMPLATE.format(\
        SUBDAILY_RAW_FCST_DIR, INIT_FCST_YEAR, ens+1, MONTH_NAME, \
        FCST_YEAR, FCST_MONTH)
        logger.info(f"Reading raw sub-daily forecast {SUBDAILY_INFILE}", subtask=task_label)
        INPUT_RAW_DATAG = xr.open_dataset(SUBDAILY_INFILE)
        INPUT_RAW_DATA = INPUT_RAW_DATAG.sel(lon=slice(LON1,LON2),lat=slice(LAT1,LAT2))

        # Bias corrected monthly value
        MON_BC_VALUE = MON_BC_DATA[FCST_VAR][INIT_FCST_YEAR-BC_FCST_SYR, LEAD_NUM, ens,:,:]

        # make sure lat/lon are aligned.
        if (not np.array_equal(MONTHLY_INPUT_RAW_DATA["lat"].values,
                               MON_BC_VALUE["lat"].values)) or \
                               (not np.array_equal(MONTHLY_INPUT_RAW_DATA["lon"].values,
                                                   MON_BC_VALUE["lon"].values)):
            MONTHLY_INPUT_RAW_DATA({"lon": MON_BC_VALUE["lon"].values,
                                    "lat": MON_BC_VALUE["lat"].values})
        if (not np.array_equal(INPUT_RAW_DATA["lat"].values,
                               MON_BC_VALUE["lat"].values)) or \
                               (not np.array_equal(INPUT_RAW_DATA["lon"].values,
                                                   MON_BC_VALUE["lon"].values)):
            INPUT_RAW_DATA({"lon": MON_BC_VALUE["lon"].values,
                            "lat": MON_BC_VALUE["lat"].values})

        correct = xr.apply_ufunc(
            scale_forcings,
            MON_BC_VALUE.chunk({"lat": "auto", "lon": "auto"}).compute(),
            MONTHLY_INPUT_RAW_DATA.chunk({"lat": "auto", "lon": "auto"}).compute(),
            INPUT_RAW_DATA[FCST_VAR].chunk({"lat": "auto", "lon": "auto"}).compute(),
            input_core_dims=[[],[],['time']],
            exclude_dims=set(('time',)),
            output_core_dims=[['time']],
            vectorize=True,
            dask="forbidden",
            output_dtypes=[np.float64],
            kwargs={'bc_var': BC_VAR},
        )

        correct2 = np.moveaxis(correct.values,2,0)
        OUTPUT_BC_DATA[:,JJ1:JJ2+1, II1:II2+1] = correct2[:,:,:]

        ### Finish correcting values for all timesteps in the given
        ### month and ensemble member
        # clip limits
        if OBS_VAR == 'PRECTOT':
            OUTPUT_BC_DATA = limits.clip_array(OUTPUT_BC_DATA, var_name=CF2VAR.get(OBS_VAR), precip=True)
        else:
            OUTPUT_BC_DATA = limits.clip_array(OUTPUT_BC_DATA, var_name=CF2VAR.get(OBS_VAR))

        logger.info(f"Writing {OUTFILE}", subtask=task_label)
        OUTPUT_BC_DATA = np.ma.masked_array(OUTPUT_BC_DATA, \
                                            mask=OUTPUT_BC_DATA == -9999.)
        date = [FCST_DATE+relativedelta(hours=n*6) for n in \
        range(NUM_TIMESTEPS)]

        write_bc_netcdf(OUTFILE, OUTPUT_BC_DATA, OBS_VAR, \
                        'Bias corrected forecasts', 'MODEL:'  +   MODEL_NAME, UNIT, \
                        OBS_VAR, LONS, LATS, FCST_DATE, date, 8, np.max(LAT_LDT), np.max(LON_LDT), \
                        np.min(LAT_LDT), np.min(LON_LDT), LON_LDT[1] - LON_LDT[0], \
                        LAT_LDT[1] - LAT_LDT[0], 21600)
        
logger.info("Starting parallel processing of ensemmbles")
for MON in [int(sys.argv[4])]:
    num_workers = int(os.environ.get('NUM_WORKERS', int(sys.argv[8])))
    # ProcessPoolExecutor parallel processing
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for ens in range(int(sys.argv[8])):
            logger.info(f"Submitting disaggregation job for ens {ens:02d}", subtask=subtask + f'-ens{ens:02d}')
            future = executor.submit(process_ensemble, MON, ens)
            futures.append(future)
    
        for future in futures:
            result = future.result()
            
    logger.info(f"CFSv2 6-hourly disaggregation completed successfully")
