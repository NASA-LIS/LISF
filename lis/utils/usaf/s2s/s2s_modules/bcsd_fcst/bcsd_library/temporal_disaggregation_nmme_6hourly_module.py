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
from numpy import ma
import xarray as xr
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
# pylint: enable=no-name-in-module
# pylint: disable=import-error
from bcsd_function import VarLimits as lim
from bcsd_stats_functions import get_domain_info
# pylint: enable=import-error

limits = lim()
PRECIP_THRES = limits.precip_thres

def use_neighbors_diurnal_cycle (precip_data):
    ''' a function to get sub-monthly distribution from a neighboring grid cell when NNMME has precip and CFSv2 doesn't'''
    global PRECIP_THRES

    # make a copy, update, and return
    arrayb = np.copy(precip_data)

    # total precipitation across time axis
    total_precip = np.mean(precip_data, axis=0)

    # Find grid cells with total precipitation greater than PRECIP_THRES
    positive_total_mask = total_precip >= PRECIP_THRES

    # Find grid cells with uniform precipitation values along the time axis
    uniform_mask = np.all(precip_data == precip_data[0,:], axis=0)

    # When both maks are TRUE, precip is positive and uniform
    target_cells_mask = uniform_mask & positive_total_mask

    # non-uniform and positive mask is where useful 6-hourly distributions are available
    useful_cells_mask = ~uniform_mask & positive_total_mask

    # mask out the cells that need to be massaged and where 6-hourly distributions are available
    current_precip = ma.masked_array(total_precip, target_cells_mask)
    useful_precip = ma.masked_array(total_precip, ~useful_cells_mask)
    tgt_cells_beg = np.sum(target_cells_mask)

    for zoom in range (1,3):
        # starting from the 3x3 window from the location gradually increase upto 39x39 window
        for direction in (-1,1):
            shift = direction * zoom

            # North-South direction [|]
            axis = 0
            if not np.any(current_precip.mask): break
            a_shifted = np.roll(useful_precip, shift=shift, axis=axis)
            arrayb_shifted = np.roll(arrayb ,shift=shift,axis=axis+1)
            shifted_ratio = np.mean(arrayb, axis=0)/a_shifted
            idx=~a_shifted.mask * current_precip.mask
            # update return array and target_cells_mask
            for t in range (0,arrayb.shape[0]):
                arrayb[t,idx]=arrayb_shifted[t,idx]*shifted_ratio[idx]
            target_cells_mask = target_cells_mask & ~idx
            current_precip = ma.masked_array(total_precip, target_cells_mask)

            #  East-west direction [-]
            axis = 1
            if not np.any(current_precip.mask): break
            a_shifted = np.roll(useful_precip, shift=shift, axis=axis)
            arrayb_shifted = np.roll(arrayb ,shift=shift,axis=axis+1)
            shifted_ratio = np.mean(arrayb, axis=0)/a_shifted
            idx=~a_shifted.mask * current_precip.mask
            # update return array and target_cells_mask
            for t in range (0,arrayb.shape[0]):
                arrayb[t,idx]=arrayb_shifted[t,idx]*shifted_ratio[idx]

            target_cells_mask = target_cells_mask & ~idx
            current_precip = ma.masked_array(total_precip, target_cells_mask)

            # Diagonal [\]
            axis = (0,1)
            if not np.any(current_precip.mask): break
            a_shifted = np.roll(useful_precip, shift=shift, axis=axis)
            arrayb_shifted = np.roll(arrayb ,shift=shift,axis=(1,2))
            shifted_ratio = np.mean(arrayb, axis=0)/a_shifted
            idx=~a_shifted.mask * current_precip.mask
            # update return array and target_cells_mask
            for t in range (0,arrayb.shape[0]):
                arrayb[t,idx]=arrayb_shifted[t,idx]*shifted_ratio[idx]
            target_cells_mask = target_cells_mask & ~idx
            current_precip = ma.masked_array(total_precip, target_cells_mask)

            # Diagonal [/]
            axis = (0,1)
            if not np.any(current_precip.mask): break
            useful_precip_flip = np.flip(useful_precip,axis=1)
            arrayb_flip = np.flip(arrayb, axis=2)
            a_shifted_flip = np.roll(useful_precip_flip, shift=shift, axis=axis)
            arrayb_shifted_flip = np.roll(arrayb_flip ,shift=shift,axis=(1,2))
            shifted_ratio_flip = np.mean(arrayb_flip, axis=0)/a_shifted_flip
            #unflip
            arrayb_shifted = np.flip(arrayb_shifted_flip, axis=2)
            shifted_ratio = np.flip(shifted_ratio_flip, axis=1)
            a_shifted = np.flip(a_shifted_flip, axis=1)
            idx=~a_shifted.mask * current_precip.mask
            # update return array and target_cells_mask
            for t in range (0,arrayb.shape[0]):
                arrayb[t,idx]=arrayb_shifted[t,idx]*shifted_ratio[idx]
            target_cells_mask = target_cells_mask & ~idx
            current_precip = ma.masked_array(total_precip, target_cells_mask)

    for t in range (0,arrayb.shape[0]):
        arrayb[t,target_cells_mask]= 0.

    tgt_cells_end = np.sum(target_cells_mask)
    return arrayb, tgt_cells_beg, tgt_cells_end

def scale_forcings (mon_bc_value, mon_raw_value, input_raw_data, bc_var = None):
    ''' perform scaling '''
    output_bc_data = np.ones(len(input_raw_data))*-9999.

    if bc_var == 'PRCP':
        if mon_bc_value == -9999:
            return output_bc_data
        if mon_raw_value == 0.:
            output_bc_data[:] = mon_bc_value
        else:
            correction_factor = mon_bc_value/mon_raw_value
            output_bc_data[:] = input_raw_data[:]*correction_factor
#    else:
#        correction_factor = mon_bc_value - mon_raw_value
#        output_bc_data[:] = input_raw_data[:] + correction_factor

    return output_bc_data

def write_bc_netcdf(outfile, var, obs_var, description, source, var_units, \
var_standard_name, lons, lats, sdate, dates, sig_digit, north_east_corner_lat, \
north_east_corner_lon, south_west_corner_lat, south_west_corner_lon, \
resolution_x, resolution_y, time_increment):
    """write netcdf"""
    rootgrp = nc4_dataset(outfile, 'w', format='NETCDF4_CLASSIC')
    time = rootgrp.createDimension('time', None)
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))

    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    times = rootgrp.createVariable('time', 'f4', ('time', ))

    # two dimensions unlimited.
    varname = rootgrp.createVariable(obs_var, 'f4', ('time', 'latitude', \
    'longitude',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    rootgrp.missing_value = -9999.
    rootgrp.description = description
    rootgrp.zenith_interp = "true,false,"
    rootgrp.MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL"
    rootgrp.conventions = "CF-1.6"
    rootgrp.SOUTH_WEST_CORNER_LAT = float(south_west_corner_lat)
    rootgrp.SOUTH_WEST_CORNER_LON = float(south_west_corner_lon)
    rootgrp.NORTH_EAST_CORNER_LAT = float(north_east_corner_lat)
    rootgrp.NORTH_EAST_CORNER_LON = float(north_east_corner_lon)
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

## Usage: <Name of variable in observed climatology> <Name of variable in
## reforecast climatology (same as the name in target forecast>
## <forecast model number>
CMDARGS = str(sys.argv)
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

print(f"*** LEAD FINAL: {LEAD_FINAL}")
BC_FCST_SYR, BC_FCST_EYR = int(sys.argv[10]), int(sys.argv[11])
CONFIG_FILE = str(sys.argv[12])
LAT1, LAT2, LON1, LON2 = get_domain_info(CONFIG_FILE, extent=True)

MONTHLY_BC_FCST_DIR = str(sys.argv[13])
SUBDAILY_RAW_FCST_DIR = str(sys.argv[14])
BASE_OUTDIR = str(sys.argv[15])
OUTDIR_TEMPLATE = '{}/{:04d}/ens{:01d}'
DOMAIN=str(sys.argv[16])

# All file formats
MONTHLY_BC_INFILE_TEMPLATE = '{}/{}.{}.{}_{:04d}_{:04d}.nc'
SUBDAILY_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.cfsv2.{:04d}{:02d}.nc'
SUBDAILY_OUTFILE_TEMPLATE = '{}/{}.{:04d}{:02d}.nc4'

for MON in [INIT_FCST_MON]:
    MONTH_NAME = MONTH_NAME_TEMPLATE.format((calendar.month_abbr[MON]).lower())
    ## This provides abbrevated version of the name of a month: (e.g. for
    ## January (i.e. Month number = 1) it will return "Jan"). The abbrevated
    ## name is used in the forecasts file name
    print(f"Forecast Initialization month is {MONTH_NAME}")
    ### First read bias corrected monthly forecast data
    BC_INFILE = MONTHLY_BC_INFILE_TEMPLATE.format(MONTHLY_BC_FCST_DIR,\
    FCST_VAR, MODEL_NAME, MONTH_NAME, BC_FCST_SYR, BC_FCST_EYR)

    print(f"Reading bias corrected monthly forecasts {BC_INFILE}")
    MON_BC_DATAG = xr.open_dataset(BC_INFILE)
    LONS = MON_BC_DATAG['longitude'].values
    LATS = MON_BC_DATAG['latitude'].values
    II1 = np.min(np.where (LONS >= LON1))
    II2 = np.max(np.where (LONS <= LON2))
    JJ1 = np.min(np.where (LATS >= LAT1))
    JJ2 = np.max(np.where (LATS <= LAT2))

    MON_BC_DATAG = MON_BC_DATAG.rename({"longitude": "lon", "latitude" : "lat"})
    MON_BC_DATA = MON_BC_DATAG.sel(lon=slice(LON1,LON2),lat=slice(LAT1,LAT2))

    ## Shape of the above dataset time, Lead, Ens, latitude, longitude
    for ens in range(ENS_NUM):
        OUTDIR = OUTDIR_TEMPLATE.format(BASE_OUTDIR, INIT_FCST_YEAR, ens+1)
        if os.path.isdir(OUTDIR):
            pass
        else:
            os.makedirs(OUTDIR, exist_ok=True)
        print(f"OUTDIR is {OUTDIR}")
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

            # Sub-Daily raw data
            SUBDAILY_INFILE = SUBDAILY_INFILE_TEMPLATE.format(\
            SUBDAILY_RAW_FCST_DIR, INIT_FCST_YEAR, ens+1, MONTH_NAME, \
            FCST_YEAR, FCST_MONTH)
            print(f"Reading raw sub-daily forecast {SUBDAILY_INFILE}")
            MONTHLY_INPUT_RAW_DATAG = xr.open_dataset(SUBDAILY_INFILE)
            INPUT_RAW_DATA = MONTHLY_INPUT_RAW_DATAG.sel(lon=slice(LON1,LON2),lat=slice(LAT1,LAT2))
            MONTHLY_INPUT_RAW_DATA = INPUT_RAW_DATA[FCST_VAR].mean(dim = 'time')

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

            # Find neighboring OUTPUT_BC_DATA to add sub-monthly distribution
            OUTPUT_BC_REVISED, cnt_beg, cnt_end = use_neighbors_diurnal_cycle (OUTPUT_BC_DATA)
            print (f'NOF cells without precip diurnal cycle : {cnt_beg} (before) {cnt_end} (after)')

            # clip limits - 6hr NMME bcsd files:
            OUTPUT_BC_REVISED = limits.clip_array(OUTPUT_BC_REVISED, var_name="PRECTOT", precip=True)

            ### Finish correcting values for all timesteps in the given
            ### month and ensemble member
            print(f"Now writing {OUTFILE}")
            OUTPUT_BC_REVISED = np.ma.masked_array(OUTPUT_BC_REVISED, \
                                                   mask=OUTPUT_BC_REVISED == -9999.)
            date = [FCST_DATE+relativedelta(hours=n*6) for n in \
            range(NUM_TIMESTEPS)]

            if DOMAIN == 'AFRICOM':
                write_bc_netcdf(OUTFILE, OUTPUT_BC_REVISED, OBS_VAR, \
                                'Bias corrected forecasts', 'MODEL:'  +   MODEL_NAME, UNIT, \
                                OBS_VAR, LONS, LATS, FCST_DATE, date, 8, 39.875, 59.875, -39.875, \
                                -19.875, 0.25, 0.25, 21600)
            if DOMAIN == 'GLOBAL':
                write_bc_netcdf(OUTFILE, OUTPUT_BC_REVISED, OBS_VAR, \
                                'Bias corrected forecasts', 'MODEL:'  +   MODEL_NAME, UNIT, \
                                OBS_VAR, LONS, LATS, FCST_DATE, date, 8, 89.875, 179.875, -89.875, \
                                -179.875, 0.25, 0.25, 21600)
