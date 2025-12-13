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
import calendar
import gc
import subprocess
from time import ctime as t_ctime
from time import time as t_time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
from numpy import ma
import xarray as xr
import yaml
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
# pylint: enable=no-name-in-module
from ghis2s.bcsd.bcsd_library.bcsd_functions import VarLimits as lim
from ghis2s.shared.utils import get_domain_info, load_ncdata, get_chunk_sizes
from ghis2s.shared.logging_utils import TaskLogger

limits = lim()
PRECIP_THRES = limits.precip_thres

## Usage: <Name of variable in observed climatology> <Name of variable in
## reforecast climatology (same as the name in target forecast>
## <forecast model number>
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
MONTH_NAME = MONTH_NAME_TEMPLATE.format(calendar.month_abbr[INIT_FCST_MON].lower())

BC_FCST_SYR, BC_FCST_EYR = int(sys.argv[10]), int(sys.argv[11])
CONFIG_FILE = str(sys.argv[12])
LAT1, LAT2, LON1, LON2 = get_domain_info(CONFIG_FILE, extent=True)
LAT_LDT, LON_LDT = get_domain_info(CONFIG_FILE, coord=True)
RESOL = f'{round((LAT_LDT[1] - LAT_LDT[0])*100)}km'
with open(CONFIG_FILE, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)
FCST_MODEL = config['BCSD']['metforce_source']
MONTHLY_BC_FCST_DIR = str(sys.argv[13])
SUBDAILY_RAW_FCST_DIR = str(sys.argv[14])
BASE_OUTDIR = str(sys.argv[15])

if len(LAT_LDT) == 1800 and len(LON_LDT) == 3600:
    # Medium resolution: ~0.1° grid
    LAT_CHUNK, LON_CHUNK = get_chunk_sizes(None, dim_in=[1800, 3600])
elif len(LAT_LDT) == 3600 and len(LON_LDT) == 7200:
    # High resolution: ~0.05° grid
    LAT_CHUNK, LON_CHUNK = get_chunk_sizes(None, dim_in=[3600, 7200])
else:
    LAT_CHUNK, LON_CHUNK = get_chunk_sizes(None, dim_in=[len(LAT_LDT), len(LON_LDT)])

# All file templates
MONTHLY_BC_INFILE_TEMPLATE = '{}/{}.{}.{}_{:04d}{:02d}.nc'
SUBDAILY_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.{}.{:04d}{:02d}.nc'
SUBDAILY_OUTFILE_TEMPLATE = '{}/{}.{:04d}{:02d}.nc4'
OUTDIR_TEMPLATE = '{}/{:04d}/ens{:01d}'

task_name = os.environ.get('SCRIPT_NAME')
SUBTASK = f'{sys.argv[7]}'
logger = TaskLogger(task_name,
                    os.getcwd(),
                    'bcsd/bcsd_library/temporal_disaggregation_nmme_6hourly_module.py processing '
                    f'{MODEL_NAME} for month {INIT_FCST_MON:02d}')

def use_neighbors_diurnal_cycle (precip_data):
    ''' 
    a function to get sub-monthly distribution from a neighboring grid cell when NNMME has 
    precip and CFSv2 doesn't
    '''
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

    return output_bc_data

def write_bc_netcdf(outfile, var, obs_var, description, source, var_units, \
var_standard_name, lons, lats, sdate, dates, sig_digit, north_east_corner_lat, \
north_east_corner_lon, south_west_corner_lat, south_west_corner_lon, \
resolution_x, resolution_y, time_increment):
    """write netcdf"""
    n_times = len(dates)
    complevel = 6
    # Determine optimal chunking based on grid size
    if len(lats) == 1800 and len(lons) == 3600:
        # Medium resolution: ~0.1° grid
        time_chunk = 24

    elif len(lats) == 3600 and len(lons) == 7200:
        # High resolution: ~0.05° grid
        time_chunk = 12

    else:
        # Fallback for other grid sizes
        time_chunk = min(24, n_times)

    rootgrp = nc4_dataset(outfile, 'w', format='NETCDF4_CLASSIC')
    time = rootgrp.createDimension('time', None)
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))

    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    times = rootgrp.createVariable('time', 'f4', ('time', ))

    # two dimensions unlimited.
    varname = rootgrp.createVariable(obs_var, 'f4', ('time', 'latitude', 'longitude',),
                                     fill_value=-9999, zlib=True,
                                     complevel=complevel,
                                     shuffle=True,
                                     chunksizes=(time_chunk, LAT_CHUNK, LON_CHUNK),
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
    times[:] = nc4_date2num(dates, units=times.units, calendar=times.calendar)
    for i in range(0, n_times, time_chunk):
        end_idx = min(i + time_chunk, n_times)
        varname[i:end_idx, :, :] = var[i:end_idx, :, :]

    rootgrp.close()
    del rootgrp, varname

FORCE_DT = 21600
if FCST_MODEL == 'CFSv2':
    FORCE_DT = 21600
if FCST_MODEL == 'GEOSv3':
    FORCE_DT = 10800

def process_ensemble(ens):
    ''' process ensemble members '''
    task_label = f'{SUBTASK}-ens{ens+1:02d}'
    logger.info(f"Forecast Initialization month is {MONTH_NAME}", subtask=task_label)

    ## Shape of the above dataset time, Lead, Ens, latitude, longitude
    outdir = OUTDIR_TEMPLATE.format(BASE_OUTDIR, INIT_FCST_YEAR, ens+1)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    for lead_num in range(0, LEAD_FINAL): ## Loop from lead =0 to Final Lead
        fcst_date = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1, 6) + \
        relativedelta(months=lead_num)
        fcst_year, fcst_month = fcst_date.year, fcst_date.month

        # First read bias corrected monthly forecast data
        bc_infile = MONTHLY_BC_INFILE_TEMPLATE.format(MONTHLY_BC_FCST_DIR,\
                                                      FCST_VAR, MODEL_NAME, MONTH_NAME,
                                                      fcst_year, fcst_month)

        logger.info(f"Reading bias corrected monthly forecasts {bc_infile}", subtask=task_label)
        mon_bc_datag = load_ncdata(bc_infile, [logger, task_label])
        lons = mon_bc_datag['longitude'].values
        lats = mon_bc_datag['latitude'].values
        ii1 = np.min(np.where (lons >= LON1))
        ii2 = np.max(np.where (lons <= LON2))
        jj1 = np.min(np.where (lats >= LAT1))
        jj2 = np.max(np.where (lats <= LAT2))

        mon_bc_datag = mon_bc_datag.rename({"longitude": "lon", "latitude" : "lat"})
        mon_bc_data = mon_bc_datag[FCST_VAR].sel(lon=slice(LON1,LON2),lat=slice(LAT1,LAT2))
        mon_bc_datag.close()
        del mon_bc_datag

        # Number of subdaily time steps in the target forecast month
        num_timesteps = (86400//FORCE_DT)*calendar.monthrange(fcst_year, fcst_month)[1]

        # Using number of days above to read input daily forecasts
        # and define array to store output file
        outfile = SUBDAILY_OUTFILE_TEMPLATE.format(outdir, OBS_VAR, \
                                                   fcst_year, fcst_month)
        output_bc_data = np.ones((num_timesteps, len(lats), len(lons)))*-9999.

        # Sub-Daily raw data
        subdaily_infile = SUBDAILY_INFILE_TEMPLATE.format(SUBDAILY_RAW_FCST_DIR,
                                                          INIT_FCST_YEAR, (ens % 12) + 1,
                                                          MONTH_NAME, FCST_MODEL.lower(),
                                                          fcst_year, fcst_month)
        logger.info(f"Reading raw sub-daily forecast {subdaily_infile}", subtask=task_label)

        monthly_input_raw_datag = load_ncdata(subdaily_infile, [logger, task_label],
                                              chunks={'lat': "auto", 'lon': "auto"})
        input_raw_data = monthly_input_raw_datag.sel(lon=slice(LON1,LON2),lat=slice(LAT1,LAT2))
        input_raw_data = input_raw_data.chunk({"time": -1, "lat": LAT_CHUNK, "lon": LON_CHUNK})
        monthly_input_raw_data = input_raw_data[FCST_VAR].mean(dim = 'time')

        # Bias corrected monthly value
        mon_bc_value = mon_bc_data[ens , :,:]

        # make sure lat/lon are aligned.
        if (not np.array_equal(monthly_input_raw_data["lat"].values,
                               mon_bc_value["lat"].values)) or \
                               (not np.array_equal(monthly_input_raw_data["lon"].values,
                                                       mon_bc_value["lon"].values)):
            monthly_input_raw_data = monthly_input_raw_data.assign_coords \
                ({"lon": mon_bc_value["lon"].values,
                  "lat": mon_bc_value["lat"].values})

        if (not np.array_equal(input_raw_data["lat"].values,
                               mon_bc_value["lat"].values)) or \
                               (not np.array_equal(input_raw_data["lon"].values,
                                                       mon_bc_value["lon"].values)):

            input_raw_data = input_raw_data.assign_coords({"lon": mon_bc_value["lon"].values,
                                                           "lat": mon_bc_value["lat"].values})

        correct = xr.apply_ufunc(
            scale_forcings,
            mon_bc_value.chunk({"lat": LAT_CHUNK, "lon": LON_CHUNK}),
            monthly_input_raw_data.chunk({"lat": LAT_CHUNK, "lon": LON_CHUNK}),
            input_raw_data[FCST_VAR].chunk({"lat": LAT_CHUNK, "lon": LON_CHUNK}),
            input_core_dims=[[],[],['time']],
            exclude_dims=set(('time',)),
            output_core_dims=[['time']],
            dask_gufunc_kwargs={"output_sizes": {'time': num_timesteps}},
            vectorize=True,
            dask="parallelized",
            output_dtypes=[np.float32],
            kwargs={'bc_var': BC_VAR},
        )

        output_bc_data[:, jj1:jj2+1, ii1:ii2+1] = correct.transpose('time', 'lat', 'lon').values
        logger.info("Calling use_neighbors_diurnal_cycle", subtask=task_label)
        # Find neighboring OUTPUT_BC_DATA to add sub-monthly distribution
        output_bc_revised, cnt_beg, cnt_end = use_neighbors_diurnal_cycle (output_bc_data)
        logger.info(f"NOF cells without precip diurnal cycle: {cnt_beg} (before) {cnt_end} (after)",
                    subtask=task_label)

        # clip limits - 6hr NMME bcsd files:
        output_bc_revised = limits.clip_array(output_bc_revised, var_name="PRECTOT", precip=True)

        ### Finish correcting values for all timesteps in the given
        ### month and ensemble member
        logger.info(f"Writing {outfile}", subtask=task_label)
        output_bc_revised = np.ma.masked_array(output_bc_revised, \
                                               mask=output_bc_revised == -9999.)
        date = [fcst_date+relativedelta(hours=n*6)
                for n in range(num_timesteps)]

        write_bc_netcdf(outfile, output_bc_revised, OBS_VAR,
                        'Bias corrected forecasts', 'MODEL:'  +   MODEL_NAME, UNIT,
                        OBS_VAR, lons, lats, fcst_date, date, 8, np.max(LAT_LDT),
                        np.max(LON_LDT), np.min(LAT_LDT), np.min(LON_LDT), LON_LDT[1] - LON_LDT[0],
                        LAT_LDT[1] - LAT_LDT[0], FORCE_DT)

        del output_bc_revised, output_bc_data, correct
        monthly_input_raw_datag.close()
        input_raw_data.close()
        monthly_input_raw_data.close()
        del monthly_input_raw_datag, input_raw_data, monthly_input_raw_data
        gc.collect()

def links():
    ''' create links for lead_months + 1'''
    init_datetime = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1)

    src_yyyymm = []
    dst_yyyymm = []
    for mon in range(LEAD_FINAL):
        src_yyyymm.append((init_datetime + relativedelta(months=mon)).strftime("%Y%m"))
        dst_yyyymm.append((init_datetime + relativedelta(months=mon)).strftime("%Y%m"))

    src_yyyymm.append((init_datetime + relativedelta(months=mon)).strftime("%Y%m"))
    dst_yyyymm.append((init_datetime + relativedelta(months=mon+1)).strftime("%Y%m"))
    last_yyyymm = len(src_yyyymm) -1

    logger.info(f"Creating symbolic links for month {LEAD_FINAL +1}")

    for iens, _ in enumerate(range(ENS_NUM)):
        outdir_ = OUTDIR_TEMPLATE.format(BASE_OUTDIR, INIT_FCST_YEAR, iens + 1)
        src_file = f"{outdir_}/{OBS_VAR}.{src_yyyymm[last_yyyymm]}.nc4"
        dst_file = f"{outdir_}/{OBS_VAR}.{dst_yyyymm[last_yyyymm]}.nc4"
        cmd = f"ln -sfn {src_file} {dst_file}"
        rc = subprocess.call(cmd, shell=True)
        if rc != 0:
            logger.error(f"Problem calling creating last precip symbolic link to {dst_file}!")

logger.info("Starting parallel processing of number of ensemble members")
num_workers = int(os.environ.get('NUM_WORKERS', ENS_NUM))

if RESOL != '25km':
    start_ens = int(sys.argv[16])
    end_ens = int(sys.argv[17])
    num_workers = end_ens - start_ens + 1  # Number of workers = number of ensembles in this job
    logger.info(f"Standard resolution grid ({RESOL}): "
                f"processing ensembles {start_ens}-{end_ens} with {num_workers} workers")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for _ens in range(start_ens, end_ens + 1):
            logger.info(f"Submitting disaggregation job for ens {_ens+1:02d}",
                        subtask=SUBTASK + f'ens{_ens+1:02d}')
            future = executor.submit(process_ensemble, _ens)
            futures.append(future)

        for future in futures:
            result = future.result()

else:
    # Standard resolution: process all at once
    logger.info(f"Standard resolution grid ({len(LAT_LDT)}x{len(LON_LDT)}): "
                f"processing all {ENS_NUM} ensembles with {num_workers} workers")

    with ProcessPoolExecutor(max_workers=num_workers) as _executor:
        _futures = []
        for _ens in range(ENS_NUM):
            logger.info(f"Submitting disaggregation job for ens {_ens+1:02d}",
                        subtask=SUBTASK + f'-ens{_ens+1:02d}')
            _future = _executor.submit(process_ensemble, _ens)
            _futures.append(_future)

        for _future in _futures:
            result_ = _future.result()

        links()

logger.info("RAN SUCCESSFULLY !")
