#!/usr/bin/env python3
#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
#------------------------------------------------------------------------------
#
# SCRIPT: compute_weekly_anom.py
#
# PURPOSE: Calculates weekly anomaly metrics to be written out.
#
# REVISION HISTORY:
# 28 Jul 2025: S. Mahanama, initial code
#
#------------------------------------------------------------------------------
"""

from datetime import date
import time
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dateutil.relativedelta import relativedelta
import yaml
import numpy as np
import xarray as xr
from ghis2s.shared.logging_utils import TaskLogger
from ghis2s.shared.utils import write_ncfile, load_ncdata
from ghis2s.s2smetric.metricslib import (sel_var, compute_anomaly, compute_sanomaly,
                                         merged_metric_filename, LONG_NAMES_ANOM,
                                         LONG_NAMES_SANOM, UNITS_ANOM, UNITS_SANOM)
# pylint: disable=too-many-arguments,too-many-locals

LEAD_WEEKS = 6
FCST_INIT_MON = int(sys.argv[1])
TARGET_YEAR = int(sys.argv[2])
NMME_MODEL = sys.argv[3]
CONFIGFILE = sys.argv[4]
BASEOUTDIR = sys.argv[5]

# Load CONFIG file
with open(CONFIGFILE, 'r', encoding="utf-8") as file:
    CONFIG = yaml.safe_load(file)

HYD_MODEL = CONFIG["EXP"]["lsmdir"]
DOMAIN_NAME = CONFIG["EXP"]["DOMAIN"]
CLIM_SYR = int(CONFIG["BCSD"]["clim_start_year"])
CLIM_EYR = int(CONFIG["BCSD"]["clim_end_year"])
HINDCASTS = CONFIG["SETUP"]["E2ESDIR"] + '/hindcast/s2spost/' + f'{FCST_INIT_MON:02d}/'
FORECASTS = "./s2spost/"
CURRENTDATE = date(TARGET_YEAR, FCST_INIT_MON, 2)
CLIM_STATFILE_TEMPLATE = \
    '{}/{}_PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.{}_' \
    'PA.ALL_DD.YYYY{:02d}01_DT.0000_FP.YYYY{:02d}{:02d}-YYYY{:02d}{:02d}_DF.NC'

TARGET_INFILE_TEMPLATE = \
    '{}/{:04d}{:02d}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.{}_' \
    'PA.ALL_DD.{:04d}{:02d}01_DT.0000_FP.{:04d}{:02d}{:02d}-{:04d}{:02d}{:02d}_DF.NC'
ENDDATE = CURRENTDATE
for _ in range(LEAD_WEEKS):
    ENDDATE += relativedelta(days=7)
ENDDATE -= relativedelta(days=1)
comp = {'zlib':True, 'complevel':6, 'shuffle':True, 'missing_value': -9999., '_FillValue': -9999.}

if CONFIG['SETUP']['DATATYPE'] == 'hindcast':
    CLIM_INFILE_TEMPLATE = \
        '{}/????{:02d}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.{}_' \
        'PA.ALL_DD.*{:02d}01_DT.0000_FP.*{:02d}{:02d}-*{:02d}{:02d}_DF.NC'

    num_vars = len(CONFIG["EXP"]["NMME_models"])
    num_workers = int(os.environ.get('NUM_WORKERS', num_vars))
    def process_model(nmme_model):
        """
        Process each NMME model for weekly outputs.
        """
        curdate = CURRENTDATE
        enddate = curdate
        enddate += relativedelta(days=6)
        mean_file = CLIM_STATFILE_TEMPLATE.format(HINDCASTS, "MEAN", nmme_model.upper(),
                                                 DOMAIN_NAME, FCST_INIT_MON,
                                                 CURRENTDATE.month, CURRENTDATE.day,
                                                 ENDDATE.month, ENDDATE.day)
        std_file = CLIM_STATFILE_TEMPLATE.format(HINDCASTS, "STD", nmme_model.upper(),
                                                 DOMAIN_NAME, FCST_INIT_MON,
                                                 CURRENTDATE.month, CURRENTDATE.day,
                                                 ENDDATE.month, ENDDATE.day)
        mean_xr = []
        std_xr = []
        for _ in range(LEAD_WEEKS):
            infile = CLIM_INFILE_TEMPLATE.format(HINDCASTS, FCST_INIT_MON,
                                                 nmme_model,nmme_model.upper(),
                                                 DOMAIN_NAME, FCST_INIT_MON,
                                                 curdate.month, curdate.day,
                                                 enddate.month, enddate.day)
            curdate += relativedelta(days=7)
            enddate += relativedelta(days=7)
            week_xr = xr.open_mfdataset(infile, combine='by_coords')
            mean_xr.append(week_xr.mean(dim = ['time','ensemble']))
            std_xr.append(week_xr.std(dim = ['time','ensemble']))

        # Concatenate the lists and add lead time dimension
        mean_concat = xr.concat(mean_xr, dim='lead_time')
        std_concat = xr.concat(std_xr, dim='lead_time')

        # Add lead time coordinates (weeks 1, 2, 3, etc.)
        lead_time_coords = range(1, LEAD_WEEKS + 1)
        mean_concat = mean_concat.assign_coords(lead_time=lead_time_coords)
        std_concat = std_concat.assign_coords(lead_time=lead_time_coords)

        # Add attributes
        mean_concat.lead_time.attrs['units'] = 'weeks'
        mean_concat.lead_time.attrs['long_name'] = 'Lead time in weeks'
        std_concat.lead_time.attrs['units'] = 'weeks'
        std_concat.lead_time.attrs['long_name'] = 'Lead time in weeks'

        # save stat files
        encoding = {var: comp for var in mean_concat.data_vars}
        mean_concat.to_netcdf(mean_file, format="NETCDF4", encoding=encoding)
        encoding = {var: comp for var in std_concat.data_vars}
        std_concat.to_netcdf(std_file, format="NETCDF4", encoding=encoding)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for nmme_model in CONFIG["EXP"]["NMME_models"]:
            future = executor.submit(process_model, nmme_model)
            futures.append(future)

        for future in futures:
            result = future.result()
    sys.exit()

WEEKLY_VARS = CONFIG["POST"]["weekly_vars"]

OUTDIR = BASEOUTDIR + '/'
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR, exist_ok=True)

OUTFILE = merged_metric_filename(OUTDIR, CURRENTDATE, ENDDATE, NMME_MODEL, DOMAIN_NAME, weekly=True)
task_name = os.environ.get('SCRIPT_NAME')
logger = TaskLogger(task_name,
                    os.getcwd(),
                    f'{NMME_MODEL} running s2smetric/compute_weekly_anom.py')

num_vars = 2*len(WEEKLY_VARS)
num_workers = int(os.environ.get('NUM_WORKERS', num_vars))

def process_variable(var_name, anom):
    """
    Function to process each variable and metric (e.g., anomalies)
    """
    logger.info(f"Starting processing {var_name} for metric {anom}", subtask=var_name)
    if anom == "ANOM":
        long_names = LONG_NAMES_ANOM
        units = UNITS_ANOM
    else:
        long_names = LONG_NAMES_SANOM
        units = UNITS_SANOM

    # Create array to store anomaly data (initialized here for each variable)
    all_anom = None
    curdate = CURRENTDATE
    enddate = curdate
    enddate += relativedelta(days=6)
    logger.info("Reading Hindcast statistics", subtask=var_name)
    mean_file = CLIM_STATFILE_TEMPLATE.format(HINDCASTS, "MEAN", NMME_MODEL.upper(),
                                              DOMAIN_NAME, FCST_INIT_MON,
                                              CURRENTDATE.month, CURRENTDATE.day,
                                              ENDDATE.month, ENDDATE.day)
    std_file = CLIM_STATFILE_TEMPLATE.format(HINDCASTS, "STD", NMME_MODEL.upper(),
                                             DOMAIN_NAME, FCST_INIT_MON,
                                             CURRENTDATE.month, CURRENTDATE.day,
                                             ENDDATE.month, ENDDATE.day)
    logger.info(f"Reading forecast climatological mean {mean_file}", subtask=var_name)
    logger.info(f"Reading forecast climatological std {std_file}", subtask=var_name)
    mean_xr = load_ncdata(mean_file, [logger, var_name])
    std_xr = load_ncdata(std_file, [logger, var_name])

    for lead in range(LEAD_WEEKS):
        logger.info(f"Processing var_name: {var_name}, lead: {lead}", subtask=var_name)

        fcst_file = TARGET_INFILE_TEMPLATE.format(FORECASTS, CURRENTDATE.year, CURRENTDATE.month,
                                                  NMME_MODEL, NMME_MODEL.upper(), DOMAIN_NAME,
                                                  TARGET_YEAR, FCST_INIT_MON,
                                                  curdate.year, curdate.month, curdate.day,
                                                  enddate.year, enddate.month, enddate.day)
        curdate += relativedelta(days=7)
        enddate += relativedelta(days=7)
        logger.info(f"Reading target {fcst_file}", subtask=var_name)
        fcst_xr = load_ncdata(fcst_file, [logger, var_name])
        fcst_data = sel_var(fcst_xr, var_name, HYD_MODEL)
        mean_data = sel_var(mean_xr.isel(lead_time=lead), var_name, HYD_MODEL)
        std_data = sel_var(std_xr.isel(lead_time=lead), var_name, HYD_MODEL)
        this_anom = None

        # Initialize the anomaly array on first lead
        if lead == 0:
            ens_count = fcst_xr.sizes['ensemble']
            lat_count = fcst_xr.sizes['lat']
            lon_count = fcst_xr.sizes['lon']
            all_anom = np.ones((ens_count, LEAD_WEEKS, lat_count, lon_count))*-9999.

        if anom == 'ANOM':
            this_anom = xr.apply_ufunc(
                compute_anomaly,
                fcst_data.chunk({"lat": "auto", "lon": "auto"}).compute(),
                mean_data.chunk({"lat": "auto", "lon": "auto"}).compute(),
                input_core_dims=[['ensemble','time',],[]],
                exclude_dims=set(('ensemble','time',)),
                output_core_dims=[['ensemble','time',]],
                vectorize=True,  # loop over non-core dims
                dask="forbidden",
                output_dtypes=[np.float64])

        if anom == 'SANOM':
            this_anom = xr.apply_ufunc(
                compute_sanomaly,
                fcst_data.chunk({"lat": "auto", "lon": "auto"}).compute(),
                mean_data.chunk({"lat": "auto", "lon": "auto"}).compute(),
                std_data.chunk({"lat": "auto", "lon": "auto"}).compute(),
                input_core_dims=[['ensemble','time',],[],[]],
                exclude_dims=set(('ensemble','time',)),
                output_core_dims=[['ensemble','time',]],
                vectorize=True,  # loop over non-core dims
                dask="forbidden",
                output_dtypes=[np.float64])

        for ens in range(ens_count):
            all_anom[ens, lead, :, :] = this_anom[:,:,ens,0]
            lats = np.array(fcst_xr.lat.values, dtype=np.float64)
            lons = np.array(fcst_xr.lon.values, dtype=np.float64)
        del fcst_xr

    ### Step-4 Writing output file
    all_anom = np.ma.masked_array(all_anom, mask=all_anom == -9999.)

    ## Creating latitude and longitude arrays
    anom_xr = xr.Dataset()
    anom_xr[var_name.replace('-','_') + '_' + anom] = \
        (('ens', 'time', 'latitude', 'longitude'), all_anom)
    anom_xr.coords['latitude'] = (('latitude'), lats)
    anom_xr.coords['longitude'] = (('longitude'), lons)
    anom_xr.coords['time'] = (('time'), np.arange(0, LEAD_WEEKS, dtype=int))
    anom_xr.coords['ens'] = (('ens'), np.arange(0, ens_count, dtype=int))
    # Add attributes
    anom_xr.attrs['Conventions'] = 'CF-1.8'
    anom_xr['latitude'].attrs = {
        'long_name': 'latitude',
        'standard_name': 'latitude',
        'units': 'degree_north',
        'axis': 'Y'
    }
    anom_xr['longitude'].attrs = {
        'long_name': 'longitude',
        'standard_name': 'longitude',
        'units': 'degree_east',
        'axis': 'X'
    }
    anom_xr['ens'].attrs = {
        'long_name': 'Ensemble members',
        'axis': 'E',
        'units': '1'
    }
    anom_xr['time'].attrs = {
        'long_name': 'Forecast month',
        'units': 'months'
    }
    anom_xr[var_name + '_' + anom].attrs = {
        'long_name': long_names[var_name],
        'units': units[var_name]
    }
    logger.info(f"Processed variable: {var_name} for metric {anom}", subtask=var_name)
    return anom_xr

logger.info("Starting parallel processing of variables")
# ProcessPoolExecutor parallel processing
# ---------------------------------------
MAX_RETRIES = 3
datasets = []

# Track job queue
job_queue = []
for _var_name in WEEKLY_VARS:
    job_queue.append((_var_name, "ANOM"))
    job_queue.append((_var_name, "SANOM"))

RETRY_COUNT = 0

while job_queue and RETRY_COUNT <= MAX_RETRIES:
    logger.info(f"Processing batch (attempt {RETRY_COUNT + 1}), {len(job_queue)} jobs remaining")

    # Process current batch
    current_batch = job_queue.copy()
    job_queue.clear()
    failed_jobs = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_job = {}
        for _var_name, proc_type in current_batch:
            logger.info(f"Submitting {proc_type} processing job for {_var_name}", subtask=_var_name)
            future = executor.submit(process_variable, _var_name, proc_type)
            future_to_job[future] = (_var_name, proc_type)

        start_time = time.time()
        COMPLETED_COUNT = 0

        for future in as_completed(future_to_job, timeout=20*60):
            _var_name, proc_type = future_to_job[future]
            try:
                result = future.result(timeout=1)
                datasets.append(result)
                COMPLETED_COUNT += 1
                logger.info(f"✓ ({COMPLETED_COUNT}/{len(current_batch)}) {proc_type}",
                            subtask=_var_name)

            except Exception as e:
                logger.error(f"✗ Failed {proc_type} - {e}", subtask=_var_name)
                failed_jobs.append((_var_name, proc_type))

        # Handle any remaining futures
        for future, job in future_to_job.items():
            if not future.done():
                _var_name, proc_type = job
                logger.warning(f"Canceling hung job: {proc_type}", subtask=_var_name)
                future.cancel()
                failed_jobs.append((_var_name, proc_type))

    # Add failed jobs back to queue for retry
    if failed_jobs and RETRY_COUNT < MAX_RETRIES:
        job_queue.extend(failed_jobs)
        RETRY_COUNT += 1
        logger.info(f"Retrying {len(failed_jobs)} failed jobs (attempt {RETRY_COUNT + 1})")
        time.sleep(5)
    elif failed_jobs:
        logger.error(f"Max retries reached. {len(failed_jobs)} jobs permanently failed: "
                     f"{failed_jobs}")
    else:
        logger.info("All jobs completed successfully")
        break

# Merge all datasets
merged_dataset = xr.merge(datasets)
ENCODING = {var: comp for var in merged_dataset.data_vars}

# Write the merged dataset to a single file
logger.info(f"Writing merged output to {OUTFILE}")
write_ncfile(merged_dataset, OUTFILE, ENCODING, [logger, ''])
