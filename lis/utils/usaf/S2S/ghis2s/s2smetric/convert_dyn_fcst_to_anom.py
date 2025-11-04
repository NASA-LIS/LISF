#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: convert_dyn_fcst_to_anom.py
#
# PURPOSE: Calculates anomaly of NMME-forced LIS FORECASTS.
#
# REVISION HISTORY:
# ?? Mar 2017: Shrad Shukla/UCSB, first version.
# ?? ??? 2019: Abheera Hazra/UMD, second version.
# 22 Oct 2021: Eric Kemp/SSAI, updated for 557WW.
# 30 Oct 2021: Eric Kemp/SSAI, updated to use s2smetric CONFIG file.
# 02 Jun 2023: K. Arsenault + S. Mahanama, updated 557 WW file conventions.
#
#------------------------------------------------------------------------------
"""

# Standard modules
from datetime import datetime, date
import glob
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import yaml

# Third-party modules
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
# Local modules
from ghis2s.shared.utils import write_ncfile, load_ncdata
from ghis2s.shared.logging_utils import TaskLogger
from ghis2s.s2smetric.metricslib import (sel_var, compute_anomaly, compute_sanomaly,
                                         merged_metric_filename, LONG_NAMES_ANOM,
                                         LONG_NAMES_SANOM, UNITS_ANOM, UNITS_SANOM)
# pylint: disable=too-many-positional-arguments
# pylint: disable=too-many-arguments,too-many-locals

# Start reading from command line.
FCST_INIT_MON = int(sys.argv[1])
TARGET_YEAR = int(sys.argv[2])
NMME_MODEL = sys.argv[3]
CONFIGFILE = sys.argv[4]
BASEOUTDIR = sys.argv[5]

# Load CONFIG file
with open(CONFIGFILE, 'r', encoding="utf-8") as file:
    CONFIG = yaml.safe_load(file)
HYD_MODEL = CONFIG["EXP"]["lsmdir"]
LEAD_NUM = int(CONFIG["EXP"]["lead_months"])
DOMAIN_NAME = CONFIG["EXP"]["DOMAIN"]
CLIM_SYR = int(CONFIG["BCSD"]["clim_start_year"])
CLIM_EYR = int(CONFIG["BCSD"]["clim_end_year"])
METRIC_VARS = CONFIG["POST"]["metric_vars"]
HINDCASTS = CONFIG["SETUP"]["E2ESDIR"] + '/hindcast/s2spost/' + f'{FCST_INIT_MON:02d}/'
FORECASTS = "./s2spost/"
CURRENTDATE = date(TARGET_YEAR, FCST_INIT_MON, 1)
ENDDATE = CURRENTDATE
ENDDATE += relativedelta(months=int(CONFIG["EXP"]["lead_months"]))
FCST_INIT_DAY = 1
OUTDIR = BASEOUTDIR + '/'
os.makedirs(OUTDIR, exist_ok=True)

# Name of variable, forecast initial month and forecast year is in the file

TARGET_INFILE_TEMPLATE1 = ""
TARGET_INFILE_TEMPLATE2 = ""
CLIM_INFILE_TEMPLATE1 = ""
CLIM_INFILE_TEMPLATE2 = ""

if DOMAIN_NAME == 'AFRICOM':
    TARGET_INFILE_TEMPLATE1 = \
        '{}/{:04d}{:02d}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_'
    CLIM_INFILE_TEMPLATE1 = \
        '{}/????{:02d}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_'
elif DOMAIN_NAME == 'GLOBAL':
    TARGET_INFILE_TEMPLATE1 = \
        '{}/{:04d}{:02d}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.GLOBAL_'
    CLIM_INFILE_TEMPLATE1 = \
        '{}/????{:02d}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.GLOBAL_'

TARGET_INFILE_TEMPLATE2 = \
    'PA.ALL_DD.{:04d}{:02d}01_DT.0000_FP.{:04d}{:02d}{:02d}-{:04d}{:02d}01_DF.NC'

CLIM_INFILE_TEMPLATE2 = 'PA.ALL_DD.*{:02d}01_DT.0000_FP.*{:02d}{:02d}-*{:02d}01_DF.NC'

TARGET_INFILE_TEMPLATE = TARGET_INFILE_TEMPLATE1 + TARGET_INFILE_TEMPLATE2
CLIM_INFILE_TEMPLATE = CLIM_INFILE_TEMPLATE1 + CLIM_INFILE_TEMPLATE2
## String in this format allows to select all the files for the given month

num_vars = 2*len(METRIC_VARS)
num_workers = int(os.environ.get('NUM_WORKERS', num_vars))

OUTFILE = merged_metric_filename(OUTDIR, CURRENTDATE, ENDDATE, NMME_MODEL, DOMAIN_NAME)
task_name = os.environ.get('SCRIPT_NAME')
logger = TaskLogger(task_name,
                    os.getcwd(),
                    f'{NMME_MODEL} running s2smetric/convert_dyn_fcast_to_anom.py')

def process_variable(var_name, metric_name):
    """
    This routine processes each variable and metric (e.g., anomalies).
    """
    logger.info(f"Starting processing {var_name} for metric {metric_name}", subtask=var_name)
    if metric_name == "ANOM":
        long_names = LONG_NAMES_ANOM
        units = UNITS_ANOM
    else:
        long_names = LONG_NAMES_SANOM
        units = UNITS_SANOM

    for lead in range(LEAD_NUM):
        logger.info(f"Reading output from Hindcast runs - lead: {lead}", subtask=var_name)
        logger.info(f"Processing var_name: {var_name}, lead: {lead}", subtask=var_name)

        ## Step-1: Read and process the climatology
        smon = datetime(CLIM_SYR, FCST_INIT_MON, FCST_INIT_DAY) + \
                    relativedelta(months=lead)
        emon = datetime(CLIM_SYR, FCST_INIT_MON, FCST_INIT_DAY) + \
                    relativedelta(months=lead+1)
        ## Adding 1 to lead to make sure the file read is from the month after
        smon1 = datetime(TARGET_YEAR, FCST_INIT_MON, FCST_INIT_DAY) + \
            relativedelta(months=lead)
        emon1 = datetime(TARGET_YEAR, FCST_INIT_MON, FCST_INIT_DAY) + \
            relativedelta(months=lead+1)

        if lead == 0:
            infile = CLIM_INFILE_TEMPLATE.format(HINDCASTS,
                                                 FCST_INIT_MON, NMME_MODEL,
                                                 NMME_MODEL.upper(),
                                                 FCST_INIT_MON,
                                                 smon.month, 2, emon.month)
            tinfile = TARGET_INFILE_TEMPLATE.format(FORECASTS,
                                                    TARGET_YEAR, FCST_INIT_MON, NMME_MODEL,
                                                    NMME_MODEL.upper(),
                                                    TARGET_YEAR, FCST_INIT_MON,
                                                    smon1.year, smon1.month, 2,
                                                    emon1.year, emon1.month)
        else:
            infile = CLIM_INFILE_TEMPLATE.format(HINDCASTS,
                                                 FCST_INIT_MON, NMME_MODEL,
                                                 NMME_MODEL.upper(),
                                                 FCST_INIT_MON,
                                                 smon.month, 1, emon.month)
            tinfile = TARGET_INFILE_TEMPLATE.format(FORECASTS,
                                                    TARGET_YEAR, FCST_INIT_MON, NMME_MODEL,
                                                    NMME_MODEL.upper(),
                                                    TARGET_YEAR, FCST_INIT_MON,
                                                    smon1.year, smon1.month, 1,
                                                    emon1.year, emon1.month)

        logger.info(f"Reading forecast climatology {infile}", subtask=var_name)
        infile1 = glob.glob(infile)

        # First reading all available years for the given
        # forecast initialization month
        all_clim_data1 = load_ncdata(infile1, [logger, var_name], combine='by_coords')

        # Now selecting only the years that are within the climatology
        sel_cim_data = all_clim_data1.sel(time= \
                       (all_clim_data1.coords['time.year'] >= \
                       CLIM_SYR) & (all_clim_data1.coords['time.year'] <= \
                       CLIM_EYR))

        # Now selecting the climatology further for the given variable
        all_clim_data = sel_var(sel_cim_data, var_name, HYD_MODEL)

        # all_clim_data has all the needed climatology for a given variable
        # To-DO: In future we may want to save the climatology all
        # together in a file so we don't have to read climatologies every time

        ####### Step-2: Read the target forecast which needs to be converted
        ## into anomaly

        logger.info(f"Reading target {tinfile}", subtask=var_name)

        # Note target will always have only one time step
        target_data = load_ncdata(tinfile, [logger, var_name])

        ## Now selecting the desired variable
        target_fcst_data = sel_var(target_data, var_name, HYD_MODEL)
        target_fcst_data = target_fcst_data.load()
        all_clim_data = all_clim_data.load()

        ## Step-3 loop through each grid cell and convert data into anomaly
        # Defining array to store anomaly data
        lat_count, lon_count, ens_count = \
            len(target_data.coords['lat']), \
            len(target_data.coords['lon']), \
            len(target_data.coords['ensemble'])

        ## Note that ens_count is coming from the Target FORECASTS,
        ## so if the target_FORECASTS have 4 members there will be
        ## 4 members in anomaly output and so on.
        if lead == 0:
            all_anom = np.ones((ens_count, LEAD_NUM, lat_count, lon_count))*-9999.

        logger.info('Converting data into anomaly', subtask=var_name)
        all_clim_mean = all_clim_data.mean (dim = ['time','ensemble'], skipna = True)

        if metric_name == "SANOM":
            all_clim_std  = all_clim_data.std (dim = ['time','ensemble'], skipna = True)

        if (not np.array_equal(all_clim_mean.lat.values, target_fcst_data.lat.values)) or \
           (not np.array_equal(all_clim_mean.lon.values, target_fcst_data.lon.values)):
            all_clim_mean = all_clim_mean.assign_coords({"lon": target_fcst_data.lon.values,
                                                         "lat": target_fcst_data.lat.values})
            if metric_name == "SANOM":
                all_clim_std = \
                    all_clim_std.assign_coords({"lon": target_fcst_data.lon.values,
                                                "lat": target_fcst_data.lat.values})
        if metric_name == "ANOM":
            this_anom = xr.apply_ufunc(
                compute_anomaly,
                target_fcst_data.chunk({"lat": "auto", "lon": "auto"}).compute(),
                all_clim_mean.chunk({"lat": "auto", "lon": "auto"}).compute(),
                input_core_dims=[['ensemble','time',],[]],
                exclude_dims=set(('ensemble','time',)),
                output_core_dims=[['ensemble','time',]],
                vectorize=True,  # loop over non-core dims
                dask="forbidden",
                output_dtypes=[np.float64])
        else:
            this_anom = xr.apply_ufunc(
                compute_sanomaly,
                target_fcst_data.chunk({"lat": "auto", "lon": "auto"}).compute(),
                all_clim_mean.chunk({"lat": "auto", "lon": "auto"}).compute(),
                all_clim_std.chunk({"lat": "auto", "lon": "auto"}).compute(),
                input_core_dims=[['ensemble','time',],[],[]],
                exclude_dims=set(('ensemble','time',)),
                output_core_dims=[['ensemble','time',]],
                vectorize=True,  # loop over non-core dims
                dask="forbidden",
                output_dtypes=[np.float64])

        for ens in range(ens_count):
            all_anom[ens, lead, :, :] = this_anom [:,:,ens,0]

        del all_clim_data, target_fcst_data, all_clim_mean

    ### Step-4 Writing output file
    all_anom = np.ma.masked_array(all_anom, mask=all_anom == -9999.)

    ## Creating an latitude and longitude array based on locations of corners
    lats = np.arange(target_data.attrs['SOUTH_WEST_CORNER_LAT'], \
                     target_data.attrs['SOUTH_WEST_CORNER_LAT'] + \
                     (lat_count*0.25), 0.25)
    lons = np.arange(target_data.attrs['SOUTH_WEST_CORNER_LON'], \
                     target_data.attrs['SOUTH_WEST_CORNER_LON'] + \
                     (lon_count*0.25), 0.25)

    anom_xr = xr.Dataset()
    anom_xr[var_name.replace('-','_') + '_' + metric_name] = \
        (('ens', 'time', 'latitude', 'longitude'), all_anom)
    anom_xr.coords['latitude'] = (('latitude'), lats)
    anom_xr.coords['longitude'] = (('longitude'), lons)
    anom_xr.coords['time'] = (('time'), np.arange(0, LEAD_NUM, dtype=int))
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
    anom_xr[var_name + '_' + metric_name].attrs = {
        'long_name': long_names[var_name],
        'units': units[var_name]
    }
    logger.info(f"Processed variable: {var_name} for metric {metric_name}", subtask=var_name)
    return anom_xr

logger.info("Starting parallel processing of variables")

# ProcessPoolExecutor parallel processing
# ---------------------------------------
MAX_RETRIES = 3
datasets = []

# Track job queue
job_queue = []
for _var_name in METRIC_VARS:
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
logger.info("Merging all datasets")
merged_dataset = xr.merge(datasets)
comp = {'zlib':True, 'complevel':6, 'shuffle':True, '_FillValue': -9999.}
encoding = {var: comp for var in merged_dataset.data_vars}

# Write the merged dataset to a single file
logger.info(f"Writing merged output to {OUTFILE}")
write_ncfile(merged_dataset, OUTFILE, encoding, [logger, ''])
