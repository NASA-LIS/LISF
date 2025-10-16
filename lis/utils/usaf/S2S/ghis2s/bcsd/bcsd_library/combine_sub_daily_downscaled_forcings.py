#!/usr/bin/env python
"""
# Author: Shrad Shukla
# coding: utf-8
#Author: Shrad Shukla
#Usage: This is a module for the BCSD code.
#This module bias corrects a forecasts following probability
#mapping approach as described in Wood et al. 2002
#Date: August 06, 2015
# In[28]:
"""

import os
import subprocess
import sys
import calendar
from time import ctime as t_ctime
from time import time as t_time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
# pylint: enable=no-name-in-module
from ghis2s.shared.utils import load_ncdata, get_domain_info
from ghis2s.shared.logging_utils import TaskLogger


def write_bc_netcdf(outfile, var, varname, description, source, var_units, \
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
    varname1 = rootgrp.createVariable(varname[0], 'f4', ('time', \
    'latitude', 'longitude',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    varname2 = rootgrp.createVariable(varname[1], 'f4', ('time', \
    'latitude', 'longitude',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    varname3 = rootgrp.createVariable(varname[2], 'f4', ('time', \
    'latitude', 'longitude',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    varname4 = rootgrp.createVariable(varname[3], 'f4', ('time', \
    'latitude', 'longitude',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    varname5 = rootgrp.createVariable(varname[4], 'f4', ('time', \
    'latitude', 'longitude',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    varname6 = rootgrp.createVariable(varname[5], 'f4', ('time', \
    'latitude', 'longitude',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    varname7 = rootgrp.createVariable('V10M', 'f4', ('time', \
    'latitude', 'longitude',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)

    rootgrp.missing_value = -9999
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
    #rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.history = 'Created ' + t_ctime(t_time())
    rootgrp.source = source
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    ### Assigning units for each variables
    varname1.units = var_units[0]
    varname2.units = var_units[1]
    varname3.units = var_units[2]
    varname4.units = var_units[3]
    varname5.units = var_units[4]
    varname6.units = var_units[5]
    varname7.units = 'm/s'

    ### Assigning standard names for each variables
    varname1.standard_name = var_standard_name[0]
    varname2.standard_name = var_standard_name[1]
    varname3.standard_name = var_standard_name[2]
    varname4.standard_name = var_standard_name[3]
    varname5.standard_name = var_standard_name[4]
    varname6.standard_name = var_standard_name[5]
    varname7.standard_name = 'V10M'

    string_date = datetime.strftime(sdate, "%Y-%m-%d %H:%M:%S")
    times.units = 'minutes since ' + string_date
    times.time_increment = time_increment
    times.begin_date = datetime.strftime(sdate, "%Y%m%d")
    times.begin_time = '000000'
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    ## Passing on values
    varname1[:, :, :] = var[0, ]
    varname2[:, :, :] = var[1, ]
    varname3[:, :, :] = var[2, ]
    varname4[:, :, :] = var[3, ]
    varname5[:, :, :] = var[4, ]
    varname6[:, :, :] = var[5, ]
    varname7[:, :, :] = np.zeros_like(var[5,])
    times[:] = nc4_date2num(dates, units=times.units, calendar=times.calendar)
    rootgrp.close()

## Usage: <Name of variable in observed climatology> <Name of variable in
## reforecast climatology (same as the name in target forecast> <forecast
## model number>
task_name = os.environ.get('SCRIPT_NAME')
logger = TaskLogger(task_name,
                    os.getcwd(),
                    f'bcsd/bcsd_library/combine_sub_daily_downscaled_forcings.py processing {sys.argv[3]} for month {int(sys.argv[1]):04d}{int(sys.argv[2]):02d}')

INIT_FCST_YEAR = int(sys.argv[1])
## initial forecast year for which to downscale the data
INIT_FCST_MON = int(sys.argv[2])
## initial forecast month for which to downscale the data

MODEL_NAME = str(sys.argv[3])
ENS_NUM = int(sys.argv[4])
LEAD_FINAL = int(sys.argv[5])
MONTH_NAME_TEMPLATE = '{}01'
MONTH_NAME = MONTH_NAME_TEMPLATE.format(calendar.month_abbr[INIT_FCST_MON].lower())

#Directory and file addresses
BASEDIR = str(sys.argv[6])
CONFIG_FILE = str(sys.argv[7])
INDIR_TEMPLATE = '{}/bcsd/6-Hourly/{}/{:04d}/ens{:01d}'
INFILE_TEMPLATE = '{}/{}.{:04d}{:02d}.nc4'
OUTDIR_TEMPLATE = '{}/final/6-Hourly/{}/{:04d}/ens{:01d}'
OUTFILE_TEMPLATE = '{}/{}.{:04d}{:02d}.nc4'

VAR_NAME_LIST = ['LWGAB', 'SWGDN', 'PS', 'QV2M', 'T2M', 'U10M']
UNITS = ['W/m^2', 'W/m^2', 'Pa', 'kg/kg', 'K', 'm/s']

latsg, _ = get_domain_info(CONFIG_FILE, coord=True)
resol = round((latsg[1] - latsg[0])*100.)/100.

def process_ensemble(_ens):
    ''' process each ensemble member '''
    subtask = f'ens{_ens:02d}'

    ## This provides abbrevated version of the name of a month: (e.g. for
    ## January (i.e. Month number = 1) it will return "Jan"). The abbrevated
    ## name is used in the forecasts file name
    logger.info(f"Forecast Initialization month is {MONTH_NAME}", subtask=subtask)
    ## Shape of the above dataset time, Lead, Ens, latitude, longitude
    #for ens in range(ENS_NUM):
    indir = INDIR_TEMPLATE.format(BASEDIR, MONTH_NAME, \
                                  INIT_FCST_YEAR, _ens+1)
    outdir = OUTDIR_TEMPLATE.format(BASEDIR, MONTH_NAME, \
                                    INIT_FCST_YEAR, _ens+1)
    if os.path.isdir(outdir):
        pass
    else:
        os.makedirs(outdir, exist_ok=True)

    for lead_num in range(0, LEAD_FINAL):
        ## Loop from lead =0 to Final Lead
        fcst_date = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1) + \
            relativedelta(months=lead_num)
        fcst_year, fcst_month = fcst_date.year, fcst_date.month
        for var_num, _ in  enumerate(VAR_NAME_LIST):
            var = VAR_NAME_LIST[var_num]
            infile = INFILE_TEMPLATE.format(indir, var, fcst_year, \
                                            fcst_month)
            logger.info(f"Reading {infile}", subtask=subtask)
            temp_da = load_ncdata(infile, [logger, subtask],  var_name=var, decode_cf=False)
            temp = temp_da.values

            if var == VAR_NAME_LIST[0]:
                lats, lons = temp_da.lat.values, temp_da.lon.values
                in_data = np.empty((len(VAR_NAME_LIST), temp.shape[0], \
                                    len(lats), len(lons)))
            in_data[var_num, ] = temp
            temp_da.close()
            del temp_da

        ### Finished reading all files now writing combined output
        outfile = OUTFILE_TEMPLATE.format(outdir, MODEL_NAME, fcst_year, fcst_month)
        logger.info(f"Writing {outfile}", subtask=subtask)
        sdate = datetime(fcst_year, fcst_month, 1, 6)
        num_days = temp.shape[0]
        dates = [sdate+relativedelta(hours=n*6) for n in range(num_days)]
        force_dt = 21600
        if MODEL_NAME == 'CFSv2':
            force_dt = 21600
        if MODEL_NAME == 'GEOSv3':
            force_dt = 10800
        write_bc_netcdf(outfile, in_data, VAR_NAME_LIST, \
                        'Bias corrected forecasts', 'MODEL:'  + MODEL_NAME, \
                        UNITS, VAR_NAME_LIST, lons, lats, sdate, dates, 8, lats[-1], \
                        lons[-1], lats[0], lons[0], resol, resol, force_dt)

logger.info("Starting parallel processing of ensembles")
num_workers = int(sys.argv[4])
# ProcessPoolExecutor parallel processing
with ProcessPoolExecutor(max_workers=num_workers) as executor:
    futures = []
    for ens in range(int(sys.argv[4])):
        logger.info(f"Submitting disaggregation job for ens {ens:02d}", subtask=f'ens{ens:02d}')
        future = executor.submit(process_ensemble, ens)
        futures.append(future)

    for future in futures:
        result = future.result()

# Create ens13, ens14, ens15 links
OUTDIR_TEMPLATE = '{}/final/6-Hourly/{}/{:04d}/'
OUTDIR = OUTDIR_TEMPLATE.format(BASEDIR, MONTH_NAME.lower(), INIT_FCST_YEAR)

os.chdir(OUTDIR)
logger.info(f"Creating ens13, ens14, ens15 links in {OUTDIR}")
CMD = "ln -sfn ens1 ens13"
RC = subprocess.call(CMD, shell=True)
CMD = "ln -sfn ens2 ens14"
RC = subprocess.call(CMD, shell=True)
CMD = "ln -sfn ens3 ens15"
RC = subprocess.call(CMD, shell=True)

logger.info(f"Creating symbolic links for month {LEAD_FINAL +1}")
init_datetime = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1)
src_yyyymm = []
dst_yyyymm = []
for mon in range(LEAD_FINAL):
    src_yyyymm.append((init_datetime + relativedelta(months=mon)).strftime("%Y%m"))
    dst_yyyymm.append((init_datetime + relativedelta(months=mon)).strftime("%Y%m"))

src_yyyymm.append((init_datetime + relativedelta(months=mon)).strftime("%Y%m"))
dst_yyyymm.append((init_datetime + relativedelta(months=mon+1)).strftime("%Y%m"))
LAST_YYYYMM = len(src_yyyymm) -1
for iens, ens_value in enumerate(range(ENS_NUM)):
    ens_nmme = iens + 1
    OUTDIR_ENS = OUTDIR + f'ens{ens_nmme}'
    src_file = f"{OUTDIR_ENS}/{MODEL_NAME}.{src_yyyymm[LAST_YYYYMM]}.nc4"
    dst_file = f"{OUTDIR_ENS}/{MODEL_NAME}.{dst_yyyymm[LAST_YYYYMM]}.nc4"
    CMD = f"ln -sfn {src_file} {dst_file}"
    RC = subprocess.call(CMD, shell=True)
    if RC != 0:
        logger.error("Problem calling creating last precip symbolic link to {dst_file}!")

logger.info("Ran SUCCESSFULLY !")
