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
import sys
from datetime import datetime
import calendar
from time import ctime as t_ctime
from time import time as t_time
from dateutil.relativedelta import relativedelta
import numpy as np
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
# pylint: enable=no-name-in-module
# pylint: disable=import-error
from shrad_modules import read_nc_files
# pylint: enable=import-error

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
    ### Assigning standard names for each variables
    varname1.standard_name = var_standard_name[0]
    varname2.standard_name = var_standard_name[1]
    varname3.standard_name = var_standard_name[2]
    varname4.standard_name = var_standard_name[3]
    varname5.standard_name = var_standard_name[4]
    varname6.standard_name = var_standard_name[5]

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
    times[:] = nc4_date2num(dates, units=times.units, calendar=times.calendar)
    rootgrp.close()

## Usage: <Name of variable in observed climatology> <Name of variable in
## reforecast climatology (same as the name in target forecast> <forecast
## model number>
CMDARGS = str(sys.argv)
INIT_FCST_YEAR = int(sys.argv[1])
## initial forecast year for which to downscale the data
INIT_FCST_MON = int(sys.argv[2])
## initial forecast month for which to downscale the data

MODEL_NAME = str(sys.argv[3])
ENS_NUM = int(sys.argv[4])
LEAD_FINAL = int(sys.argv[5])
MONTH_NAME_TEMPLATE = '{}01'
MONTH_NAME = MONTH_NAME_TEMPLATE.format(calendar.month_abbr[INIT_FCST_MON])

#Directory and file addresses
BASEDIR = str(sys.argv[6])
INDIR_TEMPLATE = '{}/bcsd/6-Hourly/{}/{:04d}/ens{:01d}'
#### Change the model name here for other models
INFILE_TEMPLATE = '{}/{}.{:04d}{:02d}.nc4'

OUTDIR_TEMPLATE = '{}/final/6-Hourly/{}/{:04d}/ens{:01d}'
OUTFILE_TEMPLATE = '{}/CFSv2.{:04d}{:02d}.nc4'

#VAR_NAME_LIST=['LWS', 'SLRSF', 'PS', 'Q2M', 'T2M', 'WIND10M']
VAR_NAME_LIST = ['LWGAB', 'SWGDN', 'PS', 'QV2M', 'T2M', 'U10M']
UNITS = ['W/m^2', 'W/m^2', 'Pa', 'kg/kg', 'K', 'm/s']

for MON in [INIT_FCST_MON]:
    MONTH_NAME = MONTH_NAME_TEMPLATE.format((calendar.month_abbr[MON]).lower())
    ## This provides abbrevated version of the name of a month: (e.g. for
    ## January (i.e. Month number = 1) it will return "Jan"). The abbrevated
    ## name is used in the forecasts file name
    print(f"Forecast Initialization month is {MONTH_NAME}")
    ## Shape of the above dataset time, Lead, Ens, latitude, longitude
    for ens in range(ENS_NUM):
        INDIR = INDIR_TEMPLATE.format(BASEDIR, MONTH_NAME, \
        INIT_FCST_YEAR, ens+1)
        OUTDIR = OUTDIR_TEMPLATE.format(BASEDIR, MONTH_NAME, \
        INIT_FCST_YEAR, ens+1)
        if os.path.isdir(OUTDIR):
            pass
        else:
            os.makedirs(OUTDIR, exist_ok=True)
        print(f"OUTDIR is {OUTDIR}")
        for LEAD_NUM in range(0, LEAD_FINAL):
            ## Loop from lead =0 to Final Lead
            FCST_DATE = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1) + \
            relativedelta(months=LEAD_NUM)
            FCST_YEAR, FCST_MONTH = FCST_DATE.year, FCST_DATE.month
            for VAR_NUM, VAR_VALUE in  enumerate(VAR_NAME_LIST):
                VAR = VAR_NAME_LIST[VAR_NUM]
                INFILE = INFILE_TEMPLATE.format(INDIR, VAR, FCST_YEAR, \
                FCST_MONTH)
                TEMP = read_nc_files(INFILE, VAR)
                if VAR == VAR_NAME_LIST[0]:
                    LATS, LONS = read_nc_files(INFILE, 'lat'), \
                    read_nc_files(INFILE, 'lon')
                    IN_DATA = np.empty((len(VAR_NAME_LIST), TEMP.shape[0], \
                    len(LATS), len(LONS)))
                IN_DATA[VAR_NUM, ] = TEMP

            ### Finished reading all files now writing combined output
            OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, FCST_YEAR, FCST_MONTH)
            print(f"Now writing {OUTFILE}")
            SDATE = datetime(FCST_YEAR, FCST_MONTH, 1, 6)
            NUM_DAYS = TEMP.shape[0]
            DATES = [SDATE+relativedelta(hours=n*6) for n in range(NUM_DAYS)]
            write_bc_netcdf(OUTFILE, IN_DATA, VAR_NAME_LIST, \
            'Bias corrected forecasts', 'MODEL:'  + MODEL_NAME, \
            UNITS, VAR_NAME_LIST, LONS, LATS, SDATE, DATES, 8, LATS[-1], \
            LONS[-1], LATS[0], LONS[0], 0.25, 0.25, 21600)
