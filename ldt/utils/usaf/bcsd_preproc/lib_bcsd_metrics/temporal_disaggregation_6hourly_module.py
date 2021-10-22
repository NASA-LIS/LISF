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

from __future__ import division
#import pandas as pdimport calendar
import os.path as op
import sys
import time
from datetime import datetime, timedelta
import calendar
from dateutil.relativedelta import relativedelta
import numpy as np
#from scipy.stats import percentileofscore
#from scipy.stats import scoreatpercentile, pearsonr
#from math import *
#import time
#from BCSD_stats_functions import *
# pylint: disable=no-name-in-module
import netCDF4 as nc
# pylint: enable=no-name-in-module
from Shrad_modules import read_nc_files, MAKEDIR

def write_bc_netcdf(outfile, var, varname, description, source, var_units, \
var_standard_name, lons, lats, sdate, dates, sig_digit, north_east_corner_lat, \
north_east_corner_lon, south_west_corner_lat, south_west_corner_lon, \
resolution_x, resolution_y, time_increment):
    """write netcdf"""
    rootgrp = nc.Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
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
    import time
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
    rootgrp.history = 'Created ' + time.ctime(time.time())
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
    times[:] = nc.date2num(dates, units=times.units, calendar=times.calendar)
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
LAT1, LAT2, LON1, LON2 = int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]), int(sys.argv[10])


MODEL_NAME = str(sys.argv[11])
ENS_NUM = int(sys.argv[12])
LEAD_FINAL = int(sys.argv[13])
MONTH_NAME_TEMPLATE = '{}01'
MONTH_NAME = MONTH_NAME_TEMPLATE.format(calendar.month_abbr[INIT_FCST_MON])

print("*** LEAD FINAL: {}".format(LEAD_FINAL))
BC_FCST_SYR, BC_FCST_EYR = int(sys.argv[14]), int(sys.argv[15])
if FCST_VAR == 'PRECTOT':
    MASK_FILE = str(sys.argv[16])
else:
    MASK_FILE = str(sys.argv[17])

MASK = read_nc_files(MASK_FILE, 'mask')
print("MASK: {}".format(MASK.shape))

MONTHLY_BC_FCST_DIR = str(sys.argv[18])
MONTHLY_RAW_FCST_DIR = str(sys.argv[19])
SUBDAILY_RAW_FCST_DIR = str(sys.argv[20])
BASE_OUTDIR = str(sys.argv[21])
OUTDIR_TEMPLATE = '{}/{:04d}/ens{:01d}'


# All file formats
MONTHLY_BC_INFILE_TEMPLATE = '{}/{}.{}.{}_{:04d}_{:04d}.nc'
MONTHLY_RAW_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.cfsv2.{:04d}{:02d}.nc'
SUBDAILY_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.cfsv2.{:04d}{:02d}.nc'
SUBDAILY_OUTFILE_TEMPLATE = '{}/{}.{:04d}{:02d}.nc4'
MONTHLY_NMME_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.nmme.monthly.{:04d}{:02d}.nc'

for MON in [INIT_FCST_MON]:
    MONTH_NAME = MONTH_NAME_TEMPLATE.format((calendar.month_abbr[MON]).lower())
    ## This provides abbrevated version of the name of a month: (e.g. for
    ## January (i.e. Month number = 1) it will return "Jan"). The abbrevated
    ## name is used in the forecasts file name
    print("Forecast Initialization month is {}".format(MONTH_NAME))
    ### First read bias corrected monthly forecast data
    BC_INFILE = MONTHLY_BC_INFILE_TEMPLATE.format(MONTHLY_BC_FCST_DIR,\
    FCST_VAR, MODEL_NAME, MONTH_NAME, BC_FCST_SYR, BC_FCST_EYR)

    print("Reading bias corrected monthly forecasts {}".format(BC_INFILE))
    LATS = read_nc_files(BC_INFILE, 'latitude')
    LONS = read_nc_files(BC_INFILE, 'longitude')
    MON_BC_DATA = read_nc_files(BC_INFILE, FCST_VAR)

    ## Shape of the above dataset time, Lead, Ens, latitude, longitude
    for ens in range(ENS_NUM):
        OUTDIR = OUTDIR_TEMPLATE.format(BASE_OUTDIR, INIT_FCST_YEAR, ens+1)
        if op.isdir(OUTDIR):
            pass
        else:
            MAKEDIR(OUTDIR)
        print("OUTDIR is {}".format(OUTDIR))
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
            OUTPUT_BC_DATA = np.ones((NUM_TIMESTEPS, len(LATS), len(LONS)))*-999
            # Monthly raw data
            if FCST_VAR != 'PRECTOT':
                MONTHLY_INFILE = MONTHLY_RAW_INFILE_TEMPLATE.format(\
                MONTHLY_RAW_FCST_DIR, INIT_FCST_YEAR, ens+1, MONTH_NAME, \
                FCST_YEAR, FCST_MONTH)
            else:
                print("Temporarily using nmme unique TEMPLATE. \
                Recode in future.")
                MONTHLY_INFILE = MONTHLY_NMME_INFILE_TEMPLATE.format(\
                MONTHLY_RAW_FCST_DIR, INIT_FCST_YEAR, ens+1, MONTH_NAME, \
                FCST_YEAR, FCST_MONTH)
            print("Reading raw monthly forecast {}".format(MONTHLY_INFILE))
            MONTHLY_INPUT_RAW_DATA = read_nc_files(MONTHLY_INFILE, FCST_VAR)[0,]
            # Sub-Daily raw data
            SUBDAILY_INFILE = SUBDAILY_INFILE_TEMPLATE.format(\
            SUBDAILY_RAW_FCST_DIR, INIT_FCST_YEAR, ens+1, MONTH_NAME, \
            FCST_YEAR, FCST_MONTH)
            print("Reading raw sub-daily forecast {}".format(SUBDAILY_INFILE))
            INPUT_RAW_DATA = read_nc_files(SUBDAILY_INFILE, FCST_VAR)
            print("MONTHLY_BC_DATA: {}".format(MON_BC_DATA.shape))
            print("MONTHLY_INPUT_RAW_DATA: {}".format(\
            MONTHLY_INPUT_RAW_DATA.shape))
            print("INPUT_RAW_DATA: {}".format(INPUT_RAW_DATA.shape))
            print("OUTPUT_BC_DATA: {}".format(OUTPUT_BC_DATA.shape))
            #for lat_num, LATS in  enumerate(LATS):
            for lat_num, lat_value in  enumerate(LATS):
                #for lon_num in LONS in enumerate(LONS):
                for lon_num, lon_value in enumerate(LONS):
                    ## Only work with grid cells that are within the given mask
                    if ((LAT1 <= LATS[lat_num]) and (LATS[lat_num] <= LAT2) \
                        and (LON1 <= LONS[lon_num]) and \
                        (LONS[lon_num] <= LON2)):
                        # Bias corrected monthly value
                        MON_BC_VALUE = MON_BC_DATA[(\
                        INIT_FCST_YEAR-BC_FCST_SYR), LEAD_NUM, ens, \
                        lat_num, lon_num]
                        # Raw Monthly value
                        MON_RAW_VALUE = MONTHLY_INPUT_RAW_DATA[lat_num, lon_num]
                        if BC_VAR == 'PRCP':
                            if MON_RAW_VALUE == 0:
                                CORRECTION_FACTOR = MON_BC_VALUE
                                ## HACK## for when input monthly value is 0
                                OUTPUT_BC_DATA[:, lat_num, lon_num] = \
                                CORRECTION_FACTOR
                            else:
                                CORRECTION_FACTOR = MON_BC_VALUE/MON_RAW_VALUE
                                OUTPUT_BC_DATA[:, lat_num, lon_num] = \
                                INPUT_RAW_DATA[:, lat_num, \
                                lon_num]*CORRECTION_FACTOR
                        else:
                            CORRECTION_FACTOR = MON_BC_VALUE-MON_RAW_VALUE
                            OUTPUT_BC_DATA[:, lat_num, lon_num] = \
                            INPUT_RAW_DATA[:, lat_num, lon_num]+\
                            CORRECTION_FACTOR
            ### Finish correcting values for all timesteps in the given
            ### month and ensemble member
            print("Now writing {}".format(OUTFILE))
            OUTPUT_BC_DATA = np.ma.masked_array(OUTPUT_BC_DATA, \
            mask=OUTPUT_BC_DATA == -999)
            date = [FCST_DATE+relativedelta(hours=n*6) for n in \
            range(NUM_TIMESTEPS)]
            write_bc_netcdf(OUTFILE, OUTPUT_BC_DATA, OBS_VAR, \
            'Bias corrected forecasts', 'MODEL:'  +   MODEL_NAME, UNIT, \
            OBS_VAR, LONS, LATS, FCST_DATE, date, 5, 39.875, 59.875, -39.875, \
            -19.875, 0.25, 0.25, 21600)
