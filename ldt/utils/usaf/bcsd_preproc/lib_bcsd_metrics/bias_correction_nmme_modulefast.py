#!/usr/bin/env python
"""
# Author: Shrad Shukla
# coding: utf-8
#Author: Shrad Shukla
#Usage: This is a module for the BCSD code.
#This module bias corrects a forecasts following
#probability mapping approach as described in Wood et al. 2002
#Date: August 06, 2015
# In[28]:
"""

from __future__ import division
#import pandas as pd
import os
import sys
import calendar
#import os.path as op
from datetime import datetime
import numpy as np
from dateutil.relativedelta import relativedelta
#from math import *
#import time
import xarray as xr
import BCSD_function
from BCSD_stats_functions import write_4d_netcdf
from Shrad_modules import read_nc_files

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

def slice_latlon(lat, lon, lat_range: list, lon_range: list):
    """
      Function for extracting a subset of Lat/Lon indices.
      Given lat and lon arrays latitudes and longitudes,
      we want to determine the arrays inxdex_lat and
      index_lon of indices where the latitudes and longitudes
      fall in provided ranges ([minLat,maxLat] and [minLon,maxLon])

            lat[:]>=minLat and lat[:]<=maxLat
            lon[:]>=minLon and lon[:]<=maxLon
    """

    indexlat = np.nonzero((lat[:] >= lat_range[0]) & (lat[:] <= lat_range[-1]))[0]
    indexlon = np.nonzero((lon[:] >= lon_range[0]) & (lon[:] <= lon_range[-1]))[0]
    return indexlat, indexlon

# Small number
EPS = 1.0e-5

## Usage: <Name of variable in observed climatology>
## <Name of variable in reforecast climatology
## (same as the name in target forecast> <forecast model number>
print("In Python Script")
CMDARGS = str(sys.argv)
OBS_VAR = str(sys.argv[1])
FCST_VAR = str(sys.argv[2])
BC_VAR = str(sys.argv[3])
## This is used to figure out if the variable is a precipitation
## variable or not
UNIT = str(sys.argv[4])
LAT1, LAT2, LON1, LON2 = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
INIT_FCST_MON = int(sys.argv[9])

# Forecast model and ensemble input arguments:
MODEL_NAME = str(sys.argv[10])
LEAD_FINAL = int(sys.argv[11])
ENS_NUMC = int(sys.argv[12])
ENS_NUMF = int(sys.argv[13])

print(LEAD_FINAL)
print(ENS_NUMC)
print(ENS_NUMF)

FCST_SYR = int(sys.argv[14])
TARGET_FCST_SYR = int(sys.argv[14])
TARGET_FCST_EYR = int(sys.argv[15])
CLIM_SYR = int(sys.argv[16])
CLIM_EYR = int(sys.argv[17])

# Directory and file addresses
FCST_CLIM_INDIR = str(sys.argv[18])
OBS_CLIM_INDIR = str(sys.argv[19])
FCST_INDIR = str(sys.argv[20])

# Observation climatology filename templates:
OBS_CLIM_FILE_TEMPLATE = '{}/{}_obs_clim.nc'
FCST_CLIM_FILE_TEMPLATE = '{}/{}/{}_fcst_clim.nc'
MONTH_NAME_TEMPLATE = '{}01'
# GEOS5 filename TEMPLATE:
FCST_INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.nmme.monthly.{:04d}{:02d}.nc'

# Input mask
MASK_FILE = str(sys.argv[21])
MASK = read_nc_files(MASK_FILE, 'mask')[0, ]
LATS = read_nc_files(MASK_FILE, 'lat')
LONS = read_nc_files(MASK_FILE, 'lon')

### Output directory
OUTFILE_TEMPLATE = '{}/{}.{}.{}_{:04d}_{:04d}.nc'
OUTDIR = str(sys.argv[22])
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

ENSS = int(sys.argv[23])
ENSF = int(sys.argv[24])


print(f"Ensemble number is {ENS_NUMF}")
NUM_YRS = (CLIM_EYR-CLIM_SYR)+1
TINY = ((1/(NUM_YRS))/ENS_NUMF)/2
# Adjust quantile, if it is out of bounds
#          This value represents 1/NYRS/NENS/2, so about
#          half the prob. interval beyond the lowest value
#          (arbitrary choice) */
## This is probably used for real-time forecasts when a
## forecasted value happened to be an outlier of the
## reforecast climatology


##### Starting bias-correction from here

# First read observed climatology for the given variable
OBS_CLIM_FILE = OBS_CLIM_FILE_TEMPLATE.format(OBS_CLIM_INDIR, OBS_VAR)
print(f"Reading observed climatology {OBS_CLIM_FILE}")
OBS_CLIM_ARRAY = xr.open_dataset(OBS_CLIM_FILE)

# Then for forecast files:
for MON in [INIT_FCST_MON]:
    MONTH_NAME = MONTH_NAME_TEMPLATE.format((calendar.month_abbr[MON]).lower())
    ## This provides abbrevated version of the name of a month:
    ## (e.g. for January (i.e. Month number = 1) it will return "Jan").
    ## The abbrevated name is used in the forecasts file name
    print(f"Forecast Initialization month is {MONTH_NAME}")
    #First read forecast climatology for the given variable and forecast
    #initialzation month
    FCST_CLIM_INFILE = FCST_CLIM_FILE_TEMPLATE.format(FCST_CLIM_INDIR, \
    MODEL_NAME, FCST_VAR)
    print(f"Reading forecast climatology {FCST_CLIM_INFILE}")
    FCST_CLIM_ARRAY = xr.open_dataset(FCST_CLIM_INFILE)
    #First read raw forecasts
    FCST_COARSE = np.empty(((TARGET_FCST_EYR-TARGET_FCST_SYR)+1, \
    LEAD_FINAL, ENS_NUMF, len(LATS), len(LONS)))
    for LEAD_NUM in range(0, LEAD_FINAL): ## Loop from lead =0 to Final Lead
        for ens in range(ENS_NUMF):
            ens1 = ens+ENSS
            for INIT_FCST_YEAR in range(TARGET_FCST_SYR, TARGET_FCST_EYR+1):
                ## Reading forecast file
                FCST_DATE = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1) + \
                relativedelta(months=LEAD_NUM)
                FCST_YEAR, FCST_MONTH = FCST_DATE.year, FCST_DATE.month
                INFILE = FCST_INFILE_TEMPLATE.format(FCST_INDIR, \
                INIT_FCST_YEAR, ens1, MONTH_NAME, FCST_YEAR, FCST_MONTH)
                print(INFILE)
                FCST_COARSE[INIT_FCST_YEAR-TARGET_FCST_SYR, LEAD_NUM, ens, ] = \
                read_nc_files(INFILE, FCST_VAR)
    LAT_RANGE = [LAT1, LAT2]
    LON_RANGE = [LON1, LON2]
    indexLat, indexLon = slice_latlon(LATS, LONS, LAT_RANGE, LON_RANGE)

    ilat_min, ilat_max = indexLat[0], indexLat[-1]
    ilon_min, ilon_max = indexLon[0], indexLon[-1]

    nlats = len(LATS)
    nlons = len(LONS)
    #print("indexLat=",indexLat)
    #xprint("indexLon=",indexLon)
    print("LAT_RANGE=", LAT_RANGE, "LON_RANGE=", LON_RANGE)
    print("latmin=", ilat_min, "latmax=", ilat_max)
    print("lonmin=", ilon_min, "lonmax=", ilon_max)
    #exit()


    # Get the values (Numpy array) for the lat/lon ranges
    np_OBS_CLIM_ARRAY = OBS_CLIM_ARRAY.clim.sel(longitude=slice(LON1, LON2), \
    latitude=slice(LAT1, LAT2)).values
    np_FCST_CLIM_ARRAY = FCST_CLIM_ARRAY.clim.sel(longitude=slice(LON1, LON2), \
    latitude=slice(LAT1, LAT2)).values

    print("Latitude:  ", nlats, ilat_min, ilat_max)
    print("Longitude: ", nlons, ilon_min, ilon_max)
    print("np_OBS_CLIM_ARRAY:", np_OBS_CLIM_ARRAY.shape, \
    type(np_OBS_CLIM_ARRAY))
    print("np_FCST_CLIM_ARRAY:", np_FCST_CLIM_ARRAY.shape, \
    type(np_FCST_CLIM_ARRAY))

    CORRECT_FCST_COARSE = BCSD_function.latlon_calculations(ilat_min, \
    ilat_max, ilon_min, ilon_max, nlats, nlons, np_OBS_CLIM_ARRAY, \
    np_FCST_CLIM_ARRAY, LEAD_FINAL, TARGET_FCST_EYR, TARGET_FCST_SYR, \
    FCST_SYR, ENS_NUMF, MON, MONTH_NAME, BC_VAR, TINY, FCST_COARSE)

    CORRECT_FCST_COARSE = np.ma.masked_array(CORRECT_FCST_COARSE, \
    mask=CORRECT_FCST_COARSE == -999)
    OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, FCST_VAR, MODEL_NAME, \
    MONTH_NAME, TARGET_FCST_SYR, TARGET_FCST_EYR)
    print(f"Now writing {OUTFILE}")
    SDATE = datetime(TARGET_FCST_SYR, MON, 1)
    dates = [SDATE+relativedelta(years=n) for n in range(CORRECT_FCST_COARSE.shape[0])]
    write_4d_netcdf(OUTFILE, CORRECT_FCST_COARSE, FCST_VAR, MODEL_NAME, \
    'Bias corrected', UNIT, 5, LONS, LATS, ENS_NUMF, LEAD_FINAL, SDATE, dates)
