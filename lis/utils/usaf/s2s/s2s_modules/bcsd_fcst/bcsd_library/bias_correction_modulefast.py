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
import calendar
import sys
from datetime import datetime
import os
from dateutil.relativedelta import relativedelta
import xarray as xr
import numpy as np
from Shrad_modules import read_nc_files
from BCSD_stats_functions import *
from BCSD_function import CALC_BCSD, get_index
import sys

## Usage: <Name of variable in observed climatology>
## <Name of variable in reforecast climatology
## (same as the name in target forecast> <forecast model number>
CMDARGS = str(sys.argv)
OBS_VAR = str(sys.argv[1])
FCST_VAR = str(sys.argv[2])
BC_VAR = str(sys.argv[3])
## This is used to figure out if the variable is a precipitation variable or not
UNIT = str(sys.argv[4])
LAT1, LAT2, LON1, LON2 = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
INIT_FCST_MON = int(sys.argv[9])

# Forecast model and ensemble input arguments:
LEAD_FINAL = int(sys.argv[10])
ENS_NUM = int(sys.argv[11])

print(LEAD_FINAL)
print(ENS_NUM)

FCST_SYR = int(sys.argv[12])
TARGET_FCST_SYR = int(sys.argv[12])
TARGET_FCST_EYR = int(sys.argv[13])
CLIM_SYR = int(sys.argv[14])
CLIM_EYR = int(sys.argv[15])

# Directory and file addresses
OBS_INDIR = str(sys.argv[16])
FCST_INDIR = str(sys.argv[17])

# Observation climatology filename templates:
OBS_CLIM_FILE_TEMPLATE = '{}/raw/Climatology/{}_obs_clim.nc'
FCST_CLIM_FILE_TEMPLATE = '{}/raw/Climatology/{}/{}_fcst_clim.nc'
MONTH_NAME_TEMPLATE = '{}01'
# GEOS5 filename template:
FCST_INFILE_TEMPLATE = '{}/raw/Monthly/{}/{:04d}/ens{:01d}/{}.cfsv2.{:04d}{:02d}.nc'

# Input mask
MASK_FILE = str(sys.argv[18])
#MASK = read_nc_files(MASK_FILE, 'mask')[0, ]
LATS = read_nc_files(MASK_FILE, 'lat')
LONS = read_nc_files(MASK_FILE, 'lon')

### Output directory
OUTFILE_TEMPLATE = '{}/{}.CFSv2.{}_{:04d}_{:04d}.nc'
OUTDIR = str(sys.argv[19])
print(OUTDIR)
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR, exist_ok=True)

print(f"Ensemble number is {ENS_NUM}")
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
    print(f"Reading forecast climatology {fcst_clim_infile}")
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
                print(infile)
                fcst_coarse[init_fcst_year-TARGET_FCST_SYR, lead_num, ens, ] = \
                read_nc_files(infile, FCST_VAR)[0,]
    # Defining array to store bias-corrected monthly forecasts
    correct_fcst_coarse = np.ones(((TARGET_FCST_EYR-TARGET_FCST_SYR)+1, \
    LEAD_FINAL, ENS_NUM, len(LATS), len(LONS)))*-999
    print("shape of fcst_coarse: ", fcst_coarse.shape)

    # Get the lat/lon indexes for the ranges
    ilat_min = get_index(LATS, LAT1)
    ilat_max = get_index(LATS, LAT2)
    ilon_min = get_index(LONS, LON1)
    ilon_max = get_index(LONS, LON2)
    nlats = len(LATS)
    nlons = len(LONS)
    fcst_clim_array2 = fcst_clim_array.rename_dims({"DIST": "DISTF", "time" : "timef" }) 
    # Get the values (Numpy array) for the lat/lon ranges
 
    np_obs_clim_xr = OBS_CLIM_ARRAY.clim.sel(longitude=slice(LON1, LON2), \
                        latitude=slice(LAT1, LAT2))
    np_fcst_clim_xr = fcst_clim_array2.clim.sel(longitude=slice(LON1, LON2), \
                         latitude=slice(LAT1, LAT2))

    print("Latitude:  ", nlats, ilat_min, ilat_max)
    print("Longitude: ", nlons, ilon_min, ilon_max)

    fcst_coarse_xr =  xr.DataArray(fcst_coarse, coords=
                                   {'fyrs': np.arange(TARGET_FCST_EYR-TARGET_FCST_SYR+1) ,'lead': np.arange(LEAD_FINAL),'ens': np.arange(ENS_NUM), 'latitude': np_fcst_clim_xr.latitude.values,'longitude': np_fcst_clim_xr.longitude.values},
                                   dims=["fyrs", "lead", "ens", "latitude", "longitude"])
    
    if (not np.array_equal(np_obs_clim_xr.latitude.values, np_fcst_clim_xr.latitude.values)) or (not np.array_equal(np_obs_clim_xr.longitude.values, np_fcst_clim_xr.longitude.values)):
        np_obs_clim_xr = np_obs_clim_xr.assign_coords({"longitude": np_fcst_clim_xr.longitude.values, "latitude": np_fcst_clim_xr.latitude.values})

    correct = xr.apply_ufunc(
        CALC_BCSD,
        np_obs_clim_xr.chunk({"latitude": "auto", "longitude": "auto"}).compute(),
        np_fcst_clim_xr.chunk({"latitude": "auto", "longitude": "auto"}).compute(),
        fcst_coarse_xr.chunk({"latitude": "auto", "longitude": "auto"}).compute(),
        input_core_dims=[['DIST','time'],['DISTF','timef'],['fyrs', 'lead', 'ens']],
        exclude_dims=set(('DIST','time','DISTF','timef','fyrs', 'lead', 'ens')),
        output_core_dims=[['fyrs', 'lead', 'ens']],
        vectorize=True,  # loop over non-core dims
        dask="forbidden",
        output_dtypes=[np.float64],
        kwargs={'MON': mon, 'MONTH_NAME': month_name, 'LEAD_FINAL' : LEAD_FINAL, 'TARGET_FCST_SYR' : TARGET_FCST_SYR, 
                'TARGET_FCST_EYR' : TARGET_FCST_EYR, 'FCST_SYR' : FCST_SYR,  'ENS_NUM' : ENS_NUM, 'BC_VAR' : BC_VAR, 'TINY' : TINY})

    correct_fcst_coarse = np.moveaxis(correct.values,[0,1],[-2,-1])
    correct_fcst_coarse = np.ma.masked_array(correct_fcst_coarse, \
                          mask=correct_fcst_coarse == -999)
    outfile = OUTFILE_TEMPLATE.format(OUTDIR, FCST_VAR, month_name, \
              TARGET_FCST_SYR, TARGET_FCST_EYR)
    print(f"Now writing {outfile}")
    sdate = datetime(TARGET_FCST_SYR, mon, 1)
    dates = [sdate+relativedelta(years=n) for n in range(correct_fcst_coarse.shape[0])]
    write_4d_netcdf(outfile, correct_fcst_coarse, FCST_VAR, 'CFSv2', \
    'Bias corrected', UNIT, 5, LONS, LATS, ENS_NUM, LEAD_FINAL, sdate, dates)

# Then for forecast files:
for MON in [INIT_FCST_MON]:
    monthly_calculations(MON)
