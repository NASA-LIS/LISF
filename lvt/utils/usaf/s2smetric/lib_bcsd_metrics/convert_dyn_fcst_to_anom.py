#!/usr/bin/env python
"""
## Shrad Shukla March 2017
## Adapted and updated by Abheera Hazra 2019
## This script calculates anomaly of NMME forecasts
"""

# coding: utf-8

# In[1]:

from __future__ import division
import os
import sys
from datetime import datetime
import glob
import xarray as xr
from dateutil.relativedelta import relativedelta
import numpy as np
import matplotlib
from All_functions import Sel_var
matplotlib.use('Agg')
#from Shrad_modules import *

CMDARGS = str(sys.argv)
FCST_INIT_MON = int(sys.argv[1])
HYD_MODEL = sys.argv[2]
LEAD_NUM = int(sys.argv[3])
DOMAIN = str(sys.argv[4])
TARGET_YEAR = int(sys.argv[5])
MODEL = str(sys.argv[6])
FCST_INIT_DAY = 1
CSYR = int(sys.argv[7])
CEYR = int(sys.argv[8])

## Hardwired output directories where output anomaly data goes to
BASEDIR1 = '/discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/'
BASEDIR2 = 'NMME_FCST_DATA_AF/DYN_ANOM'
BASEDIR = BASEDIR1 + BASEDIR2
print(BASEDIR)

## It is assumed that the following directory is already created
OUTDIR = BASEDIR + '/' + DOMAIN + '/' + HYD_MODEL
OUTFILE_TEMPLATE = '{}/{}_{}_ANOM_init_monthly_{:02d}_{:04d}.nc'
# name of variable, forecast initial month and forecast year is in the file name

print(OUTDIR)
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

## Hardwired directory addresses
FAME_MODEL_RUN = '/discover/nobackup/projects/lis_aist17/emkemp/AFWA/lis74_s2s_cf/'

## Climatology years
CLIM_SYR, CLIM_EYR = CSYR, CEYR

TARGET_INFILE_TEMPLATE1 = '{}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_'
TARGET_INFILE_TEMPLATE2 = 'PA.LIS-S2S_DP.{:04d}{:02d}??-{:04d}{:02d}??_TP.0000-0000_DF.NC'
TARGET_INFILE_TEMPLATE = TARGET_INFILE_TEMPLATE1 + TARGET_INFILE_TEMPLATE2

CLIM_INFILE_TEMPLATE1 = '{}/{:02d}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_'
CLIM_INFILE_TEMPLATE2 = 'PA.LIS-S2S_DP.*{:02d}??-*{:02d}??_TP.0000-0000_DF.NC'
CLIM_INFILE_TEMPLATE = CLIM_INFILE_TEMPLATE1 + CLIM_INFILE_TEMPLATE2
## String in this format allows the select all the files for the given month

print('Now reading output from Hindcast runs')
for VAR_NAME in ['RootZone-SM', 'Surface-SM', 'Streamflow']:
    ## Counter for LEAD
    for LEAD in range(LEAD_NUM):
        ## Step-1: Read and process the climatology
        SMON = datetime(CLIM_SYR, FCST_INIT_MON, FCST_INIT_DAY) + \
                    relativedelta(months=LEAD)
        EMON = datetime(CLIM_SYR, FCST_INIT_MON, FCST_INIT_DAY) + \
                    relativedelta(months=LEAD+1)
        ## Adding 1 to lead to make sure the file read is from the month after
        INFILE = CLIM_INFILE_TEMPLATE.format(FAME_MODEL_RUN +'for_abheera', \
                 FCST_INIT_MON, MODEL, SMON.month, EMON.month)
        print("Reading forecast climatology {}".format(INFILE))
        test = glob.glob(INFILE)
        #print(test)
        INFILE1 = test
        print("Reading forecast climatology {}".format(INFILE1))
        # First reading all available years for the given
        # forecast initializatio month
        ALL_CLIM_DATA1 = xr.open_mfdataset(INFILE1, combine='by_coords')
        # Now selecting only the years that are within the climatology
        SEL_CIM_DATA = ALL_CLIM_DATA1.sel(time=\
                       (ALL_CLIM_DATA1.coords['time.year'] >= \
                       CLIM_SYR) & (ALL_CLIM_DATA1.coords['time.year'] <= \
                       CLIM_EYR))
        # Now Selecting the climatology further for the given variable
        #print(SEL_CIM_DATA)
        ALL_CLIM_DATA = Sel_var(SEL_CIM_DATA, VAR_NAME, HYD_MODEL)
        # ALL_CLIM_DATA has all the needed climatology for a given variable
        # To-DO: In future we may want to save the climatology all
        # together in a file so we don't have to read climatologies every time

        ####### Step-2: Read the target forecast which needs to be converted
        ## into anomaly
        SMON1 = datetime(TARGET_YEAR, FCST_INIT_MON, FCST_INIT_DAY) + \
        relativedelta(months=LEAD)
        EMON1 = datetime(TARGET_YEAR, FCST_INIT_MON, FCST_INIT_DAY) + \
        relativedelta(months=LEAD+1)
        #print("Syear=", SMON1.year, "Smonth=", SMON1.month)
        #print("Eyear=", EMON1.year, "Emonth=", EMON1.month)
        INFILE = TARGET_INFILE_TEMPLATE.format(FAME_MODEL_RUN, \
                 'nmme_forecast_proxy/monthly_files', \
                 MODEL, SMON1.year, SMON1.month, \
                 EMON1.year, EMON1.month)
        print("Reading Target {}".format(INFILE))
        # Note target will always have only one time step
        TARGET_DATA = xr.open_mfdataset(INFILE, combine='by_coords')
        ## Now selecting the desired variable
        TARGET_FCST_DATA = Sel_var(TARGET_DATA, VAR_NAME, HYD_MODEL)
        TARGET_FCST_DATA = TARGET_FCST_DATA.load()
        ALL_CLIM_DATA = ALL_CLIM_DATA.load()

        ## Step-3 loop through each grid cell and convert data into anomaly
        # Defining array to store anomaly data
        LAT_COUNT, LON_COUNT, ENS_COUNT = \
        len(TARGET_DATA.coords['lat']), \
        len(TARGET_DATA.coords['lon']), \
        len(TARGET_DATA.coords['ensemble'])

        ## Note that ENS_COUNT is coming from the Target forecasts,
        ## so if the target_forecasts have 4 members there will be
        ## 4 members in anomaly output and so on.
        if LEAD == 0:
            ALL_ANOM = np.ones((ENS_COUNT, LEAD_NUM, LAT_COUNT, LON_COUNT))*-99
        print('Now converting data into anomaly')
        for lat in range(LAT_COUNT):
            for lon in range(LON_COUNT):
                FCST_CLIM_TS = ALL_CLIM_DATA.isel(lat=lat, \
                               lon=lon).values.flatten()
                ## Note that this step combines all ensemble members
                ## to make one climatology
                FCST_CLIM = np.mean(FCST_CLIM_TS, axis=None)
                for ENS in range(ENS_COUNT):
                    TARGET_VAL = TARGET_FCST_DATA.isel(ensemble=ENS, \
                                 lat=lat, lon=lon).values
                    ALL_ANOM[ENS, LEAD, lat, lon] = TARGET_VAL-FCST_CLIM
        del ALL_CLIM_DATA, TARGET_FCST_DATA
    ### Step-4 Writing output file
    ALL_ANOM = np.ma.masked_array(ALL_ANOM, mask=ALL_ANOM == -99)
    ## Creating an latitude and longitude array based on locations of corners
    LATS = np.arange(TARGET_DATA.attrs['SOUTH_WEST_CORNER_LAT'], \
                     TARGET_DATA.attrs['SOUTH_WEST_CORNER_LAT'] + \
                     (LAT_COUNT*0.25), 0.25)
    LONS = np.arange(TARGET_DATA.attrs['SOUTH_WEST_CORNER_LON'], \
                     TARGET_DATA.attrs['SOUTH_WEST_CORNER_LON'] + \
                     (LON_COUNT*0.25), 0.25)
    OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, MODEL, \
                                      VAR_NAME, FCST_INIT_MON, TARGET_YEAR)
    ANOM_XR = xr.Dataset()
    ANOM_XR['ANOM'] = (('ens', 'lead', 'latitude', 'longitude'), ALL_ANOM)
    ANOM_XR.coords['latitude'] = (('latitude'), LATS)
    ANOM_XR.coords['longitude'] = (('longitude'), LONS)
    ANOM_XR.coords['lead'] = (('lead'), np.arange(0, LEAD_NUM, dtype=np.int))
    ANOM_XR.coords['ens'] = (('ens'), np.arange(0, ENS_COUNT, dtype=np.int))
    print("Writing {}".format(OUTFILE))
    ANOM_XR.to_netcdf(OUTFILE)
