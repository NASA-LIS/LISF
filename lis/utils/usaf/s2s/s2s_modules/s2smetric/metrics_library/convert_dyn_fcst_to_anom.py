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
import subprocess
import sys
import yaml

# Third-party modules
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr

# Local modules
# pylint: disable=import-error
from metricslib import sel_var, compute_anomaly
# pylint: enable=import-error
# pylint: disable=consider-using-f-string
#
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
BASEDIR = BASEOUTDIR + "/DYN_ANOM/"
METRIC_VARS = CONFIG["POST"]["metric_vars"]
HINDCASTS = CONFIG["SETUP"]["E2ESDIR"] + '/hindcast/s2spost/' + '{:02d}/'.format(FCST_INIT_MON)
FORECASTS = "./s2spost/"
CURRENTDATE = date(TARGET_YEAR, FCST_INIT_MON, 1)

FCST_INIT_DAY = 1
OUTDIR = BASEDIR + '/' + HYD_MODEL
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR, exist_ok=True)

OUTFILE_TEMPLATE = '{}/{}_{}_ANOM_init_monthly_{:02d}_{:04d}.nc'
# name of variable, forecast initial month and forecast year is in the file

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
    'PA.ALL_DD.{:04d}{:02d}01_DT.0000_FP.{:04d}{:02d}??-{:04d}{:02d}??_DF.NC'

CLIM_INFILE_TEMPLATE2 = 'PA.ALL_DD.*{:02d}01_DT.0000_FP.*{:02d}??-*{:02d}??_DF.NC'

TARGET_INFILE_TEMPLATE = TARGET_INFILE_TEMPLATE1 + TARGET_INFILE_TEMPLATE2
CLIM_INFILE_TEMPLATE = CLIM_INFILE_TEMPLATE1 + CLIM_INFILE_TEMPLATE2
## String in this format allows to select all the files for the given month

for var_name in METRIC_VARS:
    for lead in range(LEAD_NUM):
        print('[INFO] Reading output from Hindcast runs')
        print(f'[INFO] var_name, lead: {var_name} {lead}')

        ## Step-1: Read and process the climatology
        smon = datetime(CLIM_SYR, FCST_INIT_MON, FCST_INIT_DAY) + \
                    relativedelta(months=lead)
        emon = datetime(CLIM_SYR, FCST_INIT_MON, FCST_INIT_DAY) + \
                    relativedelta(months=lead+1)
        ## Adding 1 to lead to make sure the file read is from the month after

        INFILE = CLIM_INFILE_TEMPLATE.format(HINDCASTS, \
                                             FCST_INIT_MON, NMME_MODEL, \
                                             NMME_MODEL.upper(), \
                                             FCST_INIT_MON, \
                                             smon.month, emon.month)

        print(f"[INFO] Reading forecast climatology {INFILE}")
        infile1 = glob.glob(INFILE)

        # First reading all available years for the given
        # forecast initialization month
        all_clim_data1 = xr.open_mfdataset(infile1, combine='by_coords')

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
        smon1 = datetime(TARGET_YEAR, FCST_INIT_MON, FCST_INIT_DAY) + \
            relativedelta(months=lead)
        emon1 = datetime(TARGET_YEAR, FCST_INIT_MON, FCST_INIT_DAY) + \
            relativedelta(months=lead+1)

        INFILE = TARGET_INFILE_TEMPLATE.format(FORECASTS, \
                                               TARGET_YEAR, FCST_INIT_MON, NMME_MODEL, \
                                               NMME_MODEL.upper(), \
                                               TARGET_YEAR, FCST_INIT_MON, \
                                               smon1.year, smon1.month, \
                                               emon1.year, emon1.month)

        print(f"[INFO] Reading target {INFILE}")

        # Note target will always have only one time step
        target_data = xr.open_mfdataset(INFILE, combine='by_coords')

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
            all_anom = np.ones((ens_count, LEAD_NUM, lat_count, lon_count))*-99
        print('[INFO] Converting data into anomaly')

        all_clim_mean = all_clim_data.mean (dim = ['time','ensemble'], skipna = True)
        if (not np.array_equal(all_clim_mean.lat.values, target_fcst_data.lat.values)) or \
           (not np.array_equal(all_clim_mean.lon.values, target_fcst_data.lon.values)):
            all_clim_mean = all_clim_mean.assign_coords({"lon": target_fcst_data.lon.values,
                                                         "lat": target_fcst_data.lat.values}) 
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

        for ens in range(ens_count):
            all_anom[ens, lead, :, :] = this_anom [:,:,ens,0]

        del all_clim_data, target_fcst_data, all_clim_mean

    ### Step-4 Writing output file
    all_anom = np.ma.masked_array(all_anom, mask=(all_anom == -99))

    ## Creating an latitude and longitude array based on locations of corners
    lats = np.arange(target_data.attrs['SOUTH_WEST_CORNER_LAT'], \
                     target_data.attrs['SOUTH_WEST_CORNER_LAT'] + \
                     (lat_count*0.25), 0.25)
    lons = np.arange(target_data.attrs['SOUTH_WEST_CORNER_LON'], \
                     target_data.attrs['SOUTH_WEST_CORNER_LON'] + \
                     (lon_count*0.25), 0.25)
    OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, NMME_MODEL, \
                                      var_name, FCST_INIT_MON, TARGET_YEAR)
    anom_xr = xr.Dataset()
    anom_xr['anom'] = (('ens', 'time', 'latitude', 'longitude'), all_anom)
    anom_xr.coords['latitude'] = (('latitude'), lats)
    anom_xr.coords['longitude'] = (('longitude'), lons)
    anom_xr.coords['time'] = (('time'), np.arange(0, LEAD_NUM, dtype=int))
    anom_xr.coords['ens'] = (('ens'), np.arange(0, ens_count, dtype=int))
    print(f"[INFO] Writing {OUTFILE}")
    anom_xr.to_netcdf(OUTFILE)

# Create file tag indicating completion
FILENAME = f"{OUTDIR}/anom.{HYD_MODEL}.{NMME_MODEL}.done"
CMD = f"touch {FILENAME}"
if subprocess.call(CMD, shell=True) != 0:
    print(f"[ERR] Cannot create {FILENAME}")
# sys.exit(1)
