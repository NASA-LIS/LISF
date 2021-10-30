#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: convert_dyn_fcst_to_sanom.py
#
# PURPOSE: Calculates standardized anomaly of NMME-forced LIS forecasts.
#
# REVISION HISTORY:
# ?? Mar 2017: Shrad Shukla/UCSB, first version.
# ?? ??? 2019: Abheera Hazra/UMD, second version.
# 22 Oct 2021: Eric Kemp/SSAI, updated for 557WW.
# 30 Oct 2021: Eric Kemp/SSAI, updated to use s2smetric config file.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
from datetime import datetime
import glob
import os
import subprocess
import sys

# Third-party modules
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr

# Local modules
from metricslib import sel_var

# Local constants. FIXME: Collect into common file for whole system

## hardwired output directories where output anomaly data goes to
#_BASEDIR = '/discover/nobackup/projects/lis_aist17/emkemp/AFWA'
#_BASEDIR += '/lis74_s2s_patches/work/POST/forecasts/DYN_ANOM'

## Hardwired directory addresses
#_FAME_MODEL_RUN = '/discover/nobackup/projects/lis_aist17/karsenau'
#_FAME_MODEL_RUN += '/E2ES_Test/29-Oct-2021/s2spost/forecasts'
#_HINDCASTS = '/discover/nobackup/projects/lis_aist17/karsenau'
#_HINDCASTS += '/E2ES_Test/29-Oct-2021/s2smetric/hindcasts'

#_FORECASTS = '/discover/nobackup/projects/lis_aist17/karsenau'
#_FORECASTS += '/E2ES_Test/29-Oct-2021/s2spost/forecasts'

# Start reading from command line.
fcst_init_mon = int(sys.argv[1])
target_year = int(sys.argv[2])
nmme_model = sys.argv[3]
configfile = sys.argv[4]

config = configparser.ConfigParser()
config.read(configfile)

hyd_model = config["s2smetric"]["lsm_model"]
lead_num = int(config["s2smetric"]["lead"])
domain_name = config["s2smetric"]["domain"]
clim_syr = int(config["s2smetric"]["csyr"])
clim_eyr = int(config["s2smetric"]["ceyr"])
basedir = config["s2smetric"]["sanom_base_dir"]
metric_vars = config["s2smetric"]["metric_vars"].split()
hindcasts = config["s2smetric"]["hindcasts"]
forecasts = config["s2smetric"]["forecasts"]

#hyd_model = sys.argv[2]
#lead_num = int(sys.argv[3])
#DOMAIN_NAME = str(sys.argv[4])
#NMME_MODEL = str(sys.argv[6])
#csyr = int(sys.argv[7])
#ceyr = int(sys.argv[8])

FCST_INIT_DAY = 1
outdir = basedir + '/' + domain_name + '/' + hyd_model
if not os.path.exists(outdir):
    os.makedirs(outdir)

OUTFILE_TEMPLATE = '{}/{}_{}_SANOM_init_monthly_{:02d}_{:04d}.nc'
# name of variable, forecast initial month and forecast year is in the file
# name


# Climatology years
#clim_syr, clim_eyr = csyr, ceyr

TARGET_INFILE_TEMPLATE1 = \
    '{}/{}/{:02d}/cf_{}_????{:02d}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_'
TARGET_INFILE_TEMPLATE2 = \
    'PA.LIS-S2S_DP.{:04d}{:02d}??-{:04d}{:02d}??_TP.0000-0000_DF.NC'
TARGET_INFILE_TEMPLATE = TARGET_INFILE_TEMPLATE1 + TARGET_INFILE_TEMPLATE2

CLIM_INFILE_TEMPLATE1 = \
    '{}/{}/{:02d}/cf_{}_????{:02d}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_'
CLIM_INFILE_TEMPLATE2 = 'PA.LIS-S2S_DP.*{:02d}??-*{:02d}??_TP.0000-0000_DF.NC'
CLIM_INFILE_TEMPLATE = CLIM_INFILE_TEMPLATE1 + CLIM_INFILE_TEMPLATE2
## String in this format allows the select all the files for the given month

#for var_name in ['RootZone-SM', 'Surface-SM', 'Streamflow']:
for var_name in metric_vars:
    for lead in range(lead_num):
        print('[INFO] Reading output from Hindcast runs')
        print(f'[INFO] var_name, lead: {var_name} {lead}')

        ## Step-1: Read and process the climatology
        smon = datetime(clim_syr, fcst_init_mon, FCST_INIT_DAY) + \
                    relativedelta(months=lead)
        emon = datetime(clim_syr, fcst_init_mon, FCST_INIT_DAY) + \
                    relativedelta(months=lead+1)
        ## Adding 1 to lead to make sure the file read is from the month after
        INFILE = CLIM_INFILE_TEMPLATE.format(hindcasts, \
                                             nmme_model, fcst_init_mon, \
                                             nmme_model.upper(), \
                                             smon.month, \
                                             nmme_model.upper(), \
                                             smon.month, emon.month)
        #print(f"[INFO] reading forecast climatology {INFILE}")
        infile1 = glob.glob(INFILE)
        #print(f"[INFO] reading forecast climatology {infile1}")
        print("[INFO] Reading forecast climatology")

        # First reading all available years for the given
        # forecast initialization month
        all_clim_data1 = xr.open_mfdataset(infile1, combine='by_coords')

        # Now selecting only the years that are within the climatology
        sel_cim_data = all_clim_data1.sel(time= \
                       (all_clim_data1.coords['time.year'] >= \
                       clim_syr) & (all_clim_data1.coords['time.year'] <= \
                       clim_eyr))

        # Now selecting the climatology further for the given variable
        all_clim_data = sel_var(sel_cim_data, var_name, hyd_model)

        # all_clim_data has all the needed climatology for a given variable
        # To-DO: In future we may want to save the climatology all
        # together in a file so we don't have to read climatologies every time

        ####### Step-2: Read the target forecast which needs to be converted
        ## into standardized anomaly
        smon1 = datetime(target_year, fcst_init_mon, FCST_INIT_DAY) + \
            relativedelta(months=lead)
        emon1 = datetime(target_year, fcst_init_mon, FCST_INIT_DAY) + \
            relativedelta(months=lead+1)

        INFILE = TARGET_INFILE_TEMPLATE.format(forecasts,
                                               nmme_model, fcst_init_mon,
                                               nmme_model.upper(),
                                               smon1.month,
                                               nmme_model.upper(),
                                               smon1.year, smon1.month, \
                                               emon1.year, emon1.month)

        print(f"[INFO] Reading target {INFILE}")

        # Note target will always have only one time step
        target_data = xr.open_mfdataset(INFILE, combine='by_coords')

        ## Now selecting the desired variable
        target_fcst_data = sel_var(target_data, var_name, hyd_model)
        target_fcst_data = target_fcst_data.load()
        all_clim_data = all_clim_data.load()

        ## Step-3 loop through each grid cell and convert data into
        # standardized anomaly
        # Defining array to store standardized anomaly data
        lat_count, lon_count, ens_count = \
            len(target_data.coords['lat']), \
            len(target_data.coords['lon']), \
            len(target_data.coords['ensemble'])

        ## Note that ens_count is coming from the Target forecasts,
        ## so if the target_forecasts have 4 members there will be
        ## 4 members in standardized anomaly output and so on.
        if lead == 0:
            all_anom = np.ones((ens_count, lead_num, lat_count, lon_count))*-99
        print('[INFO] Converting data into standardized anomaly')
        for lat in range(lat_count):
            for lon in range(lon_count):
                fcst_clim_ts = all_clim_data.isel(lat=lat, \
                                                  lon=lon).values.flatten()
                ## Note that this step combines all ensemble members
                ## to make one climatology
                fcst_clim = np.mean(fcst_clim_ts, axis=None)
                fcst_std = np.std(fcst_clim_ts, axis=None)
                for ens in range(ens_count):
                    target_val = target_fcst_data.isel(ensemble=ens, \
                                 lat=lat, lon=lon).values
                    all_anom[ens, lead, lat, lon] = \
                        (target_val - fcst_clim) / fcst_std
        del all_clim_data, target_fcst_data

    ### Step-4 Writing output file
    all_anom = np.ma.masked_array(all_anom, mask=(all_anom == -99))

    ## Creating an latitude and longitude array based on locations of corners
    lats = np.arange(target_data.attrs['SOUTH_WEST_CORNER_LAT'], \
                     target_data.attrs['SOUTH_WEST_CORNER_LAT'] + \
                     (lat_count*0.25), 0.25)
    lons = np.arange(target_data.attrs['SOUTH_WEST_CORNER_LON'], \
                     target_data.attrs['SOUTH_WEST_CORNER_LON'] + \
                     (lon_count*0.25), 0.25)
    OUTFILE = OUTFILE_TEMPLATE.format(outdir, nmme_model, \
                                      var_name, fcst_init_mon, target_year)
    anom_xr = xr.Dataset()
    anom_xr['anom'] = (('ens', 'lead', 'latitude', 'longitude'), all_anom)
    anom_xr.coords['latitude'] = (('latitude'), lats)
    anom_xr.coords['longitude'] = (('longitude'), lons)
    anom_xr.coords['lead'] = (('lead'), np.arange(0, lead_num, dtype=int))
    anom_xr.coords['ens'] = (('ens'), np.arange(0, ens_count, dtype=int))
    print(f"[INFO] Writing {OUTFILE}")
    anom_xr.to_netcdf(OUTFILE)

# Create file tag indicating completion
filename = f"{outdir}/sanom.{hyd_model}.{nmme_model}.done"
cmd = f"touch {filename}"
if subprocess.call(cmd, shell=True) != 0:
    print(f"[ERR] Cannot create {filename}")
    sys.exit(1)
