#!/usr/bin/env python

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
# Author: Abheera Hazra
#This module reorganizes
#NMME preciptation forecasts
#Date: May 06, 2021
# Update: KR Arsenault; Feb-2024; CFSv2 Jan/Feb update
"""

from datetime import datetime
import os
import sys
from time import ctime as t_ctime
from time import time as t_time
import numpy as np
import xarray as xr
import xesmf as xe
import yaml
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
# pylint: enable=no-name-in-module
# pylint: disable=import-error
from shrad_modules import read_nc_files
from bcsd_stats_functions import get_domain_info
from bcsd_function import VarLimits as lim
# pylint: enable=import-error

limits = lim()

def write_3d_netcdf(infile, var, varname, description, source, \
                    var_units, lons, lats, sdate):
    """write netcdf files"""
    rootgrp = nc4_dataset(infile, 'w', format='NETCDF4')
    longitude = rootgrp.createDimension('lon', len(lons))
    latitude = rootgrp.createDimension('lat', len(lats))
    time = rootgrp.createDimension('time', None)
    longitudes = rootgrp.createVariable('lon', 'f4', ('lon',))
    latitudes = rootgrp.createVariable('lat', 'f4', ('lat',))
    times = rootgrp.createVariable('time', 'f8', ('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname, 'f4', \
                                     ('time', 'lat', 'lon'), \
                                     fill_value=-9999., zlib=True)
    #import time
    rootgrp.description = description
    rootgrp.history = 'Created ' + t_ctime(t_time())
    rootgrp.source = source
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    varname.units = var_units
    string_date = datetime.strftime(sdate, "%Y-%m-%d")
    times.units = 'days since ' + string_date
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    varname[:, :, :] = var
    times[:] = nc4_date2num(sdate, units=times.units, calendar=times.calendar)
    rootgrp.close()

CMDARGS = str(sys.argv)
CMN = int(sys.argv[1])  ##
NMME_OUTPUT_DIR = str(sys.argv[2])
NMME_MODEL = str(sys.argv[3])
CONFIGFILE = str(sys.argv[4])
# Load config file
with open(CONFIGFILE, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)
NMME_DOWNLOAD_DIR = config['BCSD']['nmme_download_dir']
SUPPLEMENTARY_DIR = config['SETUP']['supplementarydir'] + '/bcsd_fcst/'
ensemble_sizes = config['EXP']['ensemble_sizes'][0]
ENS_NUM = ensemble_sizes[NMME_MODEL]
LEAD_MON = config['EXP']['lead_months']

MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', \
       'Sep', 'Oct', 'Nov', 'Dec']
MONTH = ['jan01', 'feb01', 'mar01', 'apr01', 'may01', 'jun01', 'jul01', \
          'aug01', 'sep01', 'oct01', 'nov01', 'dec01']
MONTHN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

LEADS1 = np.zeros((12, 12), dtype=int)
LDYR = np.zeros((12, 12), dtype=int)

##Fill LEADS1 and LDYR for calling the correct lead month for a given year
for i in range(0, 12):
    for j in range(0, 12):
        k = MONTHN[i]+j
        KY = 0
        if k >= 13:
            k = k-12
            KY = 1
        LEADS1[i, j] = k
        LDYR[i, j] = KY

INFILE_TEMP = '{}/{}/prec.{}.mon_{}_{:04d}_{:04d}.nc'
OUTDIR_TEMPLATE = '{}/{}/{}/{:04d}/ens{}/'
OUTFILE_TEMPLATE = '{}/{}.nmme.monthly.{:04d}{:02d}.nc'
if not os.path.exists(NMME_OUTPUT_DIR):
    os.makedirs(NMME_OUTPUT_DIR)

## Read in example fine spatial resolution file for lat and lon over domain
LATS, LONS = get_domain_info(CONFIGFILE, coord=True)

## Read in example coarse spatial resolution file for lat and lon over domain
EX_NMME_FILENAME = '/ex_raw_nmme_download.nc'
GE1 = SUPPLEMENTARY_DIR + EX_NMME_FILENAME
LONI = read_nc_files(GE1, 'X')
LATI = read_nc_files(GE1, 'Y')

## Read all forecast files
MM = CMN-1
XPREC = np.empty([40, LEAD_MON, ENS_NUM, 181, 360])
XPRECI = np.empty([1, LATS.size, LONS.size])

if NMME_MODEL == 'CFSv2':
    MODEL = 'NCEP-CFSv2'
    SYR1 = 1982
    EYR1 = 2010
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
    XPREC[0:29,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]

    if MON[MM] == 'Jan' or MON[MM] == 'Feb':
        SYR2 = 2011
        EYR2 = 2011
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
        XPREC[29,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]

        SYR3 = 2011
        EYR3 = 2021
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR3, EYR3)
        XPREC[30:40,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]

    else:
        SYR2 = 2011
        EYR2 = 2021
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
        XPREC[29:40,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]

elif NMME_MODEL == 'GEOSv2':
    MODEL = 'NASA-GEOSS2S'
    if MM == 0:
        SYR1 = 1982
        EYR1 = 2017
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
        XPREC[0:36,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]

        SYR2 = 2018
        EYR2 = 2021
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
        XPREC[36:40,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]
    else:
        SYR1 = 1982
        EYR1 = 2016
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
        XPREC[0:35,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]

        SYR2 = 2017
        EYR2 = 2021
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
        XPREC[35:40,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]
elif NMME_MODEL == 'CCM4':
    MODEL = 'CanSIPS-IC3'
    SYR1 = 1991
    EYR1 = 2020
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
    x = read_nc_files(INFILE, 'prec')
    x1 = np.moveaxis(x, 1, 2)
    XPREC[9:39,:,:,:,:] = x1[:,0:LEAD_MON,10:20,:,:]

    SYR2 = 2021
    EYR2 = 2021
    #INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
    #x = read_nc_files(INFILE, 'prec')
    #x1 = np.moveaxis(x, 1, 2)
    #XPREC[39:40,:,:,:,:] = x1[:,0:LEAD_MON,10:20,:,:]
elif NMME_MODEL == 'GNEMO5':
    MODEL = 'CanSIPS-IC3'
    SYR1 = 1991
    EYR1 = 2020
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
    x = read_nc_files(INFILE, 'prec')
    x1 = np.moveaxis(x, 1, 2)
    XPREC[9:39,:,:,:,:] = x1[:,0:LEAD_MON,0:10,:,:]

    SYR2 = 2021
    EYR2 = 2021
    #INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
    #x = read_nc_files(INFILE, 'prec')
    #x1 = np.moveaxis(x, 1, 2)
    #XPREC[39:40,:,:,:,:] = x1[:,0:LEAD_MON,0:10,:,:]
elif NMME_MODEL == 'CCSM4':
    MODEL = 'COLA-RSMAS-CCSM4'
    SYR1 = 1982
    EYR1 = 2021
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
    XPREC = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]
elif NMME_MODEL == 'GFDL':
    MODEL = 'GFDL-SPEAR'
    SYR1 = 1991
    EYR1 = 2020
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
    # IH = nc.Dataset(INFILE, mode='r')
    XPREC[9:39,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]

    SYR2 = 2021
    EYR2 = 2021
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
    # IH = nc.Dataset(INFILE, mode='r')
    XPREC[39:40,:,:,:,:] = read_nc_files(INFILE, 'prec')[:, 0:LEAD_MON, 0:ENS_NUM, :, :]
else:
    print(f"[ERR] Invalid argument for NMME_MODEL! Received {NMME_MODEL}")
    sys.exit(1)

print('Reading:', INFILE)

## Convert units to mm/s or kg/m2/s
XPREC = XPREC/86400
ds_in = xr.Dataset(
            {
                "lat": (["lat"], LATI),
                "lon": (["lon"], LONI),
        })

ds_in["slice"] = xr.DataArray(
    data = np.array (XPREC[0,0,0,:,:]),
    dims=["lat", "lon"],
    coords=dict(
        lat=(["lat"], LATI),
        lon=(["lon"], LONI))
    )

ds_in["XPREC"] = xr.DataArray(
    data = np.array (XPREC[:,0:9,:,:,:]),
    dims=["year", "mon","ens", "lat", "lon"],
    coords=dict(
        year=(["year"], np.arange(40)),
        mon=(["mon"], np.arange(9)),
        ens=(["ens"], np.arange(ENS_NUM)),
        lat=(["lat"], LATI),
        lon=(["lon"], LONI))
    )
ds_out_unmasked = xr.Dataset(
            {
                "lat": (["lat"], LATS),
                "lon": (["lon"], LONS),
        })
regridder = xe.Regridder(ds_in, ds_out_unmasked, "conservative", periodic=True)
ds_out_unmasked = regridder(ds_in)
ds_out = ds_out_unmasked.copy()

# LDT mask
ldt_xr = xr.open_dataset(config['SETUP']['supplementarydir'] + '/lis_darun/' + \
        config['SETUP']['ldtinputfile'])
mask_2d = np.array(ldt_xr['LANDMASK'].values)
mask_exp = mask_2d[np.newaxis, np.newaxis, np.newaxis,:,:]
darray = np.array(ds_out_unmasked['XPREC'].values)
mask = np.broadcast_to(mask_exp, darray.shape)
darray[mask == 0] = -9999.
ds_out['XPREC'].values = darray

YR = 1981
print(XPREC.shape)
## Reorganize and write
for y in range(0, 40):
    YR = YR+1
    for m in range(0, ENS_NUM):
        for l in range(0, 9):
            XPRECI[0, :, :] = ds_out["XPREC"].values[y,l,m,:,:]
            jy = YR+LDYR[MM, l]
            l1 = LEADS1[MM, l]
            print('Year:', jy, ',leads:', l1)
            OUTDIR = OUTDIR_TEMPLATE.format(NMME_OUTPUT_DIR, MONTH[MM], NMME_MODEL, YR, m+1)
            OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, MONTH[MM], jy, l1)
            if not os.path.exists(OUTDIR):
                os.makedirs(OUTDIR)

            XPRECI = np.nan_to_num(XPRECI, nan=-9999.)
            XPRECI = limits.clip_array(XPRECI, var_name="PRECTOT", max_val=0.004, precip=True)
            LATS = np.nan_to_num(LATS, nan=-9999.)
            LONS = np.nan_to_num(LONS, nan=-9999.)

            SDATE = datetime(YR, MM+1, 1)
            write_3d_netcdf(OUTFILE, XPRECI, 'PRECTOT', 'Downscaled to 0.25deg', \
            'Raw NMME at 1deg', 'kg m-2 s-1', LONS, LATS, SDATE)
            print(f"Writing {OUTFILE}")
