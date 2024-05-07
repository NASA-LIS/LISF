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
#  This module reorganizes NMME preciptation forecasts
#  Date: May 06, 2021
# Updated: Sarith Mahanama
#  Removed basemap call and added xarray and xesmf
#  module calls
#  Date: Nov 07, 2022
"""

from datetime import datetime
import os
import sys
from time import ctime as t_ctime
from time import time as t_time
import numpy as np
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
# pylint: enable=no-name-in-module
import xarray as xr
import xesmf as xe
import yaml
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
CYR = int(sys.argv[2]) ##
NMME_DOWNLOAD_DIR = str(sys.argv[3])
NMME_OUTPUT_DIR = str(sys.argv[4])
SUPPLEMENTARY_DIR = str(sys.argv[5])
NMME_MODEL = str(sys.argv[6])
ENS_NUM = int(sys.argv[7])
CONFIGFILE = str(sys.argv[8])
with open(CONFIGFILE, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)
LEAD_MONS = config['EXP']['lead_months']

MODEL = ['NCEP-CFSv2', 'NASA-GEOSS2S', 'CanSIPS-IC3', 'COLA-RSMAS-CCSM4', 'GFDL-SPEAR']
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', \
       'Sep', 'Oct', 'Nov', 'Dec']
MONTH = ['jan01', 'feb01', 'mar01', 'apr01', 'may01', 'jun01', 'jul01', \
          'aug01', 'sep01', 'oct01', 'nov01', 'dec01']
MONTHN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

print(MODEL[0])

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

INFILE_TEMP = '{}/{}/prec.{}.mon_{}.{:04d}.nc'
OUTDIR_TEMPLATE = '{}/{}/{}/{:04d}/ens{}/'
OUTFILE_TEMPLATE = '{}/{}.nmme.monthly.{:04d}{:02d}.nc'
if not os.path.exists(NMME_OUTPUT_DIR):
    os.makedirs(NMME_OUTPUT_DIR, exist_ok=True)

## Read in example fine spatial resolution file for lat and lon over domain
LATS, LONS = get_domain_info(CONFIGFILE, coord=True)

## Read in example coarse spatial resolution file for lat and lon over domain
EX_NMME_FILENAME = '/ex_raw_nmme_download.nc'
GE1 = SUPPLEMENTARY_DIR + EX_NMME_FILENAME
LONI = read_nc_files(GE1, 'X')
LATI = read_nc_files(GE1, 'Y')

## Read all forecast files
MM = CMN-1
XPRECI = np.empty([1, LATS.size, LONS.size])

if NMME_MODEL == 'CFSv2':
    MODEL = 'NCEP-CFSv2'
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], CYR)
    x = read_nc_files(INFILE, 'prec')
    XPREC = x[:, 0:LEAD_MONS, 0:ENS_NUM, :, :]
elif NMME_MODEL == 'GEOSv2':
    MODEL = 'NASA-GEOSS2S'
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], CYR)
    x = read_nc_files(INFILE, 'prec')
    XPREC = x[:, 0:LEAD_MONS, 0:ENS_NUM, :, :]
elif NMME_MODEL == 'CCM4':
    MODEL = 'CanSIPS-IC3'
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], CYR)
    x = read_nc_files(INFILE, 'prec')
    x1 = np.moveaxis(x, 1, 2)
    XPREC = x1[:,0:LEAD_MONS,10:20,:,:]
elif NMME_MODEL == 'GNEMO5':
    MODEL = 'CanSIPS-IC3'
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], CYR)
    x = read_nc_files(INFILE, 'prec')
    x1 = np.moveaxis(x, 1, 2)
    XPREC = x1[:,0:LEAD_MONS,0:10,:,:]
elif NMME_MODEL == 'CCSM4':
    MODEL = 'COLA-RSMAS-CCSM4'
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], CYR)
    x = read_nc_files(INFILE, 'prec')
    XPREC = x[:, 0:LEAD_MONS, 0:ENS_NUM, :, :]
elif NMME_MODEL == 'GFDL':
    MODEL = 'GFDL-SPEAR'
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], CYR)
    x = read_nc_files(INFILE, 'prec')
    XPREC = x[:, 0:LEAD_MONS, 0:ENS_NUM, :, :]
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
    data = np.array (XPREC[0,0:9,:,:,:]),
    dims=["mon","ens", "lat", "lon"],
    coords=dict(
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
mask_exp = mask_2d[np.newaxis, np.newaxis,:,:]
darray = np.array(ds_out_unmasked['XPREC'].values)
mask = np.broadcast_to(mask_exp, darray.shape)
darray[mask == 0] = -9999.
ds_out['XPREC'].values = darray

## Reorganize and write
YR = CYR
for m in range(0, ENS_NUM):
    for l in range(0, LEAD_MONS):
        XPRECI[0, :, :] = ds_out["XPREC"].values[l,m,:,:]
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
