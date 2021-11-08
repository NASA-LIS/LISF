#!/usr/bin/env python
"""
# Author: Abheera Hazra
#This module reorganizes
#NMME preciptation forecasts
#Date: May 06, 2021
# In[28]:
"""
from __future__ import division
from datetime import datetime
import os
import sys
from time import ctime as t_ctime
from time import time as t_time
import numpy as np
from mpl_toolkits import basemap
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
# pylint: enable=no-name-in-module
from Shrad_modules import read_nc_files

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

MODEL = ['NCEP-CFSv2', 'NASA-GEOSS2S', 'CanCM4i', 'GEM-NEMO', \
         'COLA-RSMAS-CCSM4', 'GFDL-SPEAR']
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
OUTDIR_TEMPLATE = '{}/{}/{:04d}/ens{}/'
OUTFILE_TEMPLATE = '{}/{}.nmme.monthly.{:04d}{:02d}.nc'
if not os.path.exists(NMME_OUTPUT_DIR):
    os.makedirs(NMME_OUTPUT_DIR)

## Read in example fine spatial resolution file for lat and lon over AFRICOM
EX_FCST_FILENAME = '/ex_raw_fcst_download.nc'
GE = SUPPLEMENTARY_DIR + EX_FCST_FILENAME
LONS = read_nc_files(GE, 'lon')
LATS = read_nc_files(GE, 'lat')

## Read in example coarse spatial resolution file for lat and lon over Globe
EX_NMME_FILENAME = '/ex_raw_nmme_download.nc'
GE1 = SUPPLEMENTARY_DIR + EX_NMME_FILENAME
LONI = read_nc_files(GE1, 'X')
LATI = read_nc_files(GE1, 'Y')
LON1 = LONI.copy()
for n, l in enumerate(LON1):
    if l >= 180:
        LON1[n] = LON1[n]-360.
LONI = LON1
LON1 = LONI[0:180]
LON2 = LONI[180:]
LONI = np.hstack((LON2, LON1))

## Read all forecast files
MM = CMN-1
XPREC = np.empty([1, 12, 94, 181, 360])
XPRECI = np.empty([1, 320, 320])

for i,m in enumerate(MODEL):
    if MODEL[i] == 'NCEP-CFSv2':
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL[i], MODEL[i], MON[MM], CYR)
        print('Reading:', INFILE)
        x = read_nc_files(INFILE, 'prec')
        XPREC[0, 0:10, 0:24, :, :] = x[0, 0:10, 0:24, :, :]
    if MODEL[i] == 'NASA-GEOSS2S':
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL[i], MODEL[i], MON[MM], CYR)
        print('Reading:', INFILE)
        XPREC[0, 0:9, 24:34, :, :] = read_nc_files(INFILE, 'prec')
    if MODEL[i] == 'CanCM4i':
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL[i], MODEL[i], MON[MM], CYR)
        print('Reading:', INFILE)
        x = read_nc_files(INFILE, 'prec')
        x1 = np.moveaxis(x, 1, 2)
        XPREC[0, 0:12, 34:44, :, :] = x1
    if MODEL[i] == 'GEM-NEMO':
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL[i], MODEL[i], MON[MM], CYR)
        print('Reading:', INFILE)
        x = read_nc_files(INFILE, 'prec')
        x1 = np.moveaxis(x, 1, 2)
        XPREC[0, 0:12, 44:54, :, :] = x1
    if MODEL[i] == 'COLA-RSMAS-CCSM4':
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL[i], MODEL[i], MON[MM], CYR)
        print('Reading:', INFILE)
        XPREC[0, 0:12, 54:64, :, :] = read_nc_files(INFILE, 'prec')
    if MODEL[i] == 'GFDL-SPEAR':
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL[i], MODEL[i], MON[MM], CYR)
        print('Reading:', INFILE)
        XPREC[0, 0:12, 64:94, :, :] = read_nc_files(INFILE, 'prec')

## Convert units to mm/s or kg/m2/s
XPREC = XPREC/86400

## Reorganize and write
YR = CYR
for m in range(0, 94):
    for l in range(0, 9):
        x = XPREC[0, l, m, :, :]
        print(x.shape, LATI.shape, LONI.shape, LATS.shape, LONS.shape)
        x1 = x[:, 0:180]
        x2 = x[:, 180:]
        x = np.hstack((x2, x1))
        lone, late = np.meshgrid(LONS, LATS)
        #print(lone.shape,late.shape)
        XPRECI[0, :, :] = basemap.interp(x, LONI, LATI, lone, late, \
                        checkbounds=False, masked=False, order=1)
        print(XPRECI.shape)
        print('interpolated')
        jy = YR+LDYR[MM, l]
        l1 = LEADS1[MM, l]
        print('Year:', jy, ',leads:', l1)
        OUTDIR = OUTDIR_TEMPLATE.format(NMME_OUTPUT_DIR, MONTH[MM], YR, m+1)
        OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, MONTH[MM], jy, l1)
        if not os.path.exists(OUTDIR):
            os.makedirs(OUTDIR)

        XPRECI = np.nan_to_num(XPRECI, nan=-9999.)
        LATS = np.nan_to_num(LATS, nan=-9999.)
        LONS = np.nan_to_num(LONS, nan=-9999.)

        SDATE = datetime(YR, MM+1, 1)
        write_3d_netcdf(OUTFILE, XPRECI, 'PRECTOT', 'Downscaled to 0.25deg', \
        'Raw NMME at 1deg', 'kg m-2 s-1', LONS, LATS, SDATE)
        print(f"Writing {OUTFILE}")
