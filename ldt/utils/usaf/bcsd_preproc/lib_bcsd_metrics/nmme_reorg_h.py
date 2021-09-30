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
import numpy as np
from mpl_toolkits import basemap
# pylint: disable=no-name-in-module
import netCDF4 as nc
# pylint: enable=no-name-in-module
from Shrad_modules import read_nc_files


def write_3d_netcdf(infile, var, varname, description, source, \
                    var_units, lons, lats, sdate):
    """write netcdf files"""
    rootgrp = nc.Dataset(infile, 'w', format='NETCDF4')
    #longitude = rootgrp.createDimension('longitude', len(lons))
    #latitude = rootgrp.createDimension('latitude', len(lats))
    time = rootgrp.createDimension('time', None)
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    times = rootgrp.createVariable('time', 'f8', ('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname, 'f4', \
                                     ('time', 'latitude', 'longitude'), \
                                     fill_value=-9999., zlib=True)
    rootgrp.description = description
    rootgrp.history = 'Created ' + time.ctime(time.time())
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
    times[:] = nc.date2num(sdate, units=times.units, calendar=times.calendar)
    rootgrp.close()

MODEL = ['NASA-GEOSS2S', 'CanCM4i', 'GEM-NEMO', 'COLA-RSMAS-CCSM4', 'GFDL-SPEAR']
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
       'Sep', 'Oct', 'Nov', 'Dec']
MONTH = ["jan01", "feb01", "mar01", "apr01", "may01", "jun01",
         "jul01", "aug01", "sep01", "oct01", "nov01", "dec01"]
MONTHN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

print(MODEL[0])

LEADS1 = np.zeros((12, 12), dtype=int)
LDYR = np.zeros((12, 12), dtype=int)

for i in range(0, 12):
    for j in range(0, 12):
        k = MONTHN[i]+j
        KY = 0
        if k >= 13:
            k = k-12
            KY = 1
        LEADS1[i, j] = k
        LDYR[i, j] = KY

DIRA = '/discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/'
DIRB = 'EXPERIMENTS/NMME/data/'
DIR = DIRA + DIRB
INFILE_TEMP = '{}/{}/prec.{}.mon_{}_{:04d}_{:04d}.nc'

DIRC = 'EXPERIMENTS/NMME/data/AF/PRECTOT_Monthly/'
DIR1 = DIRA + DIRC
OUTFILE_TEMPLATE = '{}/{:04d}/{}/ens{}/nmme/{}.nmme.monthly.{:04d}{:02d}.nc'
OUTDIR_TEMPLATE = '{}/{:04d}/{}/ens{}/nmme/'
if not os.path.exists(DIR1):
    os.makedirs(DIR1)

SYR1 = 1982
SYR2 = 1991
SYR3 = 2019
SYR4 = 2011
SYR5 = 2017
SYR6 = 2018
EYR1 = 2020
EYR2 = 2010
EYR3 = 2018
EYR4 = 2017
EYR5 = 2016

GEA = '/discover/nobackup/projects/nca/razamora/GHI_S2S/AFRICOM/data/CFSv2/'
GEB = 'RAW_CFSv2/Monthly/jan01/2008/ens10/jan01.cfsv2.200801.nc'
GE = GEA + GEB
LONS = read_nc_files(GE, 'lon')
LATS = read_nc_files(GE, 'lat')
print("Read fine res lat lon")

GEA1 = '/discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/EXPERIMENTS/'
GEB1 = 'NMME/data/GEM-NEMO/prec.GEM-NEMO.mon_Dec_2019_2020.nc'
GE1 = GEA1 + GEB1
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
print("Read coarse res lat lon")

for mm in range(7, len(MON)):
    xprec = np.empty([39, 12, 49, 181, 360])
    xpreci = np.empty([1, 720, 1440])
    for i in MODEL in enumerate(MODEL):
        if MODEL[i] == 'NASA-GEOSS2S':
            if mm == 0:
                INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i],
                                            MON[mm], SYR1, EYR4)
                print("Reading:", INFILE)
                xprec[0:36, 0:9, 0:4, :, :] = read_nc_files(INFILE, 'prec')
                INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i],
                                            MON[mm], SYR6, EYR1)
                print("Reading:", INFILE)
                xprec[36:39, 0:9, 0:4, :, :] = read_nc_files(INFILE, 'prec')[:, :, 0:4, :, :]
            else:
                INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i],
                                            MON[mm], SYR1, EYR5)
                print("Reading:", INFILE)
                xprec[0:35, 0:9, 0:4, :, :] = read_nc_files(INFILE, 'prec')
                INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i],
                                            MON[mm], SYR5, EYR1)
                print("Reading:", INFILE)
                xprec[35:39, 0:9, 0:4, :, :] = read_nc_files(INFILE, 'prec')[:, :, 0:4, :, :]
        if MODEL[i] == 'CanCM4i':
            INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i], MON[mm], \
                                        SYR1, EYR3)
            print("Reading:", INFILE)
            x = read_nc_files(INFILE, 'prec')
            x1 = np.moveaxis(x, 1, 2)
            xprec[0:37, 0:12, 4:14, :, :] = x1
            INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i], MON[mm], \
                                        SYR3, EYR1)
            print("Reading:", INFILE)
            x = read_nc_files(INFILE, 'prec')
            x1 = np.moveaxis(x, 1, 2)
            xprec[37:39, 0:12, 4:14, :, :] = x1
        if MODEL[i] == 'GEM-NEMO':
            INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i], MON[mm], \
                                        SYR1, EYR3)
            print("Reading:", INFILE)
            x = read_nc_files(INFILE, 'prec')
            x1 = np.moveaxis(x, 1, 2)
            xprec[0:37, 0:12, 14:24, :, :] = x1
            INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i], MON[mm], \
                                        SYR3, EYR1)
            print("Reading:", INFILE)
            x = read_nc_files(INFILE, 'prec')
            x1 = np.moveaxis(x, 1, 2)
            xprec[37:39, 0:12, 14:24, :, :] = x1
        if MODEL[i] == 'COLA-RSMAS-CCSM4':
            INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i], MON[mm], \
                                        SYR1, EYR1)
            print("Reading:", INFILE)
            xprec[0:39, 0:12, 24:34, :, :] = read_nc_files(INFILE, 'prec')
        if MODEL[i] == 'GFDL-SPEAR':
            INFILE = INFILE_TEMP.format(DIR, MODEL[i], MODEL[i], MON[mm], \
                                        SYR2, EYR1)
            print("Reading:", INFILE)
            IH = nc.Dataset(INFILE, mode='r')
            xprec[9:39, 0:12, 34:49, :, :] = read_nc_files(INFILE, 'prec')

    xprec = xprec/86400
    YR = 1981
    for y in range(0, 39):
        YR = YR+1
        for m in range(0, 49):
            for l in range(0, 12):
                x = xprec[y, l, m, :, :]
                #print(x.shape,lati.shape,loni.shape,lats.shape,lons.shape)
                x1 = x[:, 0:180]
                x2 = x[:, 180:]
                x = np.hstack((x2, x1))
                lone, late = np.meshgrid(LONS, LATS)
                #print(lone.shape,late.shape)
                xpreci[0, :, :] = basemap.interp(x, LONI, LATI, lone, late, \
                                                 checkbounds=False, \
                                                 masked=False, order=1)
                print(xpreci.shape)
                print("interpolated")
                jy = YR+LDYR[mm, l]
                l1 = LEADS1[mm, l]
                print("Year:", jy, ",LEADS:", l1)
                OUTFILE = OUTFILE_TEMPLATE.format(DIR1, YR, MONTH[mm], m+13, \
                                                  MONTH[mm], jy, l1)
                OUTDIR = OUTDIR_TEMPLATE.format(DIR1, YR, MONTH[mm], m+13)
                if not os.path.exists(OUTDIR):
                    os.makedirs(OUTDIR)

                xpreci = np.nan_to_num(xpreci, nan=-9999.)
                LATS = np.nan_to_num(LATS, nan=-9999.)
                LONS = np.nan_to_num(LONS, nan=-9999.)

                SDATE = datetime(YR, mm+1, 1)
                write_3d_netcdf(OUTFILE, xpreci, 'PRECTOT', \
                                'Downscaled to 0.25deg', 'Raw NMME at 1deg', \
                                'kg m-2 s-1', LONS, LATS, SDATE)
                print("Writing {}".format(OUTFILE))
