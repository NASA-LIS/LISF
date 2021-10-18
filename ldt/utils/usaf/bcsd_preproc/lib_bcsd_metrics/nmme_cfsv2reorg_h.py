#!/usr/bin/env python
"""
# Author: Abheera Hazra
#This module reorganizes
#NMME preciptation forecasts
#Date: May 06, 2021
# In[28]:
"""
from __future__ import division
#import sys
from datetime import datetime
import os
import numpy as np
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


MONTH = ["jan01", "feb01", "mar01", "apr01", "may01", "jun01", \
         "jul01", "aug01", "sep01", "oct01", "nov01", "dec01"]
MONTHN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]


DIRA = '/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/'
DIRB = 'AFRICOM/data/CFSv2/RAW_CFSv2/Monthly/'
DIR = DIRA + DIRB
INFILE_TEMP = '{}/{}/{}/ens{}/{}.cfsv2.{:04d}{:02d}.nc'

DIRA1 = '/discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/'
DIRB1 = 'EXPERIMENTS/NMME/data/AF/PRECTOT_Monthly/'
DIR1 = DIRA1 + DIRB1
OUTFILE_TEMPLATE = '{}/{:04d}/{}/ens{}/nmme/{}.nmme.monthly.{:04d}{:02d}.nc'
OUTDIR_TEMPLATE = '{}/{:04d}/{}/ens{}/nmme/'
if not os.path.exists(DIR1):
    os.makedirs(DIR1)

#reads sample file for lat and lon values
GEA1 = '/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/'
GEB1 = 'data/CFSv2/RAW_CFSv2/Monthly/jan01/2008/ens10/jan01.cfsv2.200801.nc'
GE1 = GEA1 + GEB1
LONS = read_nc_files(GE1, 'lon')
LATS = read_nc_files(GE1, 'lat')

print("Read coarse res lat lon")

for mm in range(2, len(MONTH)):
    YR1 = 2007
    for YR in range(0, 13):
        YR1 = YR1+1
        for ens in range(1, 13):
            for ld in range(0, 9):
                mm1 = mm+1+ld
                xpreci = np.empty([1, 720, 1440])
                print("Processing: Year=", YR1)
                print("MONTH=", MONTH[mm])
                print("ens=", ens)
                print("lead=", ld)
                INFILE = INFILE_TEMP.format(DIR, MONTH[mm], YR1, ens, \
                                            MONTH[mm], YR1, mm1)
                print("Reading:", INFILE)
                xpreci[0, :, :] = read_nc_files(INFILE, 'PRECTOT')
                #writing out the prectot files in the NMME folder
                OUTFILE = OUTFILE_TEMPLATE.format(DIR1, YR1, MONTH[mm], \
                                                  ens, MONTH[mm], YR1, mm1)
                OUTDIR = OUTDIR_TEMPLATE.format(DIR1, YR1, MONTH[mm], ens)
                if not os.path.exists(OUTDIR):
                    os.makedirs(OUTDIR)

                SDATE = datetime(YR1, mm1, 1)
                write_3d_netcdf(OUTFILE, xpreci, 'PRECTOT', \
                                'Downscaled to 0.25deg', \
                                'Raw NMME at 1deg', 'kg m-2 s-1', \
                                LONS, LATS, SDATE)
                print("Writing {}".format(OUTFILE))
