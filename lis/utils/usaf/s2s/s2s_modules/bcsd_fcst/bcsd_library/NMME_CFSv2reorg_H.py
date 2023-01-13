#!/usr/bin/env python
# Author: Shrad Shukla
# coding: utf-8
#Author: Shrad Shukla
#Usage: This is a module for the BCSD code.
#This module bias corrects a forecasts following probability mapping approach as described in Wood et al. 2002
#Date: August 06, 2015
# In[28]:

from __future__ import division
import numpy as np
from All_functions import *
from shrad_modules import read_nc_files
import calendar
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta
from scipy.stats import percentileofscore
from scipy.stats import scoreatpercentile, pearsonr
from scipy import interpolate
from math import *
from netCDF4 import Dataset, date2num, num2date
import time
import os, errno

def write_3d_netcdf(infile, var, varname, DESCRIPTION, SOURCE, VAR_UNITS, SIG_DIGIT, lons, lats, SDATE):
    from datetime import datetime, timedelta
    import netCDF4 as nc
    rootgrp = nc.Dataset(infile, 'w', format='NETCDF4')
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))
    time = rootgrp.createDimension('time', None)
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    times = rootgrp.createVariable('time','f8',('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname,'f4',('time', 'latitude','longitude'), fill_value=-9999., zlib=True)
    import time
    rootgrp.description = DESCRIPTION
    rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.source = SOURCE
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    varname.units = VAR_UNITS
    STRING_DATE = datetime.strftime(SDATE, "%Y-%m-%d")
    times.units = 'days since ' + STRING_DATE
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    varname[:,:,:] = var
    times[:] = nc.date2num(SDATE,units=times.units,calendar=times.calendar)
    rootgrp.close()


month = ["jan01","feb01","mar01","apr01","may01","jun01","jul01","aug01","sep01","oct01","nov01","dec01"]
monthn = [1,2,3,4,5,6,7,8,9,10,11,12]


dir='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/CFSv2/RAW_CFSv2/Monthly/'
infile_temp='{}/{}/{}/ens{}/{}.cfsv2.{:04d}{:02d}.nc'

dir1='/discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/EXPERIMENTS/NMME/data/AF/PRECTOT_Monthly/'
OUTFILE_template = '{}/{:04d}/{}/ens{}/nmme/{}.nmme.monthly.{:04d}{:02d}.nc'
OUTDIR_template = '{}/{:04d}/{}/ens{}/nmme/'
if not os.path.exists(dir1):
 os.makedirs(dir1, exist_ok=True)

#reads sample file for lat and lon values
ge1 = '/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/CFSv2/RAW_CFSv2/Monthly/jan01/2008/ens10/jan01.cfsv2.200801.nc'
lons=read_nc_files(ge1,'lon')
lats=read_nc_files(ge1,'lat')

print("Read coarse res lat lon")



for mm in range(2,len(month)):
 yr1=2007
 for yr in range(0,13):
   yr1=yr1+1
   for ens in range(1,13):
     for ld in range(0,9):
        mm1=mm+1+ld
        xpreci=np.empty([1,720,1440])
        print("Processing: Year=",yr1)
        print("month=",month[mm])
        print("ens=",ens)
        print("lead=",ld)
        INFILE = infile_temp.format(dir, month[mm], yr1, ens, month[mm], yr1, mm1)
        print("Reading:",INFILE)
        xpreci[0,:,:] = read_nc_files(INFILE,'PRECTOT')
        #writing out the prectot files in the NMME folder
        OUTFILE = OUTFILE_template.format(dir1, yr1, month[mm], ens, month[mm], yr1, mm1)
        OUTDIR = OUTDIR_template.format(dir1, yr1, month[mm], ens)
        if not os.path.exists(OUTDIR):
            os.makedirs(OUTDIR, exist_ok=True)

        SDATE = datetime(yr1, mm1, 1)
        write_3d_netcdf(OUTFILE, xpreci, 'PRECTOT', 'Downscaled to 0.25deg', 'Raw NMME at 1deg', 'kg m-2 s-1', 5, lons, lats, SDATE)
        print ("Writing {}".format(OUTFILE))


