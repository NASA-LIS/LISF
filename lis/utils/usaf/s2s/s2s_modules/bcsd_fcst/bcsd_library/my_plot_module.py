# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:16:48 2013

@author: shrad
"""

## This function reads netcdf files.
## It accepts netcdf file name and variable name to be read and returns data 
## to an array
from __future__ import division 
    
def write_2_netcdf(outfile, var, varname, DESCRIPTION, SOURCE, VAR_UNITS, SIG_DIGIT, lons, lats, SDATE, dates):
    from datetime import datetime, timedelta
    import netCDF4 as nc
    rootgrp = nc.Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))
    time = rootgrp.createDimension('time', None)
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    times = rootgrp.createVariable('time','f8',('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname,'f4',('time','latitude','longitude',), fill_value=nc.default_fillvals['f4'], zlib=True,least_significant_digit=SIG_DIGIT)
    #varname = rootgrp.createVariable(varname,'f4',('time','latitude','longitude',), fill_value=nc.default_fillvals['f4'])
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
    times[:] = nc.date2num(dates,units=times.units,calendar=times.calendar)
    rootgrp.close()

def write_4d_netcdf(infile, var, varname, DESCRIPTION, SOURCE, VAR_UNITS, SIG_DIGIT, lons, lats, ENS_NUM, SDATE, dates):
    from datetime import datetime, timedelta
    import netCDF4 as nc
    import numpy as np
    rootgrp = nc.Dataset(infile, 'w', format='NETCDF4')
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))
    ens = rootgrp.createDimension('Ens', ENS_NUM)
    time = rootgrp.createDimension('time', None)
    enss = rootgrp.createVariable('Ens','d',('Ens',))
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    times = rootgrp.createVariable('time','f8',('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname,'f4',('time', 'Ens', 'latitude','longitude'), fill_value=nc.default_fillvals['f4'], zlib=True)
    import time
    rootgrp.description = DESCRIPTION
    rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.source = SOURCE
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    enss.units = '-'
    varname.units = VAR_UNITS
    STRING_DATE = datetime.strftime(SDATE, "%Y-%m-%d")
    times.units = 'days since ' + STRING_DATE
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    enss[:]=np.arange(0, ENS_NUM)
    varname[:,:,:,:] = var
    times[:] = nc.date2num(dates,units=times.units,calendar=times.calendar)
    rootgrp.close()

def write_5d_netcdf(infile, var, varname, DESCRIPTION, SOURCE, VAR_UNITS, SIG_DIGIT, lons, lats, ENS_NUM, LEAD_NUM, SDATE, dates):
    from datetime import datetime, timedelta
    import netCDF4 as nc
    import numpy as np
    rootgrp = nc.Dataset(infile, 'w', format='NETCDF4')
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))
    lead = rootgrp.createDimension('Lead', LEAD_NUM)
    ens = rootgrp.createDimension('Ens', ENS_NUM)
    time = rootgrp.createDimension('time', None)
    leads = rootgrp.createVariable('Lead','d',('Lead',))
    enss = rootgrp.createVariable('Ens','d',('Ens',))
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    times = rootgrp.createVariable('time','f8',('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname,'f4',('time', 'Lead', 'Ens', 'latitude','longitude'), fill_value=nc.default_fillvals['f4'], zlib=True)
    import time
    rootgrp.description = DESCRIPTION
    rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.source = SOURCE
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    enss.units = '-'
    varname.units = VAR_UNITS
    STRING_DATE = datetime.strftime(SDATE, "%Y-%m-%d")
    times.units = 'days since ' + STRING_DATE
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    leads[:]=np.arange(0, LEAD_NUM)
    enss[:]=np.arange(0, ENS_NUM)
    varname[:,:,:,:,:] = var
    times[:] = nc.date2num(dates,units=times.units,calendar=times.calendar)
    rootgrp.close()
def write_5d_netcdf(infile, var, varname, DESCRIPTION, SOURCE, VAR_UNITS, SIG_DIGIT, lons, lats, ENS_NUM, LEAD_NUM, SDATE, dates):
    from datetime import datetime, timedelta
    import netCDF4 as nc
    import numpy as np
    rootgrp = nc.Dataset(infile, 'w', format='NETCDF4')
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))
    lead = rootgrp.createDimension('Lead', LEAD_NUM)
    ens = rootgrp.createDimension('Ens', ENS_NUM)
    time = rootgrp.createDimension('time', None)
    leads = rootgrp.createVariable('Lead','d',('Lead',))
    enss = rootgrp.createVariable('Ens','d',('Ens',))
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    times = rootgrp.createVariable('time','f8',('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname,'f4',('time', 'Lead', 'Ens', 'latitude','longitude'), fill_value=nc.default_fillvals['f4'], zlib=True)
    import time
    rootgrp.description = DESCRIPTION
    rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.source = SOURCE
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    enss.units = '-'
    varname.units = VAR_UNITS
    STRING_DATE = datetime.strftime(SDATE, "%Y-%m-%d")
    times.units = 'days since ' + STRING_DATE
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    leads[:]=np.arange(0, LEAD_NUM)
    enss[:]=np.arange(0, ENS_NUM)
    varname[:,:,:,:,:] = var
    times[:] = nc.date2num(dates,units=times.units,calendar=times.calendar)
    rootgrp.close()
