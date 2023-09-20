#!/usr/bin/env python
"""
Created on Fri Nov 15 11:16:48 2013

@author: shrad
Updated by Abheera Hazra feb2018
"""

## This function reads netcdf files.
## It accepts netcdf file name and variable name to be read and returns data
## to an array



from datetime import datetime, timedelta
# pylint: disable=no-name-in-module
import netCDF4 as nc4
# pylint: enable=no-name-in-module
# pylint: disable=no-member, too-many-arguments, unused-variable

def read_nc_files(infile, varname):
    ''' reads varname from the nc4 file '''
    nc_ = nc4.Dataset(infile, 'r')
    data = nc_.variables[varname][:]
    nc_.close()
    return data

def write_netcdf_files(file, var, varname, description, source, var_units,
                       sig_digit, lons, lats, sdate, interval):
    ''' writes netcdf4 file '''
    rootgrp = nc4.Dataset(file, 'w', format='NETCDF4_CLASSIC')
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))
    time = rootgrp.createDimension('time', None)
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    times = rootgrp.createVariable('time','f8',('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname,'f4',('time','latitude','longitude',),
                                     fill_value=nc4.default_fillvals['f4'],
                                     zlib=True,least_significant_digit=sig_digit)
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
    varname[:,:,:] = var
    dates = [sdate+n*timedelta(days=interval) for n in range(varname.shape[0])]
    ## interval is the number of days between two dates
    times[:] = nc4.date2num(dates,units=times.units,calendar=times.calendar)
    rootgrp.close()
