# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:16:48 2013

@author: shrad
Updated by Abheera Hazra feb2018
"""

## This function reads netcdf files.
## It accepts netcdf file name and variable name to be read and returns data 
## to an array
from __future__ import division 

# The function below calculates rsquare
def CALC_RSQ(DATA1, DATA2):
	from scipy.stats import pearsonr
	import numpy as np
	CORR = np.ones((DATA1.shape[1], DATA1.shape[2]))*-99
	for i in range(DATA1.shape[1]):
		for j in range(DATA1.shape[2]):
			if np.isfinite(DATA1[0, i, j]) and np.mean(DATA1[:, i, j] != 0):
				OBS_ANOM = DATA1[:, i, j]-np.mean(DATA1[:, i, j])
				FCST_ANOM = DATA2[:, i, j]-np.mean(DATA2[:, i, j])
				## Calculating correlaiton between anomaly
				CORR[i, j] = pearsonr(OBS_ANOM, FCST_ANOM)[0, ]
	return CORR

def MAKEDIR(path):
	import errno, os, sys
	try:
		os.makedirs(path, 493)
	except OSError as e:
		if e.errno == errno.EEXIST:  # file exists error?
			print('warning:',path,' exists')
		else:
			raise  # re-raise the exception
		# make sure its mode is right
		os.chmod(path, 493)

def read_nc_files(INFILE, VARNAME):
    from netCDF4 import Dataset
    nc = Dataset(INFILE, 'r')
    data = nc.variables[VARNAME][:]
    nc.close()
    return(data)
        
def write_netcdf_files(file, var, varname, DESCRIPTION, SOURCE, VAR_UNITS, SIG_DIGIT, lons, lats, SDATE, interval):
    from datetime import datetime, timedelta
    import netCDF4 as nc
    rootgrp = nc.Dataset(file, 'w', format='NETCDF4_CLASSIC')
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
    dates = [SDATE+n*timedelta(days=interval) for n in range(varname.shape[0])] ## interval is the number of days between two dates
    times[:] = nc.date2num(dates,units=times.units,calendar=times.calendar)
    rootgrp.close()
