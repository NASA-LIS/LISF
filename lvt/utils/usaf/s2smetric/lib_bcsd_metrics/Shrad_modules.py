# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:16:48 2013

@author: shrad
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
		if e.error == errno.EEXIST:  # file exists error?
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
## This Function runs cdo commands
def run_cdo (COMMAND, INFILE, OUTFILE):
    from subprocess import check_call
    check_call(['cdo', COMMAND, INFILE, OUTFILE])
    
## This function plots spatial plots
def plotmap(var, lats, lons, lat_min, lat_max, lat_int, lon_min, lon_max, lon_int, cmap, cbmin, cbmax, plottype, land_color, 
                ocean_color, levels, labels_1, labels_2, col_under, col_higher, extend, area_thresh=10000):
    #importing required modules               
               
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib.pyplot as mpl
    
    if plottype == 'mesh':
        lonres = (lons.max()-lons.min())/(len(lons)-1)
        latres = (lats.max()-lats.min())/(len(lats)-1)
        lats = lats - 0.5*latres
        lons = lons - 0.5*lonres
    elif plottype == 'contour':
        pass
    elif plottype == 'color':
        pass
    else:
        raise ValueError('{}: Not a valid option for plottype'.
                         format(plottype))

    m = Basemap(projection='cyl', llcrnrlat=lat_min, llcrnrlon=lon_min,
                urcrnrlat=lat_max, urcrnrlon=lon_max, rsphere=6371200., resolution='i',
                area_thresh=area_thresh)
    xi, yi = m(lons, lats)
    xi, yi = np.meshgrid(xi, yi)

    if plottype == 'mesh':
        cs=m.pcolormesh(xi, yi, var, vmax=cbmax, vmin=cbmin, cmap=cmap)
    elif plottype == 'contour':
        cmap.set_under(col_under)
        cmap.set_over(col_higher)
        norm = mpl.pyplot.cm.colors.Normalize(vmax=cbmax, vmin=cbmin, clip=False)
        cs=m.contourf(xi, yi, var, levels, cmap=cmap,xnorm=norm, extend=extend)
    elif plottype == 'color':
        cs=m.pcolor(xi, yi, var, vmax=cbmax, vmin=cbmin, cmap=cmap)

    #m.drawlsmask(land_color=land_color, ocean_color=ocean_color, lakes=True)
    m.drawlsmask(ocean_color=ocean_color, lakes=True)
    m.drawparallels(np.arange(-80,81,lat_int),labels=labels_1,fontsize=8, linewidth=0.3, rotation=45)
    m.drawmeridians(np.arange(0,360,lon_int),labels=labels_2,fontsize=8, linewidth=0.3)    
    m.drawcoastlines(linewidth=0.75)
    m.drawcountries(linewidth=0.75)
    mpl.rc('xtick', labelsize=11)
    mpl.rc('ytick', labelsize=11)
    #if (colorbar):
        #cbar=m.colorbar(cs, location=colorbar_location, shrink=shrink, extend=extend, pad=pad)
    #if (label_flag):
     #   cbar.set_label(cbar_label, fontsize=14)
    return m, cs  
    
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
