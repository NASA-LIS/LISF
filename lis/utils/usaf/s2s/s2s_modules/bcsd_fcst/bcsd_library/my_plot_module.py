# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:16:48 2013

@author: shrad
"""

## This function reads netcdf files.
## It accepts netcdf file name and variable name to be read and returns data 
## to an array
from __future__ import division 

## This function plots spatial plots
def my_plotmap(var, lats, lons, lat_min, lat_max, lat_int, lon_min, lon_max, lon_int, cmap, cbmin, cbmax, plottype, levels, labels_1, labels_2, col_under, col_higher, extend, area_thresh=10000):
    #importing required modules               
               
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    

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
	    cmap = cmap
	    cmap.set_under(col_under)
	    cmap.set_over(col_higher)
	    #norm = mpl.pyplot.cm.colors.Normalize(vmax=cbmax, vmin=cbmin, clip=False)
	    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
	    cs=m.contourf(xi, yi, var, levels, cmap=cmap, xnorm=norm, extend=extend)
    elif plottype == 'color':
	    cs=m.pcolor(xi, yi, var, vmax=cbmax, vmin=cbmin, cmap=cmap)

    #m.drawlsmask(ocean_color='gray', lakes=True)
    m.drawparallels(np.arange(-80,81,lat_int),labels=labels_1,fontsize=8, linewidth=0.3)
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
