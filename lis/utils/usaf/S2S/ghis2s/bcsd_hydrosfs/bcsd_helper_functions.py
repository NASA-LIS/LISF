#!/usr/bin/env python3
"""
Updated forecast_task_01.py, process_forecast_data.py, and convert_forecast_data_to_netcdf.py
"""
#
# Modules
#
import sys
import calendar
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
import xesmf as xe
import netCDF4 as nc4
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
from time import ctime as t_ctime
from time import time as t_time
import cfgrib
import yaml

#
# Custom Modules
#
from dict_variables import get_hydrosfs_name

#
# Helper Functions
#
class VarLimits:
    '''
    This function adjusts minimum and maximum values of a variable to recorded max and minimum.
    Below limits are 6h based
    '''
    def clip_array (self, data_array, var_name=None, min_val=None, max_val=None,
                    missing=None, min_thres=None, precip=None):
        ''' Below limits are 6h based'''
        
        min_limit = {
            get_hydrosfs_name('Pr')     : 1.e-7,
            get_hydrosfs_name('PS')     : 30000.,
            get_hydrosfs_name('T')      : 180.,
            get_hydrosfs_name('LWdown') : 10.,
            get_hydrosfs_name('SWdown') : 0.,
            get_hydrosfs_name('Q')      : 0.,
            get_hydrosfs_name('U')      : -70.,
            get_hydrosfs_name('V')      : -70.
        }

        max_limit = {
            get_hydrosfs_name('Pr')     : 0.04,
            get_hydrosfs_name('PS')     : 110000.,
            get_hydrosfs_name('T')      : 332.,
            get_hydrosfs_name('LWdown') : 700.,
            get_hydrosfs_name('SWdown') : 1367.,
            get_hydrosfs_name('Q')      : 0.05,
            get_hydrosfs_name('U')      : 70.,
            get_hydrosfs_name('V')      : 70.
        }

        if min_thres is not None:
            return min_limit.get(get_hydrosfs_name('Pr'))
        if min_val is None:
            min_val = min_limit.get(var_name)
        if max_val is None:
            max_val = max_limit.get(var_name)
        if missing is None:
            missing = -9999.

        if precip is None:
            clipped_array = np.where(data_array == missing, data_array,
                                     np.clip(data_array, min_val, max_val))
        else:
            # mask identifies values that are less than min_val but not equal to missing
            mask_lt_min = (data_array < min_val) & (data_array != missing)
            data_array[mask_lt_min] = 0.

            # mask identifies values that are greater than max_val but not equal to missing
            mask_gt_max = (data_array > max_val) & (data_array != missing)
            data_array[mask_gt_max] = max_val
            clipped_array = data_array

        return clipped_array

    def __init__ (self):
        self.precip_thres = self.clip_array(np.empty(1), get_hydrosfs_name('Pr'), min_thres = True)

def get_domain_info(configfile, extent=None, coord=None):
    """ gets domain info from LDTINPUT file"""
    # load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    ldtfile = cfg['SETUP']['supplementarydir'] + '/lis_darun/' + cfg['SETUP']['ldtinputfile']
    ldt = nc4_dataset(ldtfile, 'r')

    if extent is not None:
        lon = np.array(ldt['lon'])
        lat = np.array(ldt['lat'])
        return int(np.floor(np.min(lat[:,0]))), int(np.ceil(np.max(lat[:,0]))), \
            int(np.floor(np.min(lon[0,:]))), int(np.ceil(np.max(lon[0,:])))

    if coord is not None:
        lon = np.array(ldt['lon'])
        lat = np.array(ldt['lat'])
        return lat[:,0], lon[0,:]
    
    return None

def read_nc_files(infile, varname):
    ''' reads varname from the nc4 file '''
    nc_ = nc4.Dataset(infile, 'r')
    data = nc_.variables[varname][:]
    nc_.close()
    return data

def write_4d_netcdf(infile, var, varname, description, source, var_units, sig_digit, lons, lats, \
                    ens_num, lead_num, sdate, dates):
    """ writes 4-dimensional netcdf file """
    rootgrp = nc4_dataset(infile, 'w', format='NETCDF4')
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))
    lead = rootgrp.createDimension('Lead', lead_num)
    ens = rootgrp.createDimension('Ens', ens_num)
    time = rootgrp.createDimension('time', None)
    leads = rootgrp.createVariable('Lead','d',('Lead',))
    enss = rootgrp.createVariable('Ens','d',('Ens',))
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    times = rootgrp.createVariable('time','f8',('time',))
    # two dimensions unlimited.
#    varname = rootgrp.createVariable(varname,'f4',('time', 'Lead', 'Ens', 'latitude','longitude'), \
#                                     fill_value=nc4_default_fillvals['f4'], zlib=True, \
#                                     complevel=6, shuffle=True)
    varname = rootgrp.createVariable(varname,'f4',('Lead', 'Ens', 'latitude','longitude'), \
                                     fill_value=-9999., zlib=True, \
                                     complevel=6, shuffle=True)
    rootgrp.description = description
    rootgrp.history = 'Created ' + t_ctime(t_time())
    rootgrp.source = source
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    enss.units = 'unitless'
    varname.units = var_units
    string_date = datetime.strftime(sdate, "%Y-%m-%d")
    times.units = 'days since ' + string_date
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    leads[:]=np.arange(0.5, lead_num+0.5)
    enss[:]=np.arange(0, ens_num)
    varname[:,:,:,:] = var
    times[:] = nc4_date2num(dates,units=times.units,calendar=times.calendar)
    rootgrp.close()

def magnitude(_a, _b):
    """ computes wind magnitude u^2 + v^2 """
    func = lambda x, y: np.sqrt(x**2 + y**2)
    return xr.apply_ufunc(func, _a, _b)

def create_lat_lon_xr_dataset(lon_start, lon_end, lat_start, lat_end, cell_step):
    
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], np.arange(lat_start, lat_end + cell_step, cell_step)),
            "lon": (["lon"], np.arange(lon_start, lon_end + cell_step, cell_step))
        }
    )
    return ds_out

def get_index(ref_array, my_value):
    """
      Function for extracting the index of a Numpy array (ref_array)
      which value is closest to a given number.

      Input parameters:
        - ref_array: reference Numpy array
        - my_value:  floating point number

      Returned value:
        - An integer corresponding to the index
    """
    return np.abs(ref_array - my_value).argmin()

def get_lon_lat_indices(lon_refs, lon_min, lon_max, lat_refs, lat_min, lat_max):
    '''
    '''

    ilon_min = get_index(lon_refs, lon_min)
    ilon_max = get_index(lon_refs, lon_max)
    ilat_min = get_index(lat_refs, lat_min)
    ilat_max = get_index(lat_refs, lat_max)
    
    return ilon_min, ilon_max, ilat_min, ilat_max



# KEEP in bcsd_fcst_helper_functions.py
def create_regridder_weights_file(
    ds_in : xr.Dataset, ds_out : xr.Dataset, method : str, dir_out : str):
    """
    Writes regridder weights to file

    Parameters
    ----------
    ds_in  (xarray Dataset) : Dataset of input grid
    ds_out (xarray Dataset) : Dataset of output grid
    method            (str) : Regridder method
    dir_out           (str) : Directory path of output file

    Examples
    --------
    """
    # Generates regridder
    regridder = xe.Regridder(ds_in, ds_out, method, periodic=True)

    # Get Regridder shape variables
    Ny_in = regridder.shape_in[0]
    Nx_in = regridder.shape_in[1]
    Ny_out = regridder.shape_out[0]
    Nx_out = regridder.shape_out[1]

    # Output Filename
    filename_out = f"{dir_out}/{method}_{Ny_in}x{Nx_in}_{Ny_out}x{Nx_out}.nc"

    # Write to file
    regridder.to_netcdf(filename_out)

def apply_regridder_weights_file(ds_in : xr.Dataset, ds_out : xr.Dataset, variable_dict, dir_supplementary):

    # Resolution
    Ny_in = ds_in['lat'].size
    Nx_in = ds_in['lon'].size
    Ny_out = ds_out['lat'].size
    Nx_out = ds_out['lon'].size
    
    # Filename of regridder weights
    # TODO: This directory should not be hardcoded
    dir_weights = f"{dir_supplementary}/bcsd_fcst"    
    
    output = xr.Dataset()

    for varname in ds_in:
        method = variable_dict[varname]['method']
        filename_weights = f"{dir_weights}/{method}_{Ny_in}x{Nx_in}_{Ny_out}x{Nx_out}.nc"
        
        regridder = xe.Regridder(ds_in, ds_out, method, weights=filename_weights)

        output = output.assign(regridder(ds_in[[varname]]))
    return output