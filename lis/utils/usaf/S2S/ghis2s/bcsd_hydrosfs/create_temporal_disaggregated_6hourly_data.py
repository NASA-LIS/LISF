#!/usr/bin/env python
"""

"""

#
# Modules
#
import os
import sys
import yaml
from datetime import datetime
import calendar
from time import ctime as t_ctime
from time import time as t_time
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num

#
# Custom Modules
#
from bcsd_helper_functions import get_domain_info, get_lon_lat_indices
from bcsd_helper_functions import VarLimits as lim
from bcsd_filename_functions import (generate_raw_monthly_forecast_filename, generate_bcsd_monthly_forecast_filename, generate_raw_6hourly_forecast_filename, generate_bcsd_6hourly_forecast_filename)
from dict_variables import get_hydrosfs_name

def _read_cmd_args():
    """Read command line arguments."""

    with open(sys.argv[1], 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 10:
        print('[ERR] Invalid number of command line arguments!')
        _usage()
        sys.exit(1)

    args = {
        'config_filename' : str(sys.argv[1]),
        'fcst_init_year'  : int(sys.argv[2]),
        'fcst_init_month' : int(sys.argv[3]),
        'ensemble_number' : int(sys.argv[4]),
        'lead_month'      : int(sys.argv[5]),     
        'lead_day'        : int(sys.argv[6]),     
        'obs_var'         : str(sys.argv[7]),
        'bc_var'          : str(sys.argv[8]),
        'unit'            : str(sys.argv[9]),
        'dir_main'        : str(config['SETUP']['DIR_MAIN']),
        'fcst_name'       : str(config['BCSD']['fcst_data_type']),        
        'dataset_type'    : str(config['SETUP']['DATATYPE']),
        'clim_syr'        : int(config['BCSD']['clim_start_year']),
        'clim_eyr'        : int(config['BCSD']['clim_end_year'])
    }
    
    args['config'] = config

    return args


def scale_forcings (mon_bc_value, mon_raw_value, input_raw_data, bc_var = None):
    ''' perform scaling '''
    output_bc_data = np.ones(len(input_raw_data))*-9999.

    if bc_var == 'PRCP':
        if mon_raw_value == 0.:
            output_bc_data[:] = mon_bc_value
        else:
            correction_factor = mon_bc_value/mon_raw_value
            output_bc_data[:] = input_raw_data[:]*correction_factor
    else:
        correction_factor = mon_bc_value - mon_raw_value
        output_bc_data[:] = input_raw_data[:] + correction_factor

    return output_bc_data

def write_bc_netcdf(outfile, var, varname, description, source, var_units, \
var_standard_name, lons, lats, sdate, dates, sig_digit, north_east_corner_lat, \
north_east_corner_lon, south_west_corner_lat, south_west_corner_lon, \
resolution_x, resolution_y, time_increment):
    """write netcdf"""
    rootgrp = nc4_dataset(outfile, 'w', format='NETCDF4_CLASSIC')
    time = rootgrp.createDimension('time', None)
    longitude = rootgrp.createDimension('lon', len(lons))
    latitude = rootgrp.createDimension('lat', len(lats))

    longitudes = rootgrp.createVariable('lon', 'f4', ('lon',))
    latitudes = rootgrp.createVariable('lat', 'f4', ('lat',))
    times = rootgrp.createVariable('time', 'f4', ('time', ))

    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname, 'f4', ('time', 'lat', \
    'lon',), fill_value=-9999, zlib=True, \
    least_significant_digit=sig_digit)
    rootgrp.missing_value = -9999
    rootgrp.description = description
    rootgrp.zenith_interp = "true,false,"
    rootgrp.MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL"
    rootgrp.conventions = "CF-1.6"
    rootgrp.south_west_corner_lat = float(south_west_corner_lat)
    rootgrp.south_west_corner_lon = float(south_west_corner_lon)
    rootgrp.north_east_corner_lat = float(north_east_corner_lat)
    rootgrp.north_east_corner_lon = float(north_east_corner_lon)
    rootgrp.DX = resolution_x
    rootgrp.DY = resolution_y
    rootgrp.history = 'Created ' + t_ctime(t_time())
    rootgrp.source = source
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    varname.units = var_units
    varname.standard_name = var_standard_name
    string_date = datetime.strftime(sdate, "%Y-%m-%d %H:%M:%S")
    times.units = 'minutes since ' + string_date
    times.time_increment = time_increment
    times.begin_date = datetime.strftime(sdate, "%Y%m%d")
    times.begin_time = '000000'
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    varname[:, :, :] = var
    times[:] = nc4_date2num(dates, units=times.units, calendar=times.calendar)
    rootgrp.close()

    
#
# Main
#

# Read command line arguments
args = _read_cmd_args()

# Input Args
CONFIG_FILE = args['config_filename']
fcst_init_year = args['fcst_init_year']
fcst_init_month = args['fcst_init_month']
ensemble_number = args['ensemble_number']
lead_month = args['lead_month']
lead_day = args['lead_day']
OBS_VAR = args['obs_var']
BC_VAR = args['bc_var']
UNIT = args['unit']

# Config Args
config = args['config']
dir_main = args['dir_main']
fcst_name = args['fcst_name']
dataset_type = args['dataset_type']
CLIM_SYR = args['clim_syr']
CLIM_EYR = args['clim_eyr']

# Coordinate Information
LAT1, LAT2, LON1, LON2 = get_domain_info(CONFIG_FILE, extent = True)
LATS, LONS = get_domain_info(CONFIG_FILE, coord = True)

# Get the lat/lon indexes for the ranges
ilon_min, ilon_max, ilat_min, ilat_max = get_lon_lat_indices(LONS, LON1, LON2, LATS, LAT1, LAT2)

FCST_DATE = datetime(fcst_init_year, fcst_init_month, 1, 6) + relativedelta(days = lead_day)

# Number of subdaily time steps in the target forecast month
NUM_TIMESTEPS = 4

# Sub-Daily raw data
# 1 lead month, 1 ensemble member, and 4 timesteps (1 day) per file
SUBDAILY_INFILE = generate_raw_6hourly_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number, lead_day)
print(f"Reading raw sub-daily forecast: {SUBDAILY_INFILE} \n")
INPUT_RAW_DATAG = xr.open_dataset(SUBDAILY_INFILE)
INPUT_RAW_DATA = INPUT_RAW_DATAG.sel(lon=slice(LON1,LON2),lat=slice(LAT1,LAT2))
INPUT_RAW_DATA = INPUT_RAW_DATA[OBS_VAR][:,:,:]


# Monthly raw data
# 1 lead month and 1 ensemble member per file
MONTHLY_INFILE = generate_raw_monthly_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number, lead_month)

print(f"Reading raw monthly forecast: {MONTHLY_INFILE} \n")
MONTHLY_INPUT_RAW_DATAG = xr.open_dataset(MONTHLY_INFILE)
MONTHLY_INPUT_RAW_DATA = MONTHLY_INPUT_RAW_DATAG.sel(lon = slice(LON1,LON2), lat = slice(LAT1,LAT2))
MONTHLY_INPUT_RAW_DATA = MONTHLY_INPUT_RAW_DATA[OBS_VAR][:,:]


# Monthly bcsd data
# has all 4 lead months and all ensemble members in single file
BC_INFILE = generate_bcsd_monthly_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, OBS_VAR)

print(f"Reading bcsd monthly forecast: {BC_INFILE} \n")
MON_BC_DATAG = xr.open_dataset(BC_INFILE)
MON_BC_DATAG = MON_BC_DATAG.rename({"longitude": "lon", "latitude" : "lat"})
MON_BC_DATA = MON_BC_DATAG.sel(lon=slice(LON1,LON2),lat=slice(LAT1,LAT2))
MON_BC_VALUE = MON_BC_DATA[OBS_VAR][lead_month, ensemble_number - 1,:,:]


# make sure lat/lon are aligned.
if (not np.array_equal(MONTHLY_INPUT_RAW_DATA["lat"].values,
                       MON_BC_VALUE["lat"].values)) or \
                       (not np.array_equal(MONTHLY_INPUT_RAW_DATA["lon"].values,
                                           MON_BC_VALUE["lon"].values)):
    MONTHLY_INPUT_RAW_DATA["lon"] = MON_BC_VALUE["lon"].values
    MONTHLY_INPUT_RAW_DATA["lat"] = MON_BC_VALUE["lat"].values
if (not np.array_equal(INPUT_RAW_DATA["lat"].values,
                       MON_BC_VALUE["lat"].values)) or \
                       (not np.array_equal(INPUT_RAW_DATA["lon"].values,
                                           MON_BC_VALUE["lon"].values)):
    INPUT_RAW_DATA["lon"] = MON_BC_VALUE["lon"].values
    INPUT_RAW_DATA["lat"] = MON_BC_VALUE["lat"].values

    
correct = xr.apply_ufunc(
    scale_forcings,
    MON_BC_VALUE.chunk({"lat": "auto", "lon": "auto"}).compute(),
    MONTHLY_INPUT_RAW_DATA.chunk({"lat": "auto", "lon": "auto"}).compute(),
    INPUT_RAW_DATA.chunk({"lat": "auto", "lon": "auto"}).compute(),
    input_core_dims=[[],[],['time']],
    exclude_dims=set(('time',)),
    output_core_dims=[['time']],
    vectorize=True,
    dask="forbidden",
    output_dtypes=[np.float64],
    kwargs={'bc_var': BC_VAR},
)

correct2 = np.moveaxis(correct.values,2,0)

# Output Data
OUTPUT_BC_DATA = np.ones((NUM_TIMESTEPS, len(LATS), len(LONS)))*-9999.

OUTPUT_BC_DATA[:,ilat_min:ilat_max + 1, ilon_min:ilon_max + 1] = correct2[:,:,:]



# clip limits - monthly BC NMME precip:
limits = lim()

if OBS_VAR == get_hydrosfs_name('Pr'):
    OUTPUT_BC_DATA = limits.clip_array(OUTPUT_BC_DATA, var_name=OBS_VAR, precip=True)
else:
    OUTPUT_BC_DATA = limits.clip_array(OUTPUT_BC_DATA, var_name=OBS_VAR)



# apply mask
OUTPUT_BC_DATA = np.ma.masked_array(OUTPUT_BC_DATA, mask = OUTPUT_BC_DATA == -9999.)




date = [FCST_DATE + relativedelta(hours = n * 6) for n in range(NUM_TIMESTEPS)]



# Create output file directory
file_directory_out = generate_bcsd_6hourly_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number, lead_day, 
    OBS_VAR, 'file_directory')
if not os.path.exists(file_directory_out):
    os.makedirs(file_directory_out)

# Write output file
file_path_out = generate_bcsd_6hourly_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number, lead_day, OBS_VAR)

print(f"Now writing {file_path_out}")
write_bc_netcdf(
    file_path_out, OUTPUT_BC_DATA, OBS_VAR, 'Bias corrected forecasts', 'MODEL:' + fcst_name, \
                UNIT, OBS_VAR, LONS, LATS, FCST_DATE, date, 8, 89.875, 179.875, -89.875, \
                -179.875, 0.05, 0.05, 21600)
