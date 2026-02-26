#!/usr/bin/env python
"""

"""

#
# Modules
#
import sys
from datetime import datetime
import os
import yaml
import xarray as xr
import numpy as np

#
# Custom Modules
#
from bcsd_stats_functions import calc_bcsd
from bcsd_helper_functions import (
    read_nc_files,
    write_4d_netcdf,
    get_domain_info,
    get_lon_lat_indices
)
from bcsd_helper_functions import VarLimits as lim
from bcsd_filename_functions import (
    generate_raw_monthly_forecast_filename,
    generate_raw_climatology_forecast_filename,
    generate_raw_climatology_hydroscs_filename,
    generate_bcsd_monthly_forecast_filename
)

#
# Functions
#
def _usage():
    """Print command line usage."""
    txt =  f'[INFO] Usage: {sys.argv[0]}'
    txt += 'config_filename fcst_init_year fcst_init_month obs_var bc_var unit'
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""

    with open(sys.argv[1], 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 7:
        print('[ERR] Invalid number of command line arguments!')
        _usage()
        sys.exit(1)

    args = {
        'config_filename' : str(sys.argv[1]),
        'fcst_init_year'  : int(sys.argv[2]),
        'fcst_init_month' : int(sys.argv[3]),
        'obs_var'         : str(sys.argv[4]),
        'bc_var'          : str(sys.argv[5]),
        'unit'            : str(sys.argv[6]),
        'dir_main'        : str(config['SETUP']['DIR_MAIN']),
        'fcst_name'       : str(config['BCSD']['fcst_data_type']),        
        'dataset_type'    : str(config['SETUP']['DATATYPE']),
        'lead_months'     : int(config['EXP']['lead_months']),
        'ens_num'         : int(config['BCSD']['nof_raw_ens']),
        'clim_syr'        : int(config['BCSD']['clim_start_year']),
        'clim_eyr'        : int(config['BCSD']['clim_end_year'])
    }

    args['config'] = config

    return args

def latlon_calculations(ilat_min, ilat_max, ilon_min, ilon_max, \
                        np_obs_clim_array, np_fcst_clim_array, \
                        lead_final, \
                        ens_num, fcst_init_month, bc_var, \
                        tiny, fcst_coarse, correct_fcst_coarse):

    """Lat and Lon"""
    num_lats = ilat_max - ilat_min + 1
    num_lons = ilon_max - ilon_min + 1

    print("num_lats = ", num_lats, np_obs_clim_array.shape)
    print("num_lons = ", num_lons, fcst_coarse.shape)

    for ilat in range(num_lats):
        lat_num = ilat_min + ilat
        for ilon in range(num_lons):
            lon_num = ilon_min + ilon

            ## First read Observed clim data (all months available in one file)
            ## so don't have to read it again for each lead time
            obs_clim_all = np_obs_clim_array[:, :, ilat, ilon]

            ## Now read forecast climatology data too.
            fcst_clim_all = np_fcst_clim_array[:, :, ilat, ilon]

            target_fcst_val_arr = fcst_coarse[:, :, lat_num, lon_num]

            #print("shape of FCST_COARSE: ", TARGET_FCST_VAL_ARR.shape)

            correct_fcst_coarse[:, :, lat_num, lon_num] = calc_bcsd(obs_clim_all, fcst_clim_all,
                                                                    lead_final, target_fcst_val_arr,
                                                                    ens_num, fcst_init_month,
                                                                    bc_var, tiny)

#
# Main
#

args = _read_cmd_args()

# Input Args
CONFIG_FILE = args['config_filename']
fcst_init_year = args['fcst_init_year']
fcst_init_month = args['fcst_init_month']
OBS_VAR = args['obs_var']
BC_VAR = args['bc_var']
UNIT = args['unit']

# Config Args
dir_main = args['dir_main']
fcst_name = args['fcst_name']
dataset_type = args['dataset_type']
LEAD_FINAL = args['lead_months']
ENS_NUM = args['ens_num']
CLIM_SYR = args['clim_syr']
CLIM_EYR = args['clim_eyr']

# print("fcst_init_year: ", fcst_init_year)
# print("fcst_init_month: ", fcst_init_month)
# print("OBS_VAR: ", OBS_VAR)
# print("BC_VAR: ", BC_VAR)
# print("UNIT: ", UNIT)

# print("dir_main: ", dir_main)
# print("fcst_name: ", fcst_name)
# print("dataset_type: ", dataset_type)
# print("LEAD_FINAL: ", LEAD_FINAL)
# print("ENS_NUM: ", ENS_NUM)
# print("CLIM_SYR: ", CLIM_SYR)
# print("CLIM_EYR: ", CLIM_EYR)

# Coordinate Information
LAT1, LAT2, LON1, LON2 = get_domain_info(CONFIG_FILE, extent=True)
LATS, LONS = get_domain_info(CONFIG_FILE, coord=True)

TINY = ((1 / (CLIM_EYR - CLIM_SYR + 1)) / ENS_NUM) / 2   # Adjust quantile, if it is out of bounds

# First read observed climatology for the given variable
file_path_obs_climatology = generate_raw_climatology_hydroscs_filename(
    dir_main, 'hindcast', OBS_VAR)

OBS_CLIM_ARRAY = xr.open_dataset(file_path_obs_climatology)


# First read forecast climatology for the given variable and forecast initialzation month
fcst_clim_infile = generate_raw_climatology_forecast_filename(
    dir_main, fcst_name, 'hindcast', fcst_init_month, OBS_VAR)
fcst_clim_array = xr.open_dataset(fcst_clim_infile)
print(f"Reading forecast climatology {fcst_clim_infile}")


#First read raw forecasts
fcst_coarse = np.empty((LEAD_FINAL, ENS_NUM, len(LATS), len(LONS)))
for lead_num in range(0, LEAD_FINAL): ## Loop from lead =0 to Final Lead
    for ens in range(1, ENS_NUM + 1):
        infile = generate_raw_monthly_forecast_filename(
            dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ens,
            lead_num)
        print(infile)

        fcst_coarse[lead_num, ens - 1, ] = read_nc_files(infile, OBS_VAR)[:]

# Defining array to store bias-corrected monthly forecasts
correct_fcst_coarse = np.ones((LEAD_FINAL, ENS_NUM, len(LATS), len(LONS)))*-9999.
print("shape of fcst_coarse: ", fcst_coarse.shape)

# Get the lat/lon indexes for the ranges
ilon_min, ilon_max, ilat_min, ilat_max = get_lon_lat_indices(LONS, LON1, LON2, LATS, LAT1, LAT2)
nlats = len(LATS)
nlons = len(LONS)

# Get the values (Numpy array) for the lat/lon ranges
np_obs_clim_array = OBS_CLIM_ARRAY.clim.sel(longitude=slice(LON1, LON2), \
                    latitude=slice(LAT1, LAT2)).values
np_fcst_clim_array = fcst_clim_array.clim.sel(longitude=slice(LON1, LON2), \
                     latitude=slice(LAT1, LAT2)).values

print("Latitude:  ", nlats, ilat_min, ilat_max)
print("Longitude: ", nlons, ilon_min, ilon_max)
print("np_obs_clim_array:", np_obs_clim_array.shape, type(np_obs_clim_array))
print("np_fcst_clim_array:", np_fcst_clim_array.shape, type(np_fcst_clim_array))

latlon_calculations(ilat_min, ilat_max, ilon_min, ilon_max, np_obs_clim_array, np_fcst_clim_array, 
                    LEAD_FINAL, ENS_NUM, fcst_init_month, BC_VAR, TINY, fcst_coarse, 
                    correct_fcst_coarse)

print("1. correct_fcst_coarse:", correct_fcst_coarse.shape)
correct_fcst_coarse = np.ma.masked_array(correct_fcst_coarse, mask=correct_fcst_coarse == -9999.)
print("2. correct_fcst_coarse:", correct_fcst_coarse.shape)

# clip limits - monthly BC NMME precip:
limits = lim()
correct_fcst_coarse = limits.clip_array(correct_fcst_coarse, var_name=OBS_VAR)
print("3. correct_fcst_coarse:", correct_fcst_coarse.shape)


# Create output file directory
file_directory_monthly_out =  generate_bcsd_monthly_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, OBS_VAR, 'file_directory')
if not os.path.exists(file_directory_monthly_out):
    os.makedirs(file_directory_monthly_out)

# Write output file
file_path_monthly_out = generate_bcsd_monthly_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, OBS_VAR)

sdate = datetime(fcst_init_year, fcst_init_month, 1)
dates = [sdate]
write_4d_netcdf(file_path_monthly_out, correct_fcst_coarse, OBS_VAR, fcst_name, \
'Bias corrected', UNIT, 5, LONS, LATS, ENS_NUM, LEAD_FINAL, sdate, dates)
