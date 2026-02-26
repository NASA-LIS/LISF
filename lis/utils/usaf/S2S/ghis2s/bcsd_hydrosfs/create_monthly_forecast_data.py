#!/usr/bin/python3
"""
Creates monthly forecast data

Workflow
--------
1. Get grib2 filename and determine if file needs to be patched/replaced
2. Apply Regridder to file
3. [OPTIONAL] Apply lapse-rate correction
4. [OPTIONAL] Apply slope-aspect correction
5. [OPTIONAL] Apply landmask
6. [OPTIONAL] Clip unfeasable values
7. Select only time steps for the number of months wanted
8. Generate monthly mean of current lead_month data
9. Write variables to NetCDF
"""

#
# Modules
#
import sys
import os
import yaml
import xarray as xr
import pandas as pd

#
# Custom Modules
#
from bcsd_filename_functions import (
    generate_raw_monthly_forecast_filename,
    generate_raw_6hourly_forecast_filename
)

from dict_variables import get_all_hydrosfs_names

#
# Functions
#
def _usage():
    """Print command line usage."""
    txt =  f'[INFO] Usage: {sys.argv[0]}'
    txt += 'config_filename fcst_init_year fcst_init_month ensemble_number'
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""

    with open(sys.argv[1], 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 6:
        print('[ERR] Invalid number of command line arguments!')
        _usage()
        sys.exit(1)

    args = {
        'config_filename'   : sys.argv[1],
        'fcst_init_year'    : sys.argv[2], #2015
        'fcst_init_month'   : sys.argv[3], #5
        'ensemble_number'   : sys.argv[4], #1
        'lead_month'        : sys.argv[5], #0
        'dir_main'          : config['SETUP']['DIR_MAIN'],  #'/discover/nobackup/projects/..'
        'fcst_name'         : str(config['BCSD']['fcst_data_type']),        
        'dataset_type'      : config['SETUP']['DATATYPE'], #'hindcast'
    }
    args['config'] = config

    return args


def create_daily_mean(ds):
    """Create daily mean from a days worth of 6-hourly data """
    return ds.resample(time='D').mean()

def _driver():
    """Main driver."""

    # Read in command arguments
    sys.stdout.write(f'Start. {pd.Timestamp.now().time()} \n')
    args = _read_cmd_args()

    # Input Args
    fcst_init_year  = int(args['fcst_init_year'])
    fcst_init_month = int(args['fcst_init_month'])
    ensemble_number = int(args['ensemble_number'])
    lead_month      = int(args['lead_month'])

    # Config Args
    dir_main     = args['dir_main']
    fcst_name    = args['fcst_name']
    dataset_type = args['dataset_type']

    # Netcdf attributes for output
    netcdf_attribute_dict = {
        'zlib': True, 'complevel': 6, 'shuffle': True, 'missing_value': -9999.,
        '_FillValue': -9999.}

    new_netcdf_attribute_dict = {
        var_name: netcdf_attribute_dict for var_name in get_all_hydrosfs_names()}

    # --- Generate range of lead days ---
    # Create the initial forecast day
    fcst_init_day = pd.to_datetime(f"{fcst_init_year}-{fcst_init_month}-01")

    # Calculate the start of the forecast month
    fcst_start_day = fcst_init_day + pd.DateOffset(months = lead_month)

    # Calculate the end of the forecast month
    fcst_end_day = fcst_start_day + pd.offsets.MonthEnd(0)

    # Calculate the range of days of the forecast month
    fcst_days = pd.date_range(start = fcst_start_day, end = fcst_end_day, freq='D')

    # Calculate the lead days from the initial forecast day
    lead_days = (fcst_days - fcst_init_day).days


    # Get 6-hourly filenames
    sys.stdout.write(f'Generating filenames. {pd.Timestamp.now().time()} \n')
    file_path_6hourly_in = []
    for lead_day in lead_days:
        file_path_6hourly_in.append(generate_raw_6hourly_forecast_filename(
            dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number,
            lead_day))

    # Read 6-hourly data and preprocess to daily means
    sys.stdout.write(f'Reading Data. {pd.Timestamp.now().time()} \n')
    ds = xr.open_mfdataset(file_path_6hourly_in, preprocess = create_daily_mean, concat_dim="time",
        combine="nested")

    # Create monthly mean
    sys.stdout.write(f'Taking Mean. {pd.Timestamp.now().time()} \n')
    ds = ds.mean(dim = 'time', keep_attrs=True)

    # Create output file directory
    file_directory_monthly_out = generate_raw_monthly_forecast_filename(
        dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number,
        lead_month, 'file_directory')
    if not os.path.exists(file_directory_monthly_out):
        os.makedirs(file_directory_monthly_out)

    # Write output file
    file_path_monthly_out = generate_raw_monthly_forecast_filename(
        dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number,
        lead_month)

    sys.stdout.write(f'Writing File. {pd.Timestamp.now().time()} \n')
    ds.to_netcdf(file_path_monthly_out, format='NETCDF4', mode = 'w',
                               encoding = new_netcdf_attribute_dict)


    sys.stdout.write(f'Wrote File. {pd.Timestamp.now().time()} \n')

#
# Main Method
#

if __name__ == '__main__':
    _driver()
