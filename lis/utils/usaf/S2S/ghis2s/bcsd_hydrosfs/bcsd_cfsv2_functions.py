#!/usr/bin/env python3
"""
Functions to get file information from the S2S E2ES bcsd_fcst substep

Script
------
bcsd_cfsv2_functions.py

Requirements
------------
Python 3.9

Revision History
----------------
07 Mar 2024: Ryan Zamora, first version

"""
#
# Modules
#
import sys
from datetime import datetime
import numpy as np
import xarray as xr
import cfgrib

#
# Custom Modules
#
from bcsd_filename_functions import generate_file_information
from bcsd_helper_functions import create_lat_lon_xr_dataset, create_regridder_weights_file, apply_regridder_weights_file


#
# Variables
#
_refor_hps_varnames = ['dlwsfc', 'dswsfc', 'q2m', 'wnd10m']
_refor_fl_varnames = ['prate', 'pressfc', 'tmp2m']

# PSEUDO INPUTS
# TODO This naming convention should change in the future; currently only used in raw data

# TODO Should lead_number be renamed to lead_month?

# ------------------------- #
# - VALID VALUE FUNCTIONS - #
# ------------------------- #

def _check_valid_fcst_init_year(value : int):
    """ Checks for valid fcst_init_year value"""
    if not isinstance(value, int):
        raise TypeError('fcst_init_year must be an integer')
    if value not in range(1000,10000):
        raise ValueError('fcst_init_year must be a 4-digit integer')


def _check_valid_fcst_init_month(value : int):
    """Checks for valid fcst_init_month value"""
    if not isinstance(value, int):
        raise TypeError('fcst_init_month must be an integer')
    if value not in range(1,13):
        raise ValueError('fcst_init_month must be an integer with value 1-12')

def _check_valid_ensemble_number(value : int):
    """Checks for valid ensemble_number value"""
    if not isinstance(value, int):
        raise TypeError('ensemble_number must be an integer')
    if value <= 0:
        raise ValueError('ensemble_number must be a positive non-zero integer')


def _check_valid_lead_number(value : int):
    """Checks for valid lead_number value"""
    if not isinstance(value, int):
        raise TypeError('lead_number must be an integer')
    if value < 0:
        raise ValueError('lead_number must be a positive integer')


def _check_valid_download_variable_name(value : str):
    """Checks for valid download_variable_name value"""
    valid_download_variable_names = _refor_hps_varnames + _refor_fl_varnames
    if not isinstance(value, str):
        raise TypeError('download_variable_name must be a string')
    if value not in valid_download_variable_names:
        raise ValueError('download_variable_name must be a string with value: ' +
                         str(valid_download_variable_names))


def _check_valid_observation_variable_name(value :str):
    """Checks for valid observation_variable_name value"""
    valid_observation_variable_names = ['LWS', 'PRECTOT', 'PS', 'Q2M', 'SLRSF', 'T2M', 'WIND10M']
    if not isinstance(value, str):
        raise TypeError('observation_variable_name must be a string')
    if value not in valid_observation_variable_names:
        raise ValueError('observation_variable_name must be a string with value: ' +
                         str(valid_observation_variable_names))


def _check_valid_dataset_type(value : str):
    """Checks for valid dataset_type value"""
    valid_dataset_types = ['hindcast', 'forecast']
    if not isinstance(value, str):
        raise TypeError('dataset_type must be a string')
    if value not in valid_dataset_types:
        raise ValueError(f"dataset_type must be a string with value: {valid_dataset_types}")


# -------------------- #
# - HELPER FUNCTIONS - #
# -------------------- #

def get_downloaded_cfsv2_date_strings(
    fcst_init_year : int, fcst_init_month : int,
    ensemble_number : int) -> tuple[str, str, str, str]:
    """
    Returns the associated downloaded cfsv2 date strings

    Parameters
    ----------
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    ensemble_number (int) : Ensemble number

    Returns
    -------
    year  (str) : Year of the downloaded cfsv2 forecast
    month (str) : Month of the downloaded cfsv2 forecast
    day   (str) : Day of the downloaded cfsv2 forecast
    hour  (str) : Hour of the downloaded cfsv2 forecast

    Examples
    --------
    >>> year, month, day, hour = get_downloaded_cfsv2_date_strings(1991, 5, 7)
    >>> year
    1991
    >>> month
    '04'
    >>> day
    '21'
    >>> hour
    '12'
    """

    # Check for valid parameters
    _check_valid_fcst_init_year(fcst_init_year)
    _check_valid_fcst_init_month(fcst_init_month)
    _check_valid_ensemble_number(ensemble_number)

    # Download Hours
    fcst_download_hours = ['00', '06', '12', '18']

    # Download Monthdays (MMDD)
    fcst_download_monthdays = [
        ['1217', '1222', '1227'],
        ['0121', '0126', '0131'],
        ['0215', '0220', '0225'],
        ['0317', '0322', '0327'],
        ['0416', '0421', '0426'],
        ['0521', '0526', '0531'],
        ['0620', '0625', '0630'],
        ['0720', '0725', '0730'],
        ['0819', '0824', '0829'],
        ['0918', '0923', '0928'],
        ['1018', '1023', '1028'],
        ['1117', '1122', '1127']
    ]

    # Download Year
    year = fcst_init_year if fcst_init_month != 1 else fcst_init_year - 1

    # Download Monthday (MMDD)
    index_monthday = np.repeat(a=[0,1,2],repeats=4)[ensemble_number - 1]
    monthday = (fcst_download_monthdays[fcst_init_month - 1] * 4)[index_monthday]

    # Download Month
    month = monthday[0:2]

    # Download Day
    day = monthday[2:4]

    # Download Hour
    index_hour = (ensemble_number - 1) % 4
    hour = fcst_download_hours[index_hour]

    return year, month, day, hour


# ------------------------------ #
# - DOWNLOADED CFSv2 FUNCTIONS - #
# ------------------------------ #

def generate_downloaded_cfsv2_filename(
    fcst_init_year : int, fcst_init_month : int, ensemble_number : int, variable_name : str,
    dir_downloads : str = '/discover/nobackup/projects/lis/MET_FORCING/CFSv2',
    output_type: str = 'file_path') -> str:
    """
    Returns the associated downloaded cfsv2 file information

    Parameters
    ----------
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    ensemble_number (int) : Ensemble number
    variable_name   (str) : Variable name (download)
    dir_downloads   (str) : Directory of cfsv2 downloads
    output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']

    Returns
    -------
    file_information (str) : File information of downloaded cfsv2 dataset for given arguments

    Examples
    --------
    >>> x = generate_downloaded_cfsv2_filename(1991, 5, 7, 'prate', '~/CFSv2')
    >>> x
    '~/CFSv2/Refor_FL/1991/19910421/prate.1991042112.time.grb2'
    >>> y = generate_downloaded_cfsv2_filename(1991, 5, 7, 'wnd10m', '~/CFSv2')
    >>> y
    '~/CFSv2/Refor_HPS/1991/19910421/wnd10m.1991042112.time.grb2'
    """

    # Check for valid parameters
    _check_valid_fcst_init_year(fcst_init_year)
    _check_valid_fcst_init_month(fcst_init_month)
    _check_valid_ensemble_number(ensemble_number)
    _check_valid_download_variable_name(variable_name)
    # _check_valid_output_type(output_type)

    # Cutoff Date for Refor_HPS & Refor_FL Datasets
    refor_datetime = datetime(2011,4,1)

     # Get download _monthday and fcst_ic_hour
    year, month, day, hour = get_downloaded_cfsv2_date_strings(
        fcst_init_year, fcst_init_month, ensemble_number)

    # Download Datetime
    download_datetime = datetime.strptime(f"{year}{month}{day}{hour}", '%Y%m%d%H')

    # Determine filename attributes
    if download_datetime < refor_datetime:
        if variable_name in _refor_hps_varnames:
            subdir = 'Refor_HPS'
        else:
            subdir = 'Refor_FL'

        file_pfx = variable_name
        file_sfx = 'time.grb2'
    else:
        subdir = 'Oper_TS'
        file_pfx = f"{variable_name}.01"
        file_sfx = 'daily.grb2'

    # File directory and name information
    file_directory = f"{dir_downloads}/{subdir}/{year}/{year}{month}{day}"
    file_name = f"{file_pfx}.{year}{month}{day}{hour}.{file_sfx}"

    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)

    return file_information

def create_cfsv2_elevation_difference_file(dir_supplementary) -> xr.DataArray:

    # Create ds_out grid
    lon_start = -179.975
    lon_end = 179.975
    lat_start = -59.975
    lat_end = 89.975
    cell_step = 0.05
    
    ds_out = create_lat_lon_xr_dataset(lon_start, lon_end, lat_start, lat_end, cell_step)
    
    # CFSv2 Terrain Height Field
    file_path_cfsv2 = f"{dir_supplementary}/sample/cfsv2_staticfields_2015040100.01.grb2"
    ds_cfsv2 = cfgrib.open_dataset(file_path_cfsv2, indexpath ='')[['orog']]
    variable_dict = {'orog': {'method':'bilinear'}}

    ds_cfsv2 = apply_regridder_weights_file(ds_cfsv2, ds_out, variable_dict, dir_supplementary)
    
    # MERIT Elevation
    file_path_merit = f"{dir_supplementary}/lis_darun/lis_input.global_5km_hydroscs_mask_topo_merit1km_wclimppt.nc"
    ds_merit = xr.open_dataset(file_path_merit)
    
    
    # Compute Difference between the two
    da_elevdiff = xr.DataArray(
        ds_merit['ELEVATION'].values - ds_cfsv2['orog'].values,
        dims = ("lat", "lon"), coords = {"lat": ds_cfsv2.lat, "lon": ds_cfsv2.lon})
    
    return da_elevdiff

def read_grib2_file(filename, varname = None):
    ''' reads grib2 file and renames original names '''
    
    # New Variable Names
    converted_variable_names = {
        'prate':'PRECTOT', 'sp':'PS', 't2m':'T2M',
        'dlwrf':'LWS', 'dswrf':'SLRSF',
        'u10':'U10M', 'v10':'V10M',
        'q':'Q2M', 'sh2':'Q2M'
    }

    # Open file
    ds_out = cfgrib.open_dataset(filename, indexpath ='')
    
    # Rename time variable
    ds_out = ds_out.rename_vars({"time": "reference_time"})
    ds_out = ds_out.rename_vars({'valid_time': 'time'})
    ds_out = ds_out.swap_dims({'step': 'time'})
    ds_out = ds_out.drop_vars(['heightAboveGround','surface'], errors='ignore')
    # Rename variables
    for name, da_ in ds_out.data_vars.items():
        ds_out = ds_out.rename({name : converted_variable_names.get(name)})
    
    # Special case for wind components
    if varname in ['u10', 'v10']:
        ds_out = ds_out[[converted_variable_names.get(varname)]]
    
    # Special case for wind magnitude
    if varname == 'wind10':
        wnd10m = magnitude(ds_out.U10M, ds_out.V10M)
        ds_out['WIND10M'] = wnd10m
        ds_out['WIND10M'].attrs['units'] = 'm/s'
        ds_out['WIND10M'].attrs['short_name'] = 'wnd10m'
        ds_out['WIND10M'].attrs['long_name'] = 'Wind Speed'
        ds_out['WIND10M'].attrs['level'] = '10 m above ground'
        ds_out = ds_out[['WIND10M']]
    
    return ds_out

def create_cfsv2_5km_regridder_weights_file(method : str, dir_out : str):
    """
    Writes regridder weights to file

    Parameters
    ----------
    method  (str) : Regridder method
    dir_out (str) : Directory path of output file

    Examples
    --------
    """
    # Wanted Grid
    lon_start = -179.975
    lon_end = 179.975
    lat_start = -59.975
    lat_end = 89.975
    cell_step = 0.05
    
    # Example File to use
    ds_in = cfgrib.open_dataset(
        generate_downloaded_cfsv2_filename(2020, 5, 1, "prate"), indexpath='')
    
    # Create ds_out grid
    ds_out = create_lat_lon_xr_dataset(lon_start, lon_end, lat_start, lat_end, cell_step)
    
    # Write to file
    create_regridder_weights_file(ds_in, ds_out, method, dir_out)


# ------------------------------- #
# - WRAPPER INVENTORY FUNCTIONS - #
# ------------------------------- #

# def inventory_downloaded_cfsv2_filename(fcst_init_year : int, fcst_init_month,
# dir_downloads : str = '/discover/nobackup/projects/lis/MET_FORCING/CFSv2',
# output_type: str = 'file_path')  -> list[str]:
#     """
#     Returns all downloaded cfsv2 file information for a given forecast year and initial month

#     Parameters
#     ----------
#     fcst_init_year  (int) : Initial forecast year
#     fcst_init_month (int) : Initial forecast month
#     dir_downloads   (str) : Directory of cfsv2 downloads
#     output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']

#     Returns
#     -------
#     file_information (Array[str]) : File information of downloaded cfsv2 dataset for given arguments

#     Examples
#     --------

#     """

#     # Loop over variable_name and ensemble_number
#     generate_downloaded_cfsv2_filename(
#         fcst_init_year, fcst_init_month, ensemble_number, variable_name, dir_downloads, output_type)
    