#!/usr/bin/env python3
"""
Functions to get file information for GEOS-S2S-V3 input files

Script
------
bcsd_geosv3_functions.py

Requirements
------------
Python 3.9

Revision History
----------------
08 May 2025: Ryan Zamora, first version

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

from bcsd_helper_functions import create_lat_lon_xr_dataset, create_regridder_weights_file, apply_regridder_weights_file

# Custom Modules
from bcsd_filename_functions import generate_file_information


# ------------- #
# - VARIABLES - #
# ------------- #

# 3-Hourly Surface File
def generate_geosv3_3hr_surface_filename(
    fcst_init_year : int, fcst_init_month : int, ensemble_number : int, lead_number : int,
    dir_geos : str, output_type: str = 'file_path') -> str:
    """
    Returns the associated input geos-s2s-v3 surface file information

    Parameters
    ----------
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    ensemble_number (int) : Ensemble number
    lead_number     (int) : Lead number    
    dir_geos        (str) : Directory for geos-s2s-v3 data
    output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']    

    Returns
    -------
    file_information (str) : File information of downloaded cfsv2 dataset for given arguments

    Examples
    --------

    """
    
    # Generate datetime of forecast
    fcst_datetime = datetime(fcst_init_year, fcst_init_month,1) + relativedelta(months=lead_number)

    # Formatted forecast year
    fcst_year  = fcst_datetime.strftime('%Y')

    # Formatted forecast month
    fcst_month = fcst_datetime.strftime('%m')
    
    # File directory and name information
    file_directory = f"{dir_geos}/sfc_tavg_3hr_glo_L720x361_sfc/{fcst_init_year:04d}{fcst_init_month:02d}/ens{ensemble_number:02d}"
    file_name = f"geos_s2s_v3.{fcst_year}{fcst_month}.nc"

    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information


# 3-Hourly Radiation File
def generate_geosv3_3hr_radiation_filename(
    fcst_init_year : int, fcst_init_month : int, ensemble_number : int, lead_number : int,
    dir_geos : str, output_type: str = 'file_path') -> str:
    """
    Returns the associated input geos-s2s-v3 radiation file information

    Parameters
    ----------
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    ensemble_number (int) : Ensemble number
    lead_number     (int) : Lead number    
    dir_geos        (str) : Directory for geos-s2s-v3 data
    output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']    

    Returns
    -------
    file_information (str) : File information of downloaded cfsv2 dataset for given arguments

    Examples
    --------

    """
    
    # Generate datetime of forecast
    fcst_datetime = datetime(fcst_init_year, fcst_init_month,1) + relativedelta(months=lead_number)

    # Formatted forecast year
    fcst_year  = fcst_datetime.strftime('%Y')

    # Formatted forecast month
    fcst_month = fcst_datetime.strftime('%m')
    
    # File directory and name information
    file_directory = f"{dir_geos}/sfc_tavg_3hr_glo_L720x361_sfc/{fcst_init_year:04d}{fcst_init_month:02d}/ens{ensemble_number:02d}"
    file_name = f"rad_geos_s2s_v3.{fcst_year}{fcst_month}.nc"

    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

# 3-Hourly Albedo File
def generate_geosv3_3hr_albedo_filename(
    fcst_month : int, dir_geos : str, output_type: str = 'file_path') -> str:
    """
    Returns the associated input geos-s2s-v3 albedo file information

    Parameters
    ----------
    fcst_month  (int) : Forecast month
    dir_geos    (str) : Directory for geos-s2s-v3 data
    output_type (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']    

    Returns
    -------
    file_information (str) : File information of downloaded cfsv2 dataset for given arguments

    Examples
    --------

    """
    
    # File directory and name information
    file_directory = f"{dir_geos}/albedo_avg"
    file_name = f"albedo.{fcst_month:02d}all.nc4"

    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

def create_geosv3_5km_regridder_weights_file(method : str, dir_out : str, dir_supplementary : str):
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
    file_path_geosv3 = f"{supplementary_dir}/sample/geosv3_samplefile.nc"
    ds_in = xr.open_dataset(file_path_geosv3)
    
    # Create ds_out grid
    ds_out = create_lat_lon_xr_dataset(lon_start, lon_end, lat_start, lat_end, cell_step)
    
    # Write to file
    create_regridder_weights_file(ds_in, ds_out, method, dir_out)

def create_geosv3_elevation_difference_file(dir_supplementary : str) -> xr.DataArray:

    # Gravity Constant
    g_const = 9.80665

    # Create ds_out grid
    lon_start = -179.975
    lon_end = 179.975
    lat_start = -59.975
    lat_end = 89.975
    cell_step = 0.05

    ds_out = create_lat_lon_xr_dataset(lon_start, lon_end, lat_start, lat_end, cell_step)
    
    # GEOSv3 Terrain Height Field
    file_path_geosv3 = f"{dir_supplementary}/sample/geosv3_staticfields_20150401.nc"
    ds_geosv3 = xr.open_dataset(file_path_geosv3)[['PHIS']].isel(time = 0)
    ds_geosv3['PHIS'].values = ds_geosv3['PHIS'].values/g_const
    
    variable_dict = {'PHIS': {'method':'bilinear'}}

    ds_geosv3 = apply_regridder_weights_file(ds_geosv3, ds_out, variable_dict, dir_supplementary)
    
    # MERIT Elevation
    file_path_merit = f"{dir_supplementary}/lis_darun/lis_input.global_5km_hydroscs_mask_topo_merit1km_wclimppt.nc"
    ds_merit = xr.open_dataset(file_path_merit)
    
    # Compute Difference between the two
    da_elevdiff = xr.DataArray(
        ds_merit['ELEVATION'].values - ds_geosv3['PHIS'].values,
        dims = ("lat", "lon"), coords = {"lat": ds_geosv3.lat, "lon": ds_geosv3.lon})
    
    return da_elevdiff