#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
#------------------------------------------------------------------------------
#
# SCRIPT: temporal_aggregate.py
#
# PURPOSE: Read daily S2S CF-convention netCDF files, calculate monthly
# averages and accumulations, and write to new CF-convention netCDF file.
# This version uses xarray for improved performance over netCDF4.
#
# REQUIREMENTS as of 2025:
# * Python 3.9 or higher
# * xarray library
# * netCDF4 library (for xarray backend)
# * numpy
# * yaml
#
# REVISION HISTORY:
# 2025: Converted from netCDF4 to xarray for ~3-5x performance improvement
#------------------------------------------------------------------------------
"""

# Standard modules
import os
import sys
import datetime
import time
import yaml

# Third-party libraries
import xarray as xr
import numpy as np

# Private methods - keep unchanged utility functions
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} configfile input_dir output_dir"
    txt += " fcst_yyyymmdd start_yyyymmdd end_yyyymmdd model_forcing"
    print(txt)
    print("[INFO] where:")
    print("[INFO]  configfile: path to s2spost config file")
    print("[INFO]  input_dir: directory with daily S2S files in CF convention")
    print("[INFO]  output_dir: directory for output file")
    print("[INFO]  fcst_yyyymmdd: Initial forecast date/time of daily files")
    print("[INFO]  start_yyyymmdd: Starting lead date/time of daily files")
    print("[INFO]  end_yyyymmdd: Ending lead date/time of daily files")
    print("[INFO]  model_forcing: ID for atmospheric forcing for LIS")

def _read_cmd_args(argv, logger, subtask):
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(argv) != 7:
        logger.error("Invalid number of command line arguments!", subtask=subtask)
        _usage()
        sys.exit(1)

    # Check if config file exists.
    configfile = argv[0]
    if not os.path.exists(configfile):
        logger.error(f"Directory {configfile} does not exist!", subtask=subtask)
        _usage()
        sys.exit(1)

    # Check if input directory exists.
    input_dir = argv[1]
    if not os.path.exists(input_dir):
        logger.error(f"Directory {input_dir} does not exist!", subtask=subtask)
        _usage()
        sys.exit(1)

    # Create output directory if it doesn't exist.
    output_dir = argv[2]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get valid initial forecast date.
    fcst_yyyymmdd = argv[3]
    fcstdate = _proc_date(fcst_yyyymmdd)

    # Get valid lead starting and ending dates of data.
    start_yyyymmdd = argv[4]
    startdate = _proc_date(start_yyyymmdd)
    end_yyyymmdd = argv[5]
    enddate = _proc_date(end_yyyymmdd)
    if startdate > enddate:
        logger.error(f"Start date is after end date!", subtask=subtask)
        sys.exit(1)

    # Get ID for model forcing for LIS
    model_forcing = argv[6].upper()

    return configfile, input_dir, output_dir, fcstdate, startdate, enddate, model_forcing

def _make_varlists(config):
    """Build lists of variables."""
    varlists = {}
    varlists["var_acc_list"] = config["POST"]["var_acc_list"]
    varlists["var_tavg_land_list"] = \
        config["POST"]["var_tavg_land_list"]
    varlists["var_tavg_f_list"] = \
        config["POST"]["var_tavg_f_list"]
    varlists["var_tavg_twsgws_list"] = \
        config["POST"]["var_tavg_twsgws_list"]
    varlists["var_tair_max_list"] = \
        config["POST"]["var_tair_max_list"]
    varlists["var_tair_min_list"] = \
        config["POST"]["var_tair_min_list"]
    varlists["const_list"] = config["POST"]["const_list"]
    return varlists

def _proc_date(yyyymmdd):
    """Convert YYYYMMDD string to Python date object."""
    if len(yyyymmdd) != 8:
        print("[ERR] Invalid length for YYYYMMDD, must be 8 characters!")
        sys.exit(1)
    year = int(yyyymmdd[0:4])
    month = int(yyyymmdd[4:6])
    day = int(yyyymmdd[6:8])
    try:
        dateobj = datetime.date(year, month, day)
    except ValueError:
        print("[ERR] Invalid YYYYMMDD passed to script!")
        sys.exit(1)
    return dateobj

def _check_filename_size(name):
    """Make sure filename does not exceed 128 characters, per Air Force
    requirement."""
    if len(os.path.basename(name)) > 128:
        print("[ERR] Output file name is too long!")
        print(f"[ERR] {os.path.basename(name)} exceeds 128 characters!")
        sys.exit(1)

def _create_daily_s2s_filename(input_dir, fcstdate, curdate, model_forcing, domain):
    """Create path to daily S2S netCDF file."""
    name = f"{input_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing}"
    name += "_GR.C0P25DEG"

    if domain == 'AFRICOM':
        name += "_AR.AFRICA"
    if domain == 'GLOBAL':
        name += "_AR.GLOBAL"

    name += "_PA.ALL"
    name += f"_DD.{fcstdate.year:04d}{fcstdate.month:02d}{fcstdate.day:02d}"
    name += "_DT.0000"
    name += f"_FD.{curdate.year:04d}{curdate.month:02d}{curdate.day:02d}"
    name += "_DT.0000"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

def create_time_aggregated_s2s_filename(output_dir, fcstdate, startdate, enddate,
                                 model_forcing, domain):
    """Create path to monthly S2S netCDF file."""
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing}"
    name += "_GR.C0P25DEG"

    if domain == 'AFRICOM':
        name += "_AR.AFRICA"
    if domain == 'GLOBAL':
        name += "_AR.GLOBAL"

    name += "_PA.ALL"
    name += f"_DD.{fcstdate.year:04d}{fcstdate.month:02d}{fcstdate.day:02d}"
    name += "_DT.0000"
    name += f"_FP.{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
    name += f"-{enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

# New xarray-based functions replacing the netCDF4 versions
def _create_time_aggregated_file_xarray(varlists, input_dir, output_dir, fcstdate, startdate,
                                        enddate, model_forcing, config, logger, subtask):
    """
    This is good to time aggregate daily files
    1) For weekly files precise startdate and enddates (included) are passed as arguments 
    2) For monthly files dates included are Day02 of the month through Day01 of the next month, 
       while datestamps in filenames show day01 of the month to day01 of the next month except for the first forecast month
    """

    is_monthly = False
    curdate = startdate
    if (enddate - startdate).days > 25:
        is_monthly = True

    if is_monthly and (startdate.month != fcstdate.month):
        ''' the 1st day has been counted by the previous month already, so we exclude it'''
        curdate += datetime.timedelta(days=1)
        
    # Create list of all daily files
    daily_files = []
    while curdate <= enddate:
        infile = _create_daily_s2s_filename(input_dir, fcstdate, curdate, model_forcing,
                                          config["EXP"]["DOMAIN"])
        logger.info(f'Reading: {infile}', subtask=subtask)
        if os.path.exists(infile):
            daily_files.append(infile)
        else:
            logger.error(f"Missing file: {infile}", subtask=subtask)
        curdate += datetime.timedelta(days=1)
    
    if not daily_files:
        logger.error(f"No daily files found!", subtask=subtask)
        return None
        
    # Separate variables into different categories
    acc_vars = varlists["var_acc_list"]
    tavg_vars = []
    for listname in ["var_tavg_land_list", "var_tavg_f_list", 
                     "var_tavg_twsgws_list", "var_tair_max_list", "var_tair_min_list"]:
        tavg_vars += varlists[listname]
    
    const_vars = varlists["const_list"]
        
    # First, open just one file to get the constant variables
    ds_first = xr.open_dataset(daily_files[0])
    
    # Create monthly dataset starting with constant variables
    monthly_ds = xr.Dataset()
    
    # Copy constant variables (these don't have time dimension)
    for var in const_vars:
        if var in ds_first:
            monthly_ds[var] = ds_first[var]  
        else:
            logger.error(f"Constant variable {var} not found", subtask=subtask)
    
    # Copy coordinates from first file (except time)
    for coord in ds_first.coords:
        if coord != 'time':
            monthly_ds = monthly_ds.assign_coords({coord: ds_first.coords[coord]})
    
    ds_first.close()
    
    # Now open all files for time-varying variables only
    try:
        # Only open variables that have time dimension
        time_varying_vars = acc_vars + tavg_vars
        
        # Open with data_vars parameter to only load specific variables
        ds_all = xr.open_mfdataset(
            daily_files,
            concat_dim='time',
            combine='nested',
            compat='override',
            coords='minimal',
            parallel=True,
            data_vars=time_varying_vars  # Only load time-varying variables
        )
    except Exception as e:
        logger.error(f"Failed to open multiple files: {e}", subtask=subtask)
        return None
    
    # Process accumulation variables (sum over time)
    for var in acc_vars:
        if var in ds_all:
            monthly_ds[var] = ds_all[var].sum(dim='time', keepdims=True)
        else:
            logger.error(f"Variable {var} not found in files", subtask=subtask)
    
    # Process time-average variables (mean over time)
    for var in tavg_vars:
        if var in ds_all:
            monthly_ds[var] = ds_all[var].mean(dim='time', keepdims=True)
        else:
            logger.error(f"Variable {var} not found in files", subtask=subtask)
        
    # Convert date to datetime and then to the right format for netCDF
    end_datetime = datetime.datetime.combine(enddate, datetime.time())
    
    # Create time coordinate as a number (following the original files' approach)
    # Use the same reference time as the daily files
    reference_time = datetime.datetime.combine(startdate, datetime.time())
    
    # Calculate minutes since start date
    time_value = [(end_datetime - reference_time).total_seconds() / 60.0]
    
    monthly_ds = monthly_ds.assign_coords({
        'time': xr.DataArray(
            time_value,
            dims=['time'],
            attrs={
                'units': f"minutes since {startdate.strftime('%Y-%m-%d')} 00:00:00",
                'begin_date': startdate.strftime('%Y%m%d'),
                'begin_time': "000000",
                'calendar': 'standard',
                'axis': 'T',
                'bounds': 'time_bnds',
                'long_name': 'time'
            }
        )
    })
    
    # Copy global attributes from the multi-file dataset
    monthly_ds.attrs = ds_all.attrs.copy()
    
    # Create time bounds for the monthly period
    days_in_month = (enddate - startdate).days
    time_bnds_data = [[0, days_in_month * 24 * 60]]  # Start to end in minutes
    monthly_ds['time_bnds'] = xr.DataArray(
        time_bnds_data,
        dims=['time', 'nv'],
        coords={'time': monthly_ds.time}
    )
    
    # Update cell methods for different variable types
    _update_cell_methods_xarray(monthly_ds, varlists)
    
    # Update global attributes
    monthly_ds.attrs['history'] = f"created on date: {time.ctime()}"
    
    # Define encoding for efficient writing
    encoding = {}
    
    # Handle coordinates and data variables
    all_vars = list(monthly_ds.coords.keys()) + list(monthly_ds.data_vars.keys())
    
    for var in all_vars:
        if var in ['time', 'lat', 'lon', 'ensemble', 'soil_layer']:
            # Skip coordinate variables - let xarray handle them automatically
            continue
        elif var == 'time_bnds':
            encoding[var] = {'dtype': 'float32', 'zlib': True, 'complevel': 6}
        elif var in monthly_ds.data_vars and monthly_ds[var].dtype == 'float64':
            encoding[var] = {'dtype': 'float32', 'zlib': True, 'complevel': 6, 'shuffle': True, '_FillValue': -9999.0}
        elif var in monthly_ds.data_vars:
            encoding[var] = {'zlib': True, 'complevel': 6, 'shuffle': True, '_FillValue': -9999.0}
    
    # Write output file
    outfile = create_time_aggregated_s2s_filename(output_dir, fcstdate, startdate,
                                         enddate, model_forcing, config["EXP"]["DOMAIN"])
    
    try:
        monthly_ds.to_netcdf(outfile, format='NETCDF4', encoding=encoding)
    except Exception as e:
        logger.error(f"Failed to write output file: {e}", subtask=subtask)
        return None
    finally:
        # Clean up
        ds_all.close()
        monthly_ds.close()
    
    return outfile

def _update_cell_methods_xarray(ds, varlists):
    """Update cell_method attributes for variables in xarray dataset."""
    
    for var_name in ds.data_vars:
        if var_name in varlists["var_tavg_twsgws_list"]:
            ds[var_name].attrs['cell_methods'] = \
                "time: mean (interval: 1 day) area: point where land"
        elif var_name in varlists["var_tair_max_list"]:
            ds[var_name].attrs['cell_methods'] = \
                "time: mean (interval: 1 day comment: daily maxima)"
        elif var_name in varlists["var_tair_min_list"]:
            ds[var_name].attrs['cell_methods'] = \
                "time: mean (interval: 1 day comment: daily minima)"
        elif var_name in varlists["var_tavg_land_list"]:
            ds[var_name].attrs['cell_methods'] = \
                "time: mean (interval: 1 day comment: daily means) area: point where land"
        elif var_name in varlists["var_tavg_f_list"]:
            ds[var_name].attrs['cell_methods'] = \
                "time: mean (interval: 1 day comment: daily means)"
        elif var_name in varlists["var_acc_list"]:
            ds[var_name].attrs['cell_methods'] = \
                "time: sum (interval: 1 day comment: daily sums)"

def agg_driver(argv, logger, subtask):
    """Main driver using xarray approach."""
        
    # Get the directories and dates
    configfile, input_dir, output_dir, fcstdate, startdate, enddate, model_forcing = _read_cmd_args(argv, logger, subtask)
    
    # Load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)
    
    varlists = _make_varlists(config)
        
    # Process all files at once using xarray
    outfile = _create_time_aggregated_file_xarray(
        varlists, input_dir, output_dir, fcstdate, 
        startdate, enddate, model_forcing, config,
        logger, subtask
    )
        
    if outfile:
        logger.info(f'Temporal aggregation complete: {outfile}', subtask=subtask)
    else:
        logger.error(f"Temporal aggregation failed", subtask=subtask)
        sys.exit(1)


