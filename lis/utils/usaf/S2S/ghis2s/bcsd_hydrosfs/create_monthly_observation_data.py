#!/usr/bin/env python3
"""
Create Monthly Observation Data

Workflow
--------

"""

#
# Modules
#
import sys
import os
from datetime import datetime
import math
import yaml
import numpy as np
import xarray as xr
import pandas as pd
from dateutil.relativedelta import relativedelta

#
# Custom Modules
#
from bcsd_filename_functions import (generate_downloaded_hydroscs_filename, generate_raw_monthly_hydroscs_filename)
from dict_variables import get_hydroscs_name, get_all_hydroscs_names, get_all_hydrosfs_names

#
# Functions
#
def _usage():
    """ Print command line usage """
    txt = f"[INFO] Usage: {sys.argv[0]} config_filename year month"
    print(txt)

def _read_cmd_args():
    """ Read command line arguments """

    with open(sys.argv[1], "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    args = {
        "config_filename" : str(sys.argv[1]),
        "year"            : int(sys.argv[2]),
        "month"           : int(sys.argv[3]),
        "dir_main"        : str(config["SETUP"]["DIR_MAIN"]),
        "dir_hydroscs"    : str(config["SETUP"]["DIR_HYDROSCS"]),
        "fcst_name"       : str(config["BCSD"]["fcst_data_type"]),
        "dataset_type"    : str(config["SETUP"]["DATATYPE"]),
    }

    args["config"] = config

    return args

def _driver():
    """ Main driver """
    args = _read_cmd_args()
    
    # Config Args
    dir_main = args["dir_main"]
    dir_hydroscs = args["dir_hydroscs"]
    fcst_name    = args["fcst_name"]
    dataset_type = args["dataset_type"]
    year = args["year"]
    month = args["month"]
    
    netcdf_attribute_dict = {
        "zlib": True, "complevel": 6, "shuffle": True, "missing_value": -9999., "_FillValue": -9999.}
    
    new_netcdf_attribute_dict = {
        var_name: netcdf_attribute_dict for var_name in get_all_hydrosfs_names()}
    
    # Generate datetime information
    dt_daily_begin = datetime(year, month, 1)
    dt_daily_end   = dt_daily_begin + relativedelta(months = 1, days = -1)
    dt_daily_range = pd.date_range(start = dt_daily_begin, end = dt_daily_end, freq = "D")
    
    ds_out = xr.Dataset()
    
    for start_of_day in dt_daily_range:
        print(start_of_day)
        end_of_day = start_of_day + relativedelta(days = 1, hours = -1)
    
        dt_range = pd.date_range(start = start_of_day, end = end_of_day, freq = "h")
    
        # Get 6-hourly filenames
        file_path_hourly_in = [generate_downloaded_hydroscs_filename(dir_hydroscs, year, month, day_time.day, day_time.hour) for day_time in dt_range]
    
        # Read 6-hourly data
        time_name = "new_time"
        time_value = np.datetime64(f"{year:04d}-{month:02d}-01")
        ds_in = xr.open_mfdataset(file_path_hourly_in, concat_dim = time_name, combine = "nested")
    
        # Create daily mean
        ds_in = ds_in.mean(dim = time_name).assign_coords({"time": time_value}).expand_dims("time")
    
        if start_of_day == dt_daily_range[0]:
            ds_out = ds_in.copy()
        else:
            ds_out = xr.concat([ds_out, ds_in.copy()], dim = "time")
    
    # Daily to Monthly Mean
    ds_out = ds_out.mean(dim = "time")

    # Correction to HydroSCS Data
    print("[INFO] In the HydroSCS dataset (as of 09/01/2025) U-wind is multiplied by -1.")
    print("[INFO] This script is multiplying this data by -1 to correct it.")
    print("[INFO] If in the future, the dataset is corrected, remove the proceeding line.")
    print("ds_out[get_hydroscs_name("U")].values = -1*ds_out[get_hydroscs_name("U")].values")
    ds_out[get_hydroscs_name("U")].values = -1*ds_out[get_hydroscs_name("U")].values
    
    # Rename variables to final naming convention (Naming conventions are set in dict_variables.py)
    variable_dict = dict(zip(get_all_hydroscs_names(), get_all_hydrosfs_names()))
    ds_out = ds_out.rename(variable_dict)
    
    # Create output file directory
    file_directory_monthly_out = generate_raw_monthly_hydroscs_filename(
        dir_main, dataset_type, year, month, "file_directory")
    if not os.path.exists(file_directory_monthly_out):
        os.makedirs(file_directory_monthly_out)
    
    # Write output file
    file_path_monthly_out = generate_raw_monthly_hydroscs_filename(
        dir_main, dataset_type, year, month)
    
    ds_out.load().to_netcdf(file_path_monthly_out, format="NETCDF4", mode = "w",
                               encoding = new_netcdf_attribute_dict)


#
# Main
#
if __name__ == "__main__":
    _driver()
