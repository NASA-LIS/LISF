#!/usr/bin/env python3
""" Creates a variable observation climatology from the monthly averaged files
1. Generate quantile array
2. Create sorted time series
"""

#
# Modules
#
import os
import sys
import numpy as np
import xarray as xr
import yaml

#
# Custom Modules
#
from bcsd_helper_functions import get_domain_info
from bcsd_filename_functions import (
    generate_raw_monthly_hydroscs_filename,
    generate_raw_climatology_hydroscs_filename
)

#
# Functions
#
def _usage():
    """Print command line usage."""
    txt =  f'[INFO] Usage: {sys.argv[0]}'
    txt += 'config_filename variable_name'
    print(txt)

def _read_cmd_args():
    """Read command line arguments"""

    with open(sys.argv[1], "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 3:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    args = {
        "config_filename" : str(sys.argv[1]),
        "variable_name"   : str(sys.argv[2]),
        "dir_main"        : str(config["SETUP"]["DIR_MAIN"]),
        "dataset_type"    : str(config["SETUP"]["DATATYPE"]),
        "clim_syr"        : int(config["BCSD"]["clim_start_year"]),
        "clim_eyr"        : int(config["BCSD"]["clim_end_year"]),
    }
    args["config"] = config

    return args


#
# Main
#
args = _read_cmd_args()

# Input Args
variable_name  = str(args["variable_name"])

# Config Args
dir_main = args["dir_main"]
dataset_type = args["dataset_type"]
VAR_NAME = args["variable_name"]
CLIM_SYR = args["clim_syr"]
CLIM_EYR = args["clim_eyr"]

# Coordinate Information
LATS, LONS = get_domain_info(args["config_filename"], coord=True)
lat1, lat2, lon1, lon2 = get_domain_info(args["config_filename"], extent=True)

netcdf_attribute_dict = {
    "zlib": True, "complevel": 6, "shuffle": True, "missing_value": -9999., "_FillValue": -9999.}


# Defining array to store data
NYEARS = CLIM_EYR - CLIM_SYR + 1
CLIM_ARRAY = np.empty((13, NYEARS, len(LATS), len(LONS)))

# Create the initial quantile array
quant_ts_1d = np.arange(1, NYEARS + 1) / (NYEARS + 1)

# Reshape and broadcast the quantiels from (nyears) to (nyears, nlats, nlons)
quant_ts = np.broadcast_to(
    quant_ts_1d[:, np.newaxis, np.newaxis],
    (NYEARS, len(LATS), len(LONS))
)

# Set the quantiles to index 0 of the "month" dimension
CLIM_ARRAY[0, :, :, :] = quant_ts

# Read monthly data for all months and years
for month in range(1,13):

    # Generate filelist to read by month
    file_path_monthly_in = [
        generate_raw_monthly_hydroscs_filename(
            dir_main,
            dataset_type,
            year,
            month
        )
        for year in range(CLIM_SYR, CLIM_EYR + 1)
    ]

    # Read all monthly files and select the needed variable
    ds_in = xr.open_mfdataset(
        file_path_monthly_in, concat_dim = "time", combine = "nested")[VAR_NAME]

    CLIM_ARRAY[month, :, :, :] = np.sort(ds_in.to_numpy(), axis = 0)

# Create the output dataset
ds_out = xr.Dataset(
    data_vars = {
        "clim": (("DIST", "time", "latitude", "longitude"), CLIM_ARRAY),
    },
    coords = {
        "DIST": np.arange(13),
        "time": np.arange((CLIM_EYR - CLIM_SYR) + 1),
        "latitude": LATS,
        "longitude": LONS,
    }
)

# Select only the data for the lat/lon box
ds_out = ds_out.sel(longitude = slice(lon1, lon2), latitude = slice(lat1, lat2))

# Create output file directory
file_directory_climatology_out = generate_raw_climatology_hydroscs_filename(
    dir_main, dataset_type, VAR_NAME, "file_directory")
if not os.path.exists(file_directory_climatology_out):
    os.makedirs(file_directory_climatology_out)

# Write output file
file_path_climatology_out = generate_raw_climatology_hydroscs_filename(
    dir_main, dataset_type, VAR_NAME)

ds_out.to_netcdf(file_path_climatology_out, format="NETCDF4", mode = "w",
                           encoding = {"clim" : netcdf_attribute_dict})
