#!/usr/bin/env python3
""" Creates a variable forecast climatology from the monthly averaged files
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
    generate_raw_monthly_forecast_filename,
    generate_raw_climatology_forecast_filename,
)

#
# Functions
#
def _usage():
    """Print command line usage."""
    txt =  f'[INFO] Usage: {sys.argv[0]}'
    txt += 'config_filename variable_name fcst_init_month'
    print(txt)

def _read_cmd_args():
    """Read command line arguments"""

    with open(sys.argv[1], "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    args = {
        "config_filename" : str(sys.argv[1]),
        "fcst_init_month" : int(sys.argv[2]),
        "variable_name"   : str(sys.argv[3]),
        "dir_main"        : str(config["SETUP"]["DIR_MAIN"]),
        "fcst_name"       : str(config["BCSD"]["fcst_data_type"]),
        "dataset_type"    : str(config["SETUP"]["DATATYPE"]),
        "clim_syr"        : int(config["BCSD"]["clim_start_year"]),
        "clim_eyr"        : int(config["BCSD"]["clim_end_year"]),
        "num_lead_months" : int(config['EXP']['lead_months']),
        "num_ensembles"   : int(config['BCSD']['nof_raw_ens']),
    }
    args["config"] = config

    return args

def preprocess_monthly_files(ds):
    """Selects the needed variable in the preprocessor, to save memory"""
    da = ds[var_name]
    # ens_num = int(ds.encoding['source'].split('ens')[1].split('/')[0])
    # da = da.expand_dims('ens').assign_coords({'ens':[ens_num]})
    return da
#
# Main
#
args = _read_cmd_args()

# Input Args
fcst_name = args["fcst_name"]
fcst_init_month = args["fcst_init_month"]

# Config Args
dir_main = args["dir_main"]
dataset_type = args["dataset_type"]
var_name = args["variable_name"]
clim_syr = args["clim_syr"]
clim_eyr = args["clim_eyr"]
num_lead_months = args["num_lead_months"]
num_ensembles = args["num_ensembles"]

# Coordinate Information
lats, lons = get_domain_info(args["config_filename"], coord=True)
lat1, lat2, lon1, lon2 = get_domain_info(args["config_filename"], extent=True)

netcdf_attribute_dict = {
    "zlib":True, "complevel":6, "shuffle":True, "missing_value": np.nan, "_FillValue": np.nan}

# Defining array to store data
num_years = clim_eyr - clim_syr + 1
clim_array = np.empty((num_lead_months + 1, num_years * num_ensembles, len(lats), len(lons)))

# Create the initial quantile array
quant_ts_1d = np.arange(1, num_years * num_ensembles + 1) / (num_years * num_ensembles + 1)

# Reshape and broadcast the quantiels from (nyears*nens) to (nyears*nens, nlats, nlons)
quant_ts = np.broadcast_to(
    quant_ts_1d[:, np.newaxis, np.newaxis],
    (num_years * num_ensembles, len(lats), len(lons))
)

# Set the quantiles to index 0 of the time dimension
clim_array[0, :, :, :] = quant_ts

# Read monthly data for all lead months, fcst_init_years, and ensemble members
for lead_month in range(num_lead_months):

    # Generate filelist to read by month
    file_path_monthly_in = [
        generate_raw_monthly_forecast_filename(
            dir_main,
            fcst_name,
            dataset_type,
            fcst_init_year,
            fcst_init_month,
            ens_num,
            lead_month
        )
        for fcst_init_year in range(clim_syr, clim_eyr + 1)
        for ens_num in range(1, num_ensembles + 1)
    ]

    # Read all monthly files and select the needed variable
    ds_in = xr.open_mfdataset(
        file_path_monthly_in,
        concat_dim = "time",
        combine = "nested",
        preprocess = preprocess_monthly_files
    )
    clim_array[lead_month + 1, :, :, :] = np.sort(ds_in.to_numpy(), axis = 0)

# Create the output dataset
ds_out = xr.Dataset(
    data_vars = {
        "clim": (("DIST", "time", "latitude", "longitude"), clim_array),
    },
    coords = {
        "DIST": np.arange(num_lead_months + 1),
        "time": np.arange(num_years * num_ensembles),
        "latitude": lats,
        "longitude": lons,
    }
)

# Select only the data for the lat/lon box
ds_out = ds_out.sel(longitude = slice(lon1, lon2), latitude = slice(lat1, lat2))

# Create output file directory
file_directory_climatology_out = generate_raw_climatology_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_month, var_name, "file_directory")
if not os.path.exists(file_directory_climatology_out):
    os.makedirs(file_directory_climatology_out)

# Write output file
enc_dict = {var: netcdf_attribute_dict for var in ds_out.data_vars}

file_path_climatology_out = generate_raw_climatology_forecast_filename(
    dir_main, fcst_name, dataset_type, fcst_init_month, var_name)

ds_out.to_netcdf(file_path_climatology_out, format = "NETCDF4", mode = 'w', encoding = enc_dict)
