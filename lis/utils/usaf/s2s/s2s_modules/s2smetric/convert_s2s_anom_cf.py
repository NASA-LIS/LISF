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
# SCRIPT: convert_s2s_anom_cf.py
#
# PURPOSE: Updates attributes of S2S monthly anomaly/standardized anomaly
# netCDF files to conform with CF-1.8 convention.
#
# REQUIREMENTS as of 28 May 2023:
# * Python 3.9 or higher
#
# REFERENCES:
# https://cfconventions.org for specifications of NetCDF Climate and Forecast
#   (CF) Metadata Conventions.
# https://unidata.ucar.edu/software/udunits for documentation on UDUNITS2
#   library, which CF is generally consistent with for unit specifications.
#
# REVISION HISTORY:
# 26 Sep 2021: Eric Kemp (SSAI), first version.
# 30 Oct 2021: Eric Kemp/SSAI, now uses s2smetric config file.
# 02 Jun 2023: K. Arsenault + S. Mahanama, Updated the 557 WW file names.
#
#------------------------------------------------------------------------------
"""


# Standard modules

import os
import shutil
import sys
from netCDF4 import Dataset #pylint: disable=no-name-in-module

# Local constants
# Units for variable anomalies.  Standardized anomalies will be dimensionless.
_UNITS_ANOM = {
    "RZSM" : "m3 m-3",
    "SFCSM" : "m3 m-3",
    "TWS" : "mm",
    "Precip" : "kg m-2",
    "AirT" : "K",
    "ET" : "kg m-2 s-1",
    "Streamflow" : "m3 s-1",
}
_LONG_NAMES_SANOM = {
    "RZSM" : "Root zone soil moisture standardized anomaly",
    "SFCSM" : "Surface soil moisture standardized anomaly",
    "TWS" : "Terrestrial water storage standardized anomaly",
    "Precip" : "Total precipitation amount standardized anomaly",
    "AirT" : "Air temperature standardized anomaly",
    "ET" : "Total evapotranspiration standardized anomaly",
    "Streamflow" : "Streamflow standardized anomaly",
}
_LONG_NAMES_ANOM = {
    "RZSM" : "Root zone soil moisture anomaly",
    "SFCSM" : "Surface soil moisture anomaly",
    "TWS" : "Terrestrial water storage anomaly",
    "Precip" : "Total precipitation amount anomaly",
    "AirT" : "Air temperature anomaly",
    "ET" : "Total evapotranspiration anomaly",
    "Streamflow" : "Streamflow anomaly",
}

def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} anom_file output_dir"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]  anom_file: Name of netCDF file with anom or sanom metric")
    print("[INFO]  output_dir: Directory to write CF-convention file")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 3:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if netCDF file exists
    anom_filename = sys.argv[1]
    if not os.path.exists(anom_filename):
        print(f"[ERR] {anom_filename} does not exist!")
        sys.exit(1)

    # Create output directory if it doesn't exist.
    output_dir = sys.argv[2]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return anom_filename, output_dir

def _copy_anom_file(anom_filename, output_dir):
    """Copy anom file to output directory."""
    shutil.copy(anom_filename, output_dir)

def _get_metric_long_name(anom_filename):
    """Get long_name of anomaly variable."""
    basename = os.path.basename(anom_filename)
    varname = basename.split("_")[1]
    metric = basename.split("_")[2]
    if metric == "SANOM":
        long_name = _LONG_NAMES_SANOM[varname]
    else:
        long_name = _LONG_NAMES_ANOM[varname]
    return long_name

def _get_metric_units(anom_filename):
    """Get units of anomaly variable."""
    basename = os.path.basename(anom_filename)
    metric = basename.split("_")[2]
    if metric == "SANOM":
        metric_units = "1" # Standardize anomalies are dimensionless
    else:
        varname = basename.split("_")[1]
        metric_units = _UNITS_ANOM[varname]
    return metric_units

def _get_output_filename(anom_filename, output_dir):
    """Construct path to output file."""
    basename = os.path.basename(anom_filename)
    output_filename = f"{output_dir}/{basename}"
    return output_filename

def _driver():
    """Main driver."""

    # Read command line
    anom_filename, output_dir = _read_cmd_args()

    # Copy the file to the output directory, so we can edit without affecting
    # the original.
    _copy_anom_file(anom_filename, output_dir)

    # Find units of anomaly or standardized anomaly metric
    metric_units = _get_metric_units(anom_filename)
    metric_long_name = _get_metric_long_name(anom_filename)

    # Edit the attributes
    output_filename = _get_output_filename(anom_filename, output_dir)
    nc4 = Dataset(output_filename, 'r+', format='NETCDF4')
    nc4.Conventions = 'CF-1.8'
    nc4["latitude"].long_name = 'latitude'
    nc4["latitude"].standard_name = 'latitude'
    nc4["latitude"].units = 'degree_north'
    nc4["latitude"].axis = 'Y'

    nc4["longitude"].long_name = 'longitude'
    nc4["longitude"].standard_name = 'longitude'
    nc4["longitude"].units = 'degree_east'
    nc4["longitude"].axis = 'X'

    nc4["ens"].long_name = 'Ensemble members'
    nc4["ens"].axis = 'E'
    nc4["ens"].units = '1'

    nc4["time"].long_name = 'Forecast month'
    nc4["time"].units = 'months'

    nc4["anom"].long_name = metric_long_name
    nc4["anom"].units = metric_units
    nc4.close()

# Invoke the main driver
if __name__ == "__main__":
    _driver()
