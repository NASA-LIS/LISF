#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: convert_s2s_anom_cf.py
#
# PURPOSE: Updates attributes of S2S monthly anomaly/standardized anomaly
# netCDF files to conform with CF-1.8 convention.
#
# REQUIREMENTS as of 26 Sep 2021:
# * Python 3.8 or higher
# * netCDF Operator (NCO) binaries version 5.0.1 or later.
#
# REFERENCES:
# https://cfconventions.org for specifications of NetCDF Climate and Forecast
#   (CF) Metadata Conventions.
# https://nco.sourceforge.net for NCO utilities documentation and source code.
# https://unidata.ucar.edu/software/udunits for documentation on UDUNITS2
#   library, which CF is generally consistent with for unit specifications.
#
# REVISION HISTORY:
# 26 Sep 2021: Eric Kemp (SSAI), first version.
# 30 Oct 2021: Eric Kemp/SSAI, now uses s2smetric config file.
#
#------------------------------------------------------------------------------
"""

# Standard modules

import os
import shutil
import subprocess
import sys

# Local constants
# Units for variable anomalies.  Standardized anomalies will be dimensionless.
_UNITS_ANOM = {
    "RootZone-SM" : "m3 m-3",
    "Streamflow" : "m3 s-1",
    "Surface-SM" : "m3 m-3",
    "TWS" : "mm",
    "Precip" : "kg m-2",
    "Air_T" : "K",
    "ET" : "kg m-2 s-1",
    "Air" : "K",
}
_LONG_NAMES_SANOM = {
    "RootZone-SM" : "Root zone soil moisture standardized anomaly",
    "Streamflow" : "Streamflow standardized anomaly",
    "Surface-SM" : "Surface soil moisture standardized anomaly",
    "TWS" : "Terrestrial water storage standardized anomaly",
    "Precip" : "Total precipitation amount standardized anomaly",
    "Air_T" : "Air temperature standardized anomaly",
    "ET" : "Total evapotranspiration standardized anomaly",
    "Air" : "Air temperature standardized anomaly",
}
_LONG_NAMES_ANOM = {
    "RootZone-SM" : "Root zone soil moisture anomaly",
    "Streamflow" : "Streamflow anomaly",
    "Surface-SM" : "Surface soil moisture anomaly",
    "TWS" : "Terrestrial water storage anomaly",
    "Precip" : "Total precipitation amount anomaly",
    "Air_T" : "Air temperature anomaly",
    "ET" : "Total evapotranspiration anomaly",
    "Air" : "Air temperature anomaly",
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

def _run_cmd(cmd, error_msg):
    """Handle running shell command and checking error."""
    #print(cmd)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print(error_msg)
        sys.exit(1)

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

def _update_global_attrs(output_filename):
    """Update global attributes of output filename."""
    cmd = f"ncatted"
    cmd += " -a Conventions,global,c,c,'CF-1.8'"
    cmd += f" {output_filename}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_latitude_attrs(output_filename):
    """Update attributes of latitude."""
    cmd = f"ncatted"
    cmd += " -a long_name,latitude,c,c,'latitude'"
    cmd += " -a standard_name,latitude,c,c,'latitude'"
    cmd += " -a units,latitude,c,c,'degree_north'"
    cmd += f" {output_filename}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_longitude_attrs(output_filename):
    """Update attributes of longitude."""
    cmd = f"ncatted"
    cmd += " -a long_name,longitude,c,c,'longitude'"
    cmd += " -a standard_name,longitude,c,c,'longitude'"
    cmd += " -a units,longitude,c,c,'degree_east'"
    cmd += f" {output_filename}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_ens_attrs(output_filename):
    """Update attributes of ens."""
    cmd = f"ncatted"
    cmd += " -a long_name,ens,c,c,'Ensemble members'"
    cmd += " -a units,ens,c,c,'1'"
    cmd += f" {output_filename}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_lead_attrs(output_filename):
    """Update attributes of lead."""
    cmd = f"ncatted"
    cmd += " -a long_name,lead,c,c,'Forecast month'"
    cmd += " -a units,lead,c,c,'months'"
    cmd += f" {output_filename}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_anom_attrs(output_filename, metric_long_name,
                       metric_units):
    """Update attributes of anom variable."""
    cmd = f"ncatted"
    cmd += f" -a long_name,anom,c,c,'{metric_long_name}'"
    cmd += f" -a units,anom,c,c,'{metric_units}'"
    cmd += f" {output_filename}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

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
    _update_global_attrs(output_filename)
    _update_latitude_attrs(output_filename)
    _update_longitude_attrs(output_filename)
    _update_ens_attrs(output_filename)
    _update_lead_attrs(output_filename)
    _update_anom_attrs(output_filename, metric_long_name, metric_units)

# Invoke the main driver
if __name__ == "__main__":
    _driver()
