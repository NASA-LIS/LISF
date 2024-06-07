#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
#------------------------------------------------------------------------------
#
# SCRIPT: merge_s2s_anom_cf.py
#
# PURPOSE:  Merge CF-convention S2S anomaly files into single file for sharing.
#
# REQUIREMENTS as of 28 May 2023:
# * Python 3.9 or higher
#
# REFERENCES:
#
# REVISION HISTORY:
# 27 Sep 2021: Eric Kemp (SSAI), first version.
# 30 Oct 2021: Eric Kemp/SSAI, added support for s2smetric config file.
# 02 Jun 2023: K. Arsenault + S. Mahanama, updated 557 WW file names.
#
#------------------------------------------------------------------------------
"""


# Standard modules
import os
import sys
import datetime
import yaml
import xarray as xr

# Local constants
_METRIC_LIST = ["ANOM", "SANOM"]

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} input_dir output_dir"
    txt += " start_yyyymmdd end_yyyymmdd"
    txt += " model_forcing configfile"
    print(txt)
    print("[INFO] where:")
    print("[INFO] input_dir: directory with S2S metric files in CF convention")
    print("[INFO] output_dir: directory for output file")
    print("[INFO] start_yyyymmdd: Starting date/time of metrics files")
    print("[INFO] end_yyyymmdd: Starting date/time of metrics files")
    print("[INFO] model_forcing; ID for atmospheric forcing for LIS")
    print("[INFO] configfile: Path to s2smetric config file")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 7:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if input directory exists.
    input_dir = sys.argv[1]
    if not os.path.exists(input_dir):
        print(f"[ERR] Directory {input_dir} does not exist!")
        sys.exit(1)

    # Create output directory if it doesn't exist.
    output_dir = sys.argv[2]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get valid starting and ending dates of data.
    start_yyyymmdd = sys.argv[3]
    startdate = _proc_date(start_yyyymmdd)
    end_yyyymmdd = sys.argv[4]
    enddate = _proc_date(end_yyyymmdd)
    if startdate > enddate:
        print("[ERR] Start date is after end date!")
        sys.exit(1)

    # Get ID for model forcing for LIS
    model_forcing = sys.argv[5]

    configfile = sys.argv[6]
    if not os.path.exists(configfile):
        print(f"[ERR] Cannot file config file {configfile}!")
        sys.exit(1)

    return input_dir, output_dir, startdate, enddate, model_forcing, configfile

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

def _create_var_metric_filename(input_dir, model_forcing, var, metric,
                                startdate):
    """Create path to S2S metric file."""
    name = f"{input_dir}/"
    name += f"{model_forcing}_"
    name += f"{var}_"
    name += f"{metric}_"
    name += f"init_monthly_{startdate.month:02d}_{startdate.year:04d}"
    name += ".nc"
    return name

def _create_merged_metric_filename(output_dir, startdate, enddate,
                                   model_forcing, domain):
    """Create path to merged S2S metric netCDF file."""
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing.upper()}"
    name += "_GR.C0P25DEG"
    if domain == 'AFRICOM':
        name += "_AR.AFRICA"
    if domain == 'GLOBAL':
        name += "_AR.GLOBAL"
    name += "_PA.S2SMETRICS"
    name += f"_DD.{startdate.year:04d}{startdate.month:02d}01"
    name += f"_FP.{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
    name += f"-{enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

def _merge_files(config, input_dir, model_forcing, startdate, mergefile):
    """Merge individual variable metrics into single file."""

    # Copy first variable/metric file
    var_list = config["POST"]["metric_vars"]
    first_var = var_list[0]
    first_metric = _METRIC_LIST[0]
    metricfile = _create_var_metric_filename(input_dir, model_forcing,
                                             first_var, first_metric,
                                             startdate)
    if not os.path.exists(metricfile):
        print(f"[ERR] {metricfile} does not exist!")
        sys.exit(1)

    ds_out = xr.open_dataset(metricfile)
    ds_out = ds_out.rename ({"anom": first_var.replace('-','_') + '_' + first_metric})

    # Loop through remaining var metric files, copy *only* the anom variable,
    # and then rename the anom variable.
    for var in var_list:
        for metric in _METRIC_LIST:

            if (var, metric) == (first_var, first_metric):
                continue

            metricfile = _create_var_metric_filename(input_dir, model_forcing,
                                                     var, metric,
                                                     startdate)

            if not os.path.exists(metricfile):
                print(f"[ERR] {metricfile} does not exist!")
                sys.exit(1)
            _ds = xr.open_dataset(metricfile)
            ds_out = xr.merge([ds_out, _ds.rename ({"anom": var.replace('-','_') + '_' + metric})])

    comp = dict(zlib=True, complevel=6)
    encoding = {var: comp for var in ds_out.data_vars}
    ds_out.to_netcdf(mergefile, encoding=encoding)

def _driver():
    """Main driver"""
    input_dir, output_dir, startdate, enddate, model_forcing, configfile \
        = _read_cmd_args()
    # load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)
    output_filename = _create_merged_metric_filename(output_dir,
                                                     startdate, enddate,
                                                     model_forcing, config["EXP"]["DOMAIN"])
    _merge_files(config, input_dir, model_forcing, startdate, \
                 output_filename)

# Invoke driver
if __name__ == "__main__":
    _driver()
