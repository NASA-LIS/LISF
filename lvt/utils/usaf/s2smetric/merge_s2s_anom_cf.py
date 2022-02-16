#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: merge_s2s_anom_cf.py
#
# PURPOSE:  Merge CF-convention S2S anomaly files into single file for sharing.
#
# REQUIREMENTS as of 27 Sep 2021:
# * Python 3.8 or higher
# * netCDF Operator (NCO) binaries version 5.0.1 or later.
#
# REFERENCES:
# https://nco.sourceforge.net for NCO utilities documentation and source code.
#
# REVISION HISTORY:
# 27 Sep 2021: Eric Kemp (SSAI), first version.
# 30 Oct 2021: Eric Kemp/SSAI, added support for s2smetric config file.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
import datetime
import os
import subprocess
import sys

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

def _check_nco_binaries(config):
    """Check to see if necessary NCO binaries are available."""
    ncodir = config["s2smetric"]["ncodir"]
    nco_bins = ["ncks", "ncrename"]
    for nco_bin in nco_bins:
        path = f"{ncodir}/{nco_bin}"
        if not os.path.exists(path):
            print(f"[ERR] Cannot find {path} for converting LIS netCDF4 data!")
            print("[ERR] Make sure NCO package is installed on the system!")
            print("[ERR] And update ncodir in s2smetric config if necessary!")
            sys.exit(1)

def _run_cmd(cmd, error_msg):
    """Handle running shell command and checking error."""
    print(cmd)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print(error_msg)
        sys.exit(1)

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
                                   model_forcing):
    """Create path to merged S2S metric netCDF file."""
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing.upper()}-ANOM"
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S-ANOM"
    name += f"_DP.{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
    name += f"-{enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
    name += "_TP.0000-0000"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

def _merge_files(config, input_dir, model_forcing, startdate, mergefile):
    """Merge individual variable metrics into single file."""

    # Copy first variable/metric file
    var_list = config["s2smetric"]["metric_vars"].split()
    first_var = var_list[0]
    first_metric = _METRIC_LIST[0]
    metricfile =  _create_var_metric_filename(input_dir, model_forcing,
                                              first_var, first_metric,
                                              startdate)
    if not os.path.exists(metricfile):
        print(f"[ERR] {metricfile} does not exist!")
        sys.exit(1)

    ncodir = config["s2smetric"]["ncodir"]
    cmd = f"{ncodir}/ncks"
    cmd += f" {metricfile} -6 {mergefile}"
    _run_cmd(cmd, "[ERR] Problem with ncks!")

    # Rename ANOM array
    cmd = f"{ncodir}/ncrename -O"
    cmd += f" -v anom,{first_var.replace('-','_')}_{first_metric}"
    cmd += f" {mergefile}"
    _run_cmd(cmd, "[ERR] Problem with ncrename!")

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

            cmd = f"{ncodir}/ncks -A -C"
            cmd += " -v anom"
            cmd += f" {metricfile} -6 {mergefile}"
            _run_cmd(cmd, "[ERR] Problem with ncks!")

            cmd = f"{ncodir}/ncrename -O"
            cmd += f" -v anom,{var.replace('-','_')}_{metric}"
            cmd += f" {mergefile}"
            _run_cmd(cmd, "[ERR] Problem with ncrename!")

def _copy_to_final_file(config, mergefile, final_file):
    """Copy to new file, with netCDF4 compression."""
    ncodir = config["s2smetric"]["ncodir"]
    cmd = f"{ncodir}/ncks"
    cmd += f" {mergefile} -7 -L 1 {final_file}"
    _run_cmd(cmd, "[ERR] Problem with ncks!")

def _driver():
    """Main driver"""
    input_dir, output_dir, startdate, enddate, model_forcing, configfile \
        = _read_cmd_args()
    config = configparser.ConfigParser()
    config.read(configfile)
    _check_nco_binaries(config)
    output_filename = _create_merged_metric_filename(output_dir,
                                                     startdate, enddate,
                                                     model_forcing)
    tmp_output_filename = f"{output_dir}/tmp.nc"
    _merge_files(config, input_dir, model_forcing, startdate, \
                 tmp_output_filename)
    _copy_to_final_file(config, tmp_output_filename, output_filename)
    os.unlink(tmp_output_filename)

# Invoke driver
if __name__ == "__main__":
    _driver()
