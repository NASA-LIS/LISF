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
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import subprocess
import sys

# Path to NCO binaries.
_NCO_DIR = "/usr/local/other/nco/5.0.1/bin" # On Discover

# Lists of variables and metrics to process
_VAR_LIST = ["RootZone-SM", "Streamflow", "Surface-SM"]
_METRIC_LIST = ["ANOM", "SANOM"]

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s input_dir output_dir" %(sys.argv[0])
    txt += " start_yyyymmdd end_yyyymmdd"
    txt += " model_forcing"
    print(txt)
    print("[INFO] where:")
    print("[INFO] input_dir: directory with S2S metric files in CF convention")
    print("[INFO] output_dir: directory for output file")
    print("[INFO] start_yyyymmdd: Starting date/time of metrics files")
    print("[INFO] end_yyyymmdd: Starting date/time of metrics files")
    print("[INFO] model_forcing; ID for atmospheric forcing for LIS")

def _check_nco_binaries():
    """Check to see if necessary NCO binaries are available."""
    nco_bins = ["ncks", "ncrename"]
    for nco_bin in nco_bins:
        path = "%s/%s" %(_NCO_DIR, nco_bin)
        if not os.path.exists(path):
            print("[ERR] Cannot find %s for converting LIS netCDF4 data!" \
                  %(path))
            print("[ERR] Make sure NCO package is installed on the system!")
            print("[ERR] And update _NCO_DIR in this script if necessary!")
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

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 6:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if input directory exists.
    input_dir = sys.argv[1]
    if not os.path.exists(input_dir):
        print("[ERR] Directory %s does not exist!")
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

    return input_dir, output_dir, startdate, enddate, model_forcing

def _check_filename_size(name):
    """Make sure filename does not exceed 128 characters, per Air Force
    requirement."""
    if len(os.path.basename(name)) > 128:
        print("[ERR] Output file name is too long!")
        print("[ERR] %s exceeds 128 characters!" %(os.path.basename(name)))
        sys.exit(1)

def _create_var_metric_filename(input_dir, model_forcing, var, metric,
                                startdate):
    """Create path to S2S metric file."""
    name = "%s/" %(input_dir)
    name += "%s_" %(model_forcing)
    name += "%s_" %(var)
    name += "%s_" %(metric)
    name += "init_monthly_%2.2d_%4.4d" %(startdate.month,
                                         startdate.year)
    name += ".nc"
    return name

def _create_merged_metric_filename(output_dir, startdate, enddate,
                                   model_forcing):
    """Create path to merged S2S metric netCDF file."""
    name = "%s" %(output_dir)
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += "_GP.LIS-S2S-%s-ANOM" %(model_forcing.upper())
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S-ANOM"
    name += "_DP.%4.4d%2.2d%2.2d-%4.4d%2.2d%2.2d" \
        %(startdate.year, startdate.month, startdate.day,
          enddate.year, enddate.month, enddate.day)
    name += "_TP.0000-0000"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

def _merge_files(input_dir, model_forcing, startdate, mergefile):
    """Merge individual variable metrics into single file."""

    # Copy first variable/metric file
    first_var = _VAR_LIST[0]
    first_metric = _METRIC_LIST[0]
    metricfile =  _create_var_metric_filename(input_dir, model_forcing,
                                              first_var, first_metric,
                                              startdate)
    if not os.path.exists(metricfile):
        print("[ERR] %s does not exist!" %(metricfile))
        sys.exit(1)

    cmd = "%s/ncks" %(_NCO_DIR)
    cmd += " %s -6 %s" %(metricfile, mergefile)
    _run_cmd(cmd, "[ERR] Problem with ncks!")

    # Rename ANOM array
    cmd = "%s/ncrename -O" %(_NCO_DIR)
    cmd += " -v anom,%s_%s" %(first_var.replace("-","_"), first_metric)
    cmd += " %s" %(mergefile)
    _run_cmd(cmd, "[ERR] Problem with ncrename!")

    # Loop through remaining var metric files, copy *only* the anom variable,
    # and then rename the anom variable.
    for var in _VAR_LIST:
        for metric in _METRIC_LIST:

            if (var, metric) == (first_var, first_metric):
                continue

            metricfile =  _create_var_metric_filename(input_dir, model_forcing,
                                                      var, metric,
                                                      startdate)
            if not os.path.exists(metricfile):
                print("[ERR] %s does not exist!" %(metricfile))
                sys.exit(1)

            cmd = "%s/ncks -A -C" %(_NCO_DIR)
            cmd += " -v anom"
            cmd += " %s -6 %s" %(metricfile, mergefile)
            _run_cmd(cmd, "[ERR] Problem with ncks!")

            cmd = "%s/ncrename -O" %(_NCO_DIR)
            cmd += " -v anom,%s_%s" %(var.replace("-","_"), metric)
            cmd += " %s" %(mergefile)
            _run_cmd(cmd, "[ERR] Problem with ncrename!")

def _copy_to_final_file(mergefile, final_file):
    """Copy to new file, with netCDF4 compression."""
    cmd = "%s/ncks" %(_NCO_DIR)
    cmd += " %s -7 -L 1 %s" %(mergefile, final_file)
    _run_cmd(cmd, "[ERR] Problem with ncks!")

def _driver():
    """Main driver"""
    _check_nco_binaries()
    input_dir, output_dir, startdate, enddate, model_forcing = _read_cmd_args()
    output_filename = _create_merged_metric_filename(output_dir,
                                                     startdate, enddate,
                                                     model_forcing)
    tmp_output_filename = "%s/tmp.nc" %(output_dir)
    _merge_files(input_dir, model_forcing, startdate, tmp_output_filename)
    _copy_to_final_file(tmp_output_filename, output_filename)
    os.unlink(tmp_output_filename)

# Invoke driver
if __name__ == "__main__":
    _driver()
