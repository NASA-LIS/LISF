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
# SCRIPT: run_s2spost_cf_1month.py
#
# PURPOSE: Automates generation of daily and monthly CF files from one month
# of LIS output.
#
# REVISION HISTORY:
# 23 Sep 2021: Eric Kemp (SSAI), first version.
# 27 Oct 2021: Eric Kemp/SSAI, address pylint string objections.
# 29 Oct 2021: Eric Kemp/SSAI, add config file.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
import datetime
import os
import pathlib
import subprocess
import sys

# Local functions
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]}"
    txt += " configfile topdatadir YYYYMM model_forcing"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]   configfile is path to LDT parameter file")
    print("[INFO]   topdatadir is top-level directory for LIS data")
    print("[INFO]   YYYYMM is month to process")
    print("[INFO]   model_forcing is ID for atmospheric forcing for LIS")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 5:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Get path to config file
    configfile = sys.argv[1]
    if not os.path.exists(configfile):
        print(f"[ERR] Config file {configfile} does not exist!")
        sys.exit(1)

    # Get top directory of LIS data
    topdatadir = sys.argv[2]
    if not os.path.exists(topdatadir):
        print(f"[ERR] LIS data directory {topdatadir} does not exist!")
        sys.exit(1)

    # Get valid year and month
    yyyymm = sys.argv[3]
    if len(yyyymm) != 6:
        print("[ERR] Invalid length of YYYYMM, must be 6 characters!")
        sys.exit(1)
    year = int(yyyymm[0:4])
    month = int(yyyymm[4:6])
    try:
        startdate = datetime.datetime(year, month, day=1)
    except ValueError:
        print("[ERR] Invalid YYYYMM passed to script!")
        sys.exit(1)

    # Get model forcing ID
    model_forcing = sys.argv[4]

    return configfile, topdatadir, startdate, model_forcing

def _read_config(configfile):
    """Read from s2spost config file."""
    config = configparser.ConfigParser()
    config.read(configfile)
    return config

def _is_lis_output_missing(topdatadir, curdate, model_forcing):
    """Checks for missing LIS output files."""
    for model in ["SURFACEMODEL", "ROUTING"]:
        filename = f"{topdatadir}"
        filename += f"/{model_forcing}"
        filename += f"/{model}"
        filename += f"/{curdate.year:04d}{curdate.month:02d}"
        filename += "/LIS_HIST_"
        filename += f"{curdate.year:04d}{curdate.month:02d}{curdate.day:02d}"
        filename += "0000.d01.nc"
        if not os.path.exists(filename):
            return True
    return False

def _loop_daily(config, configfile, topdatadir, startdate, model_forcing):
    """Automate daily processing for given month."""

    delta = datetime.timedelta(days=1)

    scriptdir = config["s2spost"]["script_dir"]

    # The very first day may be missing. Gracefully handle this
    firstdate = startdate
    if _is_lis_output_missing(topdatadir, firstdate, model_forcing):
        firstdate += delta

    if startdate.month == 12:
        enddate = datetime.datetime(year=(startdate.year + 1),
                                    month=1,
                                    day=1)
    else:
        enddate = datetime.datetime(year=startdate.year,
                                    month=(startdate.month + 1),
                                    day=1)

    curdate = firstdate
    while curdate <= enddate:
        cmd = f"{scriptdir}/daily_s2spost_nc.py {configfile}"
        for model in ["SURFACEMODEL", "ROUTING"]:
            cmd += f" {topdatadir}/{model_forcing}/{model}/"
            cmd += f"{curdate.year:04d}{curdate.month:02d}"
            cmd += "/LIS_HIST_"
            cmd += f"{curdate.year:04d}{curdate.month:02d}{curdate.day:02d}"
            cmd += "0000.d01.nc"

        cmd += f" {topdatadir}/cf_{model_forcing}_"
        cmd += f"{startdate.year:04d}{startdate.month:02d}"

        cmd += f" {curdate.year:04d}{curdate.month:02d}{curdate.day:02d}00"

        cmd += f" {model_forcing}"

        print(cmd)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem running CF conversion!")
            sys.exit(1)

        curdate += delta

def _proc_month(config, configfile, topdatadir, startdate, model_forcing):
    """Create the monthly CF file."""

    scriptdir = config["s2spost"]["script_dir"]

    # The very first day may be missing.  Gracefully handle this.
    firstdate = startdate
    if _is_lis_output_missing(topdatadir, firstdate, model_forcing):
        firstdate += datetime.timedelta(days=1)

    if startdate.month == 12:
        enddate = datetime.datetime(year=(startdate.year + 1),
                                    month=1,
                                    day=1)
    else:
        enddate = datetime.datetime(year=startdate.year,
                                    month=(startdate.month + 1),
                                    day=1)
    cmd = f"{scriptdir}/monthly_s2spost_nc.py {configfile} "
    workdir =  f"{topdatadir}/cf_{model_forcing}_"
    workdir += f"{startdate.year:04d}{startdate.month:02d}"

    cmd += f" {workdir}"
    cmd += f" {workdir}" # Use same directory
    cmd += f" {firstdate.year:04d}{firstdate.month:02d}{firstdate.day:02d}"
    cmd += f" {enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
    cmd += f" {model_forcing}"

    print(cmd)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print("[ERR] Problem creating monthly file!")
        sys.exit(1)

def _create_done_file(topdatadir, startdate, model_forcing):
    """Create a 'done' file indicating batch job has finished."""
    path = f"{topdatadir}/cf_{model_forcing}_"
    path += f"{startdate.year:04d}{startdate.month:02d}/done"
    pathlib.Path(path).touch()

def _driver():
    """Main driver"""
    configfile, topdatadir, startdate, model_forcing = _read_cmd_args()
    config = _read_config(configfile)
    _loop_daily(config, configfile, topdatadir, startdate, model_forcing)
    _proc_month(config, configfile, topdatadir, startdate, model_forcing)
    _create_done_file(topdatadir, startdate, model_forcing)

# Invoke driver
if __name__ == "__main__":
    _driver()
