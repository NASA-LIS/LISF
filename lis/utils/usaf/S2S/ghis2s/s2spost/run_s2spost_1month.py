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
# SCRIPT: run_s2spost_cf_1month.py
#
# PURPOSE: Automates generation of daily and monthly CF files from one month
# of LIS output.
#
# REVISION HISTORY:
# 23 Sep 2021: Eric Kemp (SSAI), first version.
# 27 Oct 2021: Eric Kemp/SSAI, address pylint string objections.
# 29 Oct 2021: Eric Kemp/SSAI, add config file.
# 02 Jun 2023: K. Arsenault, updated the s2spost filenaming conventions
#
#------------------------------------------------------------------------------
"""


# Standard modules
import datetime
from dateutil.relativedelta import relativedelta
import os
import pathlib
import subprocess
import sys
import shutil
import glob
import yaml
# pylint: disable=f-string-without-interpolation
# Local functions
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]}"
    txt += " configfile topdatadir fcstyyyymm YYYYMM model_forcing"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]   configfile is path to LDT parameter file")
    print("[INFO]   topdatadir is top-level directory for LIS data")
    print("[INFO]   fcstyyyymm is the initial forecast year and month [YYYYMM]")
    print("[INFO]   YYYYMM is the lead month to process")
    print("[INFO]   model_forcing is ID for atmospheric forcing for LIS")

def read_cmd_args(argv, arg_max):
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(argv) != arg_max:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Get path to config file
    configfile = argv[1]
    if not os.path.exists(configfile):
        print(f"[ERR] Config file {configfile} does not exist!")
        sys.exit(1)

    # Get top directory of LIS data
    topdatadir = argv[2]
    if not os.path.exists(topdatadir):
        print(f"[ERR] LIS data directory {topdatadir} does not exist!")
        sys.exit(1)

    # Get initial year and month
    fcstyyyymm = argv[3]
    if len(fcstyyyymm) != 6:
        print("[ERR] Invalid length of initial forecast yyyymm, must be 6 characters!")
        sys.exit(1)
    fcstyear = int(fcstyyyymm[0:4])
    fcstmonth = int(fcstyyyymm[4:6])
    try:
        fcstdate = datetime.datetime(fcstyear, fcstmonth, day=1)
    except ValueError:
        print("[ERR] Invalid initial forecast yyyymm passed to script!")
        sys.exit(1)

    # Get valid lead year and month
    yyyymm = argv[4]
    day=1
    if len(yyyymm) != 6:
        day = int(yyyymm[6:8])
    year = int(yyyymm[0:4])
    month = int(yyyymm[4:6])
    try:
        startdate = datetime.datetime(year, month, day=day)
    except ValueError:
        print("[ERR] Invalid lead YYYYMM passed to script!")
        sys.exit(1)

    # Get model forcing ID
    model_forcing = argv[5]

    return configfile, topdatadir, fcstdate, startdate, model_forcing

def _is_lis_output_missing(curdate, model_forcing):
    """Checks for missing LIS output files."""
    for model in ["SURFACEMODEL", "ROUTING"]:
        filename = f"lis_fcst"
        filename += f"/{model_forcing}"
        filename += f"/{model}"
        filename += f"/{curdate.year:04d}{curdate.month:02d}"
        filename += "/LIS_HIST_"
        filename += f"{curdate.year:04d}{curdate.month:02d}{curdate.day:02d}"
        filename += "0000.d01.nc"
        if not os.path.exists(filename):
            return True
    return False

def _loop_daily(config, configfile, topdatadir, fcstdate, startdate, model_forcing):
    """Automate daily processing for given month."""

    delta = datetime.timedelta(days=1)
    scriptdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/S2S/ghis2s/s2spost/'

    # The very first day may be missing. Gracefully handle this
    firstdate = startdate
    if _is_lis_output_missing(firstdate, model_forcing):
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
        cmd = f"python {scriptdir}/daily_s2spost_nc.py {configfile}"
        for model in ["SURFACEMODEL", "ROUTING"]:
            cmd += f" lis_fcst/{model_forcing}/{model}/"
            cmd += f"{curdate.year:04d}{curdate.month:02d}"
            cmd += "/LIS_HIST_"
            cmd += f"{curdate.year:04d}{curdate.month:02d}{curdate.day:02d}"
            cmd += "0000.d01.nc"

        cmd += f" {topdatadir}"

        cmd += f" {fcstdate.year:04d}{fcstdate.month:02d}{fcstdate.day:02d}"

        cmd += f" {curdate.year:04d}{curdate.month:02d}{curdate.day:02d}00"

        cmd += f" {model_forcing}"

        print(cmd)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem running CF conversion!")
            sys.exit(1)

        curdate += delta

def _proc_time_period(config, configfile, topdatadir, fcstdate, startdate, model_forcing, period):
    """Create the monthly CF file."""

    scriptdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/S2S/ghis2s/s2spost/'

    firstdate = startdate
    if period == "MONTHLY":
        # The very first day may be missing.  Gracefully handle this.    
        if _is_lis_output_missing(firstdate, model_forcing):
            firstdate += datetime.timedelta(days=1)

        if startdate.month == 12:
            enddate = datetime.datetime(year=(startdate.year + 1),
                                        month=1,
                                        day=1)
        else:
            enddate = datetime.datetime(year=startdate.year,
                                        month=(startdate.month + 1),
                                        day=1)
    if period == "WEEKLY":
        enddate = firstdate
        enddate += relativedelta(days=6)

    cmd = f"python {scriptdir}/temporal_aggregate.py {configfile} "

    cmd += f" {topdatadir}"
    cmd += f" {topdatadir}" 
    cmd += f" {fcstdate.year:04d}{fcstdate.month:02d}{fcstdate.day:02d}"
    cmd += f" {firstdate.year:04d}{firstdate.month:02d}{firstdate.day:02d}"
    cmd += f" {enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
    cmd += f" {model_forcing}"

    print(cmd)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print("[ERR] Problem creating monthly file!")
        sys.exit(1)

def _driver():
    """Main driver"""

    configfile, topdatadir, fcstdate, startdate, model_forcing = read_cmd_args(sys.argv, len(sys.argv))
    # load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if len(sys.argv) == 6:
        _loop_daily(config, configfile, topdatadir, fcstdate, startdate, model_forcing)
    else:
        _proc_time_period(config, configfile, topdatadir, fcstdate, startdate, model_forcing, sys.argv[6])

# Invoke driver
if __name__ == "__main__":
    _driver()
