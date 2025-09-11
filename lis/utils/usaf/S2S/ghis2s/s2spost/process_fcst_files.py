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
# SCRIPT: process_fcst_files.py
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
# Local modules
from ghis2s.shared.logging_utils import TaskLogger
from merge_lisf_files import create_final_filename, merge_files_xarray
from temporal_aggregate import agg_driver

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

    task_name = os.environ.get('SCRIPT_NAME')
    logger = TaskLogger(task_name,
                        os.getcwd(),
                        f's2spost/process_fcst_files.py processing daily {model_forcing} for month {startdate.year:04d}{startdate.month:02d}')
    
    delta = datetime.timedelta(days=1)
    scriptdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/S2S/ghis2s/s2spost/'
    ldtfile = config['SETUP']['supplementarydir'] + '/lis_darun/' + config["SETUP"]["ldtinputfile"]

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
        subtask = f'{model_forcing} {startdate.year:04d}{startdate.month:02d}'
        model = "SURFACEMODEL"
        noahmp_file = f"lis_fcst/{model_forcing}/{model}/{curdate.year:04d}{curdate.month:02d}/LIS_HIST_{curdate.year:04d}{curdate.month:02d}{curdate.day:02d}0000.d01.nc"
        model = "ROUTING"
        hymap_file = f"lis_fcst/{model_forcing}/{model}/{curdate.year:04d}{curdate.month:02d}/LIS_HIST_{curdate.year:04d}{curdate.month:02d}{curdate.day:02d}0000.d01.nc"
        print()
        if not os.path.exists(noahmp_file):
            logger.error(f"Missing: {noahmp_file}", subtask=f'{model_forcing} {startdate.year:04d}{startdate.month:02d}')
            sys.exit()
        if not os.path.exists(hymap_file):
            logger.error(f"Missing: {hymap_file}", subtask=f'{model_forcing} {startdate.year:04d}{startdate.month:02d}')
            sys.exit()

        merge_file = create_final_filename(topdatadir, fcstdate, curdate, model_forcing.upper(), config["EXP"]["DOMAIN"])
        logger.info(f's2spost/merge_lisf_files.py processing {curdate.year:04d}{curdate.month:02d}{curdate.day:02d}',
                    subtask=subtask)
        try:
            merge_files_xarray(ldtfile, noahmp_file, hymap_file, merge_file, fcstdate, logger, subtask)
            logger.info(f'Merged file: {merge_file}', subtask=f'{model_forcing} {startdate.year:04d}{startdate.month:02d}')
        except Exception as e:
            logger.error(f"Failed processing {curdate.year:04d}{curdate.month:02d}{curdate.day:02d}: {str(e)}", 
                        subtask=f'{model_forcing} {startdate.year:04d}{startdate.month:02d}')
            sys.exit(1)

        curdate += delta

def _proc_time_period(config, configfile, topdatadir, fcstdate, startdate, model_forcing, period):
    """Create the monthly CF file."""

    task_name = os.environ.get('SCRIPT_NAME')
    logger = TaskLogger(task_name,
                        os.getcwd(),
                        f's2spost/process_fcst_files.py processing {model_forcing} {period} {startdate.year:04d}{startdate.month:02d}')

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

    argv = []
    argv.append(configfile)
    argv.append(topdatadir)
    argv.append(topdatadir)
    argv.append(f"{fcstdate.year:04d}{fcstdate.month:02d}{fcstdate.day:02d}")
    argv.append(f"{firstdate.year:04d}{firstdate.month:02d}{firstdate.day:02d}")
    argv.append(f"{enddate.year:04d}{enddate.month:02d}{enddate.day:02d}")
    argv.append(model_forcing)
    subtask = f'{model_forcing} {startdate.year:04d}{startdate.month:02d}'
    logger.info(f's2spost/temporal_aggregate.py processing {firstdate.year:04d}{firstdate.month:02d}-{enddate.year:04d}{enddate.month:02d}', subtask=subtask)
    agg_driver(argv, logger, subtask)

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
