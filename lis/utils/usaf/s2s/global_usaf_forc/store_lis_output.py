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
# SCRIPT: store_lis_output.py
#
# PURPOSE: Moves output from LIS USAF forcing run to appropriate new locations,
# purging as needed.
#
# REQUIREMENTS as of 21 Nov 2022:
# * Python 3.9 or higher.
#
# REVISION HISTORY:
# 21 Sep 2021: Eric Kemp (SSAI), first version.
# 02 Nov 2021: Eric Kemp/SSAI, tweaks to appease pylint.
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import shutil
import sys

def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} tmp_output_dir restart_dir grib_dir" + \
        " YYYYMMDD"
    print(txt)
    print("[INFO] where: ")
    print("[INFO]  tmp_output_dir: Path to output files from recent LIS run.")
    print("[INFO]  restart_dir: Path to permanently store restart files.")
    print("[INFO]  grib_dir: Path to permanently store GRIB files for S2S.")
    print("[INFO]  YYYYMMDD: year/month/day of start of recent LIS run.")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 5:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if paths exist.
    tmp_output_dir = sys.argv[1]
    if not os.path.exists(tmp_output_dir):
        print(f"[ERR] {tmp_output_dir} does not exist!")
        sys.exit(1)

    restart_dir = sys.argv[2]
    if not os.path.exists(restart_dir):
        print(f"[ERR] {restart_dir} does not exist!")
        sys.exit(1)

    grib_dir = sys.argv[3]
    if not os.path.exists(grib_dir):
        os.makedirs(grib_dir)

    # Get start date of recent LIS run
    yyyymmdd = sys.argv[4]
    if len(yyyymmdd) != 8:
        print("[ERR] Invalid length for YYYYMMDD, must be 8 characters!")
        sys.exit(1)
    year = int(yyyymmdd[0:4])
    month = int(yyyymmdd[4:6])
    day = int(yyyymmdd[6:8])
    try:
        startdate = datetime.date(year, month, day)
    except ValueError:
        print("[ERR] Invalid YYYYMMDD passed to script!")
        sys.exit(1)

    return tmp_output_dir, restart_dir, grib_dir, startdate

def _copy_restart_file(tmp_output_dir, restart_dir, startdate):
    """Copy last restart file to more permanent storage."""
    enddate = startdate + datetime.timedelta(days=1)

    restart_file = f"{tmp_output_dir}/LIS_RST_NOAH39_" + \
        f"{enddate.year:04d}{enddate.month:02d}{enddate.day:02d}" + \
        "0000.d01.nc"
    if not os.path.exists(restart_file):
        print(f"[ERR] Restart file {restart_file} does not exist!")
        sys.exit(1)

    newfile = shutil.copy(restart_file, restart_dir)
    if not os.path.exists(newfile):
        print(f"[ERR] Problem copying {restart_file} to {restart_dir}")
        sys.exit(1)

def _copy_grib_files(tmp_output_dir, grib_dir, startdate):
    """Copy grib files to more permanent storage. EXCEPTION: Skip first
    file since that has zero precipitation."""

    startdt = datetime.datetime(year=startdate.year,
                                month=startdate.month,
                                day=startdate.day,
                                hour=3)
    enddate = startdate + datetime.timedelta(days=1)
    enddt = datetime.datetime(year=enddate.year,
                              month=enddate.month,
                              day=enddate.day,
                              hour=0)
    delta = datetime.timedelta(seconds=10800)

    curdt = startdt
    while curdt <= enddt:
        grib_file = f"{tmp_output_dir}/"
        grib_file += "PS.AFWA_SC.U_DI.C_DC.ANLYS_GP.LIS_GR.C0P09DEG"
        grib_file += "_AR.GLOBAL_PA.03-HR-SUM"
        grib_file += f"_DD.{curdt.year:04d}{curdt.month:02d}{curdt.day:02d}"
        grib_file += f"_DT.{curdt.hour:02d}00"
        grib_file += "_DF.GR1"
        if not os.path.exists(grib_file):
            print(f"[ERR] {grib_file} does not exist!")
            sys.exit(1)

        newfile = shutil.copy(grib_file, grib_dir)
        if not os.path.exists(newfile):
            print(f"[ERR] Problem copying {grib_file} to {grib_dir}")
            sys.exit(1)

        curdt += delta

def _purge_tmp_output_dir(tmp_output_dir):
    """Purges temporary output directory from recent LIS run."""
    for filename in os.listdir(tmp_output_dir):
        path = os.path.join(tmp_output_dir, filename)
        try:
            os.remove(path)
        except IsADirectoryError:
            continue

def _driver():
    """Main driver."""
    tmp_output_dir, restart_dir, grib_dir, startdate = _read_cmd_args()
    _copy_restart_file(tmp_output_dir, restart_dir, startdate)
    _copy_grib_files(tmp_output_dir, grib_dir, startdate)
    _purge_tmp_output_dir(tmp_output_dir)

# Invoke driver
if __name__ == "__main__":
    _driver()
