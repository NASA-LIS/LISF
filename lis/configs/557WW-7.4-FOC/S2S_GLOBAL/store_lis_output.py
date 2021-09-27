#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: store_lis_usaf_atmos_output.py
#
# PURPOSE: Moves output from LIS USAF forcing run to appropriate new locations,
# purging as needed.
#
# REQUIREMENTS as of 21 Sep 2021:
# * Python 3.8 or higher.
#
# REVISION HISTORY:
# 21 Sep 2021: Eric Kemp (SSAI), first version.
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import shutil
import sys

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s tmp_output_dir restart_dir grib_dir YYYYMMDD" \
        %(sys.argv[0])
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
        print("[ERR] %s does not exist!" %(tmp_output_dir))
        sys.exit(1)

    restart_dir = sys.argv[2]
    if not os.path.exists(restart_dir):
        print("[ERR] %s does not exist!" %(restart_dir))
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

    restart_file = "%s/LIS_RST_NOAH39_%4.4d%2.2d%2.2d0000.d01.nc" \
        %(tmp_output_dir, enddate.year, enddate.month, enddate.day)
    if not os.path.exists(restart_file):
        print("[ERR] Restart file %s does not exist!" %(restart_file))
        sys.exit(1)

    newfile = shutil.copy(restart_file, restart_dir)
    if not os.path.exists(newfile):
        print("[ERR] Problem copying %s to %s" %(restart_file, restart_dir))
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
        grib_file = "%s/" %(tmp_output_dir)
        grib_file += "PS.AFWA_SC.U_DI.C_DC.ANLYS_GP.LIS_GR.C0P09DEG"
        grib_file += "_AR.GLOBAL_PA.03-HR-SUM"
        grib_file += "_DD.%4.4d%2.2d%2.2d" %(curdt.year,
                                             curdt.month,
                                             curdt.day)
        grib_file += "_DT.%2.2d00" %(curdt.hour)
        grib_file += "_DF.GR1"
        if not os.path.exists(grib_file):
            print("[ERR] %s does not exist!" %(grib_file))
            sys.exit(1)

        newfile = shutil.copy(grib_file, grib_dir)
        if not os.path.exists(newfile):
            print("[ERR] Problem copying %s to %s" %(grib_file, grib_dir))
            sys.exit(1)

        curdt += delta

def _purge_tmp_output_dir(tmp_output_dir):
    """Purges temporary output directory from recent LIS run."""
    for filename in os.listdir(tmp_output_dir):
        os.remove(os.path.join(tmp_output_dir, filename))

def _driver():
    """Main driver."""
    tmp_output_dir, restart_dir, grib_dir, startdate = _read_cmd_args()
    _copy_restart_file(tmp_output_dir, restart_dir, startdate)
    _copy_grib_files(tmp_output_dir, grib_dir, startdate)
    _purge_tmp_output_dir(tmp_output_dir)

# Invoke driver
if __name__ == "__main__":
    _driver()
