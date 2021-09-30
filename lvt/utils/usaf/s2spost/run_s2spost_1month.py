#!/usr/bin/env python3
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
#
#------------------------------------------------------------------------------
"""

import datetime
import os
import subprocess
import sys

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s ldt_file topdatadir YYYYMM model_forcing" \
        %(sys.argv[0])
    print(txt)
    print("[INFO]  where:")
    print("[INFO]   ldt_file is path to LDT parameter file")
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

    # Assume all scripts are bundled together in the dirname of this script.
    scriptdir = os.path.dirname(sys.argv[0])

    # Get path to LDT parameter file
    ldtfile = sys.argv[1]
    if not os.path.exists(ldtfile):
        print("[ERR] LDT paramter file %s does not exist!" %(ldtfile))
        sys.exit(1)

    # Get top directory of LIS data
    topdatadir = sys.argv[2]
    if not os.path.exists(topdatadir):
        print("[ERR] LIS data directory %s does not exist!" %(topdatadir))
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

    return scriptdir, ldtfile, topdatadir, startdate, model_forcing

def _is_lis_output_missing(topdatadir, curdate, model_forcing):
    """Checks for missing LIS output files."""
    for model in ["SURFACEMODEL", "ROUTING"]:
        filename = "%s/%s/%s/%4.4d%2.2d" %(topdatadir,
                                           model_forcing,
                                           model,
                                           curdate.year,
                                           curdate.month)
        filename += "/LIS_HIST_%4.4d%2.2d%2.2d0000.d01.nc" %(curdate.year,
                                                             curdate.month,
                                                             curdate.day)
        if not os.path.exists(filename):
            return True
    return False

def _loop_daily(scriptdir, ldtfile, topdatadir, startdate, model_forcing):
    """Automate daily processing for given month."""

    delta = datetime.timedelta(days=1)

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
        cmd = "%s/daily_s2spost_nc.py" %(scriptdir)
        cmd += " %s" %(ldtfile)
        for model in ["SURFACEMODEL", "ROUTING"]:
            cmd += " %s/%s/%s/%4.4d%2.2d" %(topdatadir,
                                            model_forcing,
                                            model,
                                            curdate.year,
                                            curdate.month)
            cmd += "/LIS_HIST_%4.4d%2.2d%2.2d0000.d01.nc" %(curdate.year,
                                                            curdate.month,
                                                            curdate.day)
        cmd += " ./cf_%s_%4.4d%2.2d" %(model_forcing.upper(),
                                       startdate.year,
                                       startdate.month)
        cmd += " %4.4d%2.2d%2.2d00" %(curdate.year,
                                      curdate.month,
                                      curdate.day)
        cmd += " %s" %(model_forcing.upper())

        print(cmd)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem running CF conversion!")
            sys.exit(1)

        curdate += delta

def _proc_month(scriptdir, topdatadir, startdate, model_forcing):
    """Create the monthly CF file."""

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
    cmd = "%s/monthly_s2spost_nc.py" %(scriptdir)
    workdir = "./cf_%s_%4.4d%2.2d" %(model_forcing.upper(),
                                     startdate.year,
                                     startdate.month)
    cmd += " %s" %(workdir)
    cmd += " %s" %(workdir)
    cmd += " %4.4d%2.2d%2.2d" %(firstdate.year,
                                firstdate.month,
                                firstdate.day)
    cmd += " %4.4d%2.2d%2.2d" %(enddate.year,
                                enddate.month,
                                enddate.day)
    cmd += " %s" %(model_forcing.upper())

    print(cmd)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print("[ERR] Problem creating monthly file!")
        sys.exit(1)

def _create_done_file(startdate, model_forcing):
    """Create a 'done' file indicating batch job has finished."""
    path = "./cf_%s_%4.4d%2.2d/done" %(model_forcing.upper(),
                                       startdate.year,
                                       startdate.month)
    fobj = open(path, "w")
    fobj.close()

def _driver():
    """Main driver"""

    scriptdir, ldtfile, topdatadir, startdate, model_forcing = _read_cmd_args()
    _loop_daily(scriptdir, ldtfile, topdatadir, startdate, model_forcing)
    _proc_month(scriptdir, topdatadir, startdate, model_forcing)
    _create_done_file(startdate, model_forcing)

# Invoke driver
if __name__ == "__main__":
    _driver()
