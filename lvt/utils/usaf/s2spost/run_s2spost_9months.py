#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: run_s2spost_9months.py
#
# PURPOSE: Loops through multiple months of a completed LIS S2S forecast, and
# submits batch jobs (one per month) to generate CF-convention netCDF files.
#
# REVISION HISTORY:
# 24 Sep 2021: Eric Kemp (SSAI), first version.
# 27 Oct 2021: Eric Kemp/SSAI, address pylint string objections.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import glob
import os
import subprocess
import shutil
import sys
import time

_TOTAL_MONTHS = 9

def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} ldt_file topdatadir startYYYYMM"
    txt += " model_forcing"
    txt += " [--collect_output]"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]   ldt_file is path to LDT parameter file")
    print("[INFO]   topdatadir is top-level directory with LIS data")
    print("[INFO]   startYYYYMM is year-month of start of LIS forecast")
    print("[INFO]   model_forcing is ID for atmospheric forcing for LIS")
    txt = "[INFO]   [--collect_output] is option to collect files into"
    txt += " single directory"
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) not in [5, 6]:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Assume all scripts are bundled together in the dirname of this script.
    scriptdir = os.path.dirname(sys.argv[0])

    # Get path to LDT parameter file
    ldtfile = sys.argv[1]
    if not os.path.exists(ldtfile):
        print(f"[ERR] LDT parameter file {ldtfile} does not exist!")
        sys.exit(1)

    # Get top directory of LIS data
    topdatadir = sys.argv[2]
    if not os.path.exists(topdatadir):
        print(f"[ERR] LIS data directory {topdatadir} does not exist!")
        sys.exit(1)

    # Get valid year and month
    yyyymm = sys.argv[3]
    if len(yyyymm) != 6:
        print("[ERR] Invalid length of startYYYYMM, must be 6 characters!")
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

    # Handle option to consolidate files
    collect_output = False
    if len(sys.argv) == 6:
        if sys.argv[5] == "--collect_output":
            collect_output = True
        else:
            print(f"[ERR] Unknown argument {sys.argv[5]}")
            _usage()
            sys.exit(1)

    return scriptdir, ldtfile, topdatadir, startdate, model_forcing, \
        collect_output

def _advance_date_by_month(curdate):
    """Calculate new date one month in advance."""
    if curdate.month == 12:
        newdate = datetime.date(year=(curdate.year + 1),
                                month=1,
                                day=1)
    else:
        newdate = datetime.date(year=curdate.year,
                                month=(curdate.month + 1),
                                day=1)
    return newdate

def _submit_batch_jobs(scriptdir, ldtfile, topdatadir, startdate,
                       model_forcing):
    """Submit batch jobs for processing LIS forecast."""

    # Loop over all months
    curdate = startdate
    for _ in range(0, _TOTAL_MONTHS):
        txt = "[INFO] Submitting batch job for"
        txt += f" cf_{model_forcing}_{curdate.year:04d}{curdate.month:04d}"
        print(txt)
        cmd = f"sbatch {scriptdir}/run_s2spost_1month.sh {ldtfile}"
        cmd += f" {topdatadir}"
        cmd += f" {curdate.year:04d}{curdate.month:02d} {model_forcing}"
        #print(cmd)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem running run_s2spost_1month.sh")
            sys.exit(1)
        time.sleep(1)

        newdate = _advance_date_by_month(curdate)
        curdate = newdate

def _check_batch_job_completion(topdatadir, startdate, model_forcing):
    """Check for markers indicating batch jobs are completed."""

    print("[INFO] Checking for completion of batch jobs...")

    # Loop over all months
    curdate = startdate
    for _ in range(0, _TOTAL_MONTHS):

        subdir = f"{topdatadir}/cf_{model_forcing.upper()}"
        subdir += f"_{curdate.year:04d}{curdate.month:02d}"
        print(f"[INFO] Waiting for {subdir}")
        while True:
            if os.path.exists(f"{subdir}/done"):
                txt = f"[INFO] {model_forcing.upper()}"
                txt += f"_{curdate.year:04d}{curdate.month:02d}"
                print(txt)
                break
        newdate = _advance_date_by_month(curdate)
        curdate = newdate

def _consolidate_files(topdatadir, startdate, model_forcing):
    """Move CF files into common directory for single model forcing."""

    newdir = f"{topdatadir}/cf_{model_forcing.upper()}"
    newdir += f"_{startdate.year:04d}{startdate.month:02d}_all"
    if not os.path.exists(newdir):
        os.makedirs(newdir)

    # Loop over all months
    curdate = startdate
    for _ in range(0, _TOTAL_MONTHS):
        subdir = f"{topdatadir}/cf_{model_forcing.upper()}"
        subdir += f"_{curdate.year:04d}{curdate.month:02d}"
        print(f"[INFO] Copying {subdir} files to {newdir}" %(subdir, newdir))
        files = glob.glob(f"{subdir}/*.NC" %(subdir))
        for filename in files:
            shutil.copy(filename, newdir)
        #shutil.rmtree(subdir)

        newdate = _advance_date_by_month(curdate)
        curdate = newdate

def _driver():
    """Main driver."""
    scriptdir, ldtfile, topdatadir, startdate, model_forcing, \
        collect_output = _read_cmd_args()
    _submit_batch_jobs(scriptdir, ldtfile, topdatadir, startdate,
                       model_forcing)
    _check_batch_job_completion(topdatadir, startdate, model_forcing)
    if collect_output:
        _consolidate_files(topdatadir, startdate, model_forcing)

# Invoke driver
if __name__ == "__main__":
    _driver()
