#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_01.py
#
# PURPOSE: Processes the CFSv2 forecast data and outputs in 6-hourly and
# monthly time resolutions.  Based on FORECAST_TASK_01.sh.
#
# REVISION HISTORY:
# 22 Oct 2021: Eric Kemp/SSAI, first version
#
#------------------------------------------------------------------------------
"""

# Standard modules
import os
import subprocess
import sys

# Local constants.  FIXME:  Put in single location for whole system

# Path of the main project directory
_PROJDIR = '/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'

# Path of the directory where all the download codes are kept:
_SRCDIR = "%s/scripts/code_library" %(_PROJDIR)

# Path of the directory to where the patch files for missing data are kept:
_PATCHDIR = "%s/data/CFSv2_25km/patch_files" %(_PROJDIR)

# Paths for the daily forecast data (input and output paths):
_FORCEDIR = '/discover/nobackup/projects/lis/MET_FORCING/CFSv2'
_OUTDIR = "%s/data/forecast/CFSv2_25km/raw" %(_PROJDIR)
_GRIDDESC = "%s/supplementary_files/CFSv2_25km_AFRICOM_grid_description.txt" \
    %(_SRCDIR)

# Log file output directory
_LOGDIR = "%s/scripts/log_files" %(_PROJDIR)

# Local methods
def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s fcst_syr fcst_eyr month_abbr" %(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] fcst_eyr: End year of forecast")
    print("[INFO] month_abbr: Abbreviated month to start forecast")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    try:
        fcst_syr = int(sys.argv[1])
    except ValueError:
        print("[ERR] Invalid argument for fcst_syr!  Received %s" \
              %(sys.argv[1]))
        _usage()
        sys.exit(1)
    if fcst_syr < 0:
        print("[ERR] Invalid argument for fcst_syr!  Received %s" \
              %(sys.argv[1]))
        _usage()
        sys.exit(1)

    try:
        fcst_eyr = int(sys.argv[2])
    except ValueError:
        print("[ERR] Invalid argument for fcst_eyr!  Received %s" \
              %(sys.argv[2]))
        _usage()
        sys.exit(1)
    if fcst_eyr < 0:
        print("[ERR] Invalid argument for fcst_eyr!  Received %s" \
              %(sys.argv[2]))
        _usage()
        sys.exit(1)

    month_abbr = sys.argv[3]

    return fcst_syr, fcst_eyr, month_abbr

def calc_ic_dates(icmon):
    """Generates forecast initialization dates based on the initialization
    month."""

    # We'll store the dates in a dictionary, and then pull the appropriate
    # selection based on the initialization month code.
    ic_dates_all = {
        "jan01" : ['1217', '1222', '1227'],
        "feb01" : ['0121', '0126', '0131'],
        "mar01" : ['0215', '0220', '0225'],
        "apr01" : ['0317', '0322', '0327'],
        "may01" : ['0416', '0421', '0426'],
        "jun01" : ['0521', '0526', '0531'],
        "jul01" : ['0620', '0625', '0630'],
        "aug01" : ['0720', '0725', '0730'],
        "sep01" : ['0819', '0824', '0829'],
        "oct01" : ['0918', '0923', '0928'],
        "nov01" : ['1018', '1023', '1028'],
        "dec01" : ['1117', '1122', '1127'],
    }
    try:
        ic_dates = ic_dates_all[icmon]
    except KeyError:
        print("[ERR] Unknown initialization month %s" %(icmon))
        sys.exit(1)
    return ic_dates

def _driver():
    """Main driver."""
    if not os.path.exists(_LOGDIR):
        os.makedirs(_LOGDIR)
    fcst_syr, fcst_eyr, month_abbr = _read_cmd_args()
    imon = "%s01" %(month_abbr)
    ic_dates = calc_ic_dates(imon)

    # Process 3-hrly CFSv2 forecasts and output in monthly and 6-hrly formats
    print("[INFO] Processing CFSv2 3-hrly forecast variables")
    for year in range(fcst_syr, (fcst_eyr + 1)):
        cmd = "sbatch"
        cmd += " {srcdir}/run_process_forecast_data.scr".format(srcdir=_SRCDIR)
        cmd += " {year}".format(year=year)
        cmd += " {year}".format(year=year)
        cmd += " {imon}".format(imon=imon)
        cmd += " {srcdir}".format(srcdir=_SRCDIR)
        cmd += " {outdir}".format(outdir=_OUTDIR)
        cmd += " {forcedir}".format(forcedir=_FORCEDIR)
        cmd += " {griddesc}".format(griddesc=_GRIDDESC)
        cmd += " {patchdir}".format(patchdir=_PATCHDIR)
        for ic_date in ic_dates:
            cmd += " %s" %(ic_date)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)
    print("[INFO] Jobs submitted to process CFSv2 forecast files for %s" \
          %(imon))

if __name__ == "__main__":
    _driver()
