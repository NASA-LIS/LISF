#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_09.py
#
# PURPOSE: Combine all non-precip 6-hourly files into one file and copy BCSD
# precip files in to the same directory. Based on FORECAST_TASK_09.sh.
#
# REVISION HISTORY:
# 24 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""

#
# Standard modules
#

import os
import subprocess
import sys

#
# Local constants.  FIXME:  Put in single location for whole system
#

# Path of the main project directory
PROJDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'

# Path of the directory where all the BC codes are kept:
SRCDIR="{}/scripts/code_library".format(PROJDIR)

#  Log file output directory
LOGDIR="{}/scripts/log_files".format(PROJDIR)

# Path for the final 6-hourly forcing data:
FORCEDIR="{}/data/forecast/CFSv2_25km".format(PROJDIR)

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} FCST_SYR FCST_EYR MONTH_ABBR MONTH_NUM LAT1 LAT2 "\
        "LON1 LON2 FCST_TYPE LEAD_MONTHS ENS_NUM".format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] FCST_SYR: Start year of forecast")
    print("[INFO] FCST_EYR: End year of forecast")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] MONTH_NUM: Integer number of the initialization month")
    print("[INFO] LAT1: Minimum latitudinal extent")
    print("[INFO] LAT2: Maximum latitudinal extent")
    print("[INFO] LON1: Minimum longitudinal extent")
    print("[INFO] LON2: Maximum longitudinal extent")
    print("[INFO] FCST_TYPE: Forecast type (should be nmme)")
    print("[INFO] LEAD_MONTHS: Number of lead months")
    print("[INFO] ENS_NUM: Integer number of ensembles")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 12:
        print("[ERR] Invalid number of command line arguments!")
        usage()
        sys.exit(1)

    # FCST_SYR
    try:
        FCST_SYR = int(sys.argv[1])
    except ValueError:
        print("[ERR] Invalid argument for FCST_SYR!  Received {}" \
            .format(sys.argv[1]))
        usage()
        sys.exit(1)
    if FCST_SYR < 0:
        print("[ERR] Invalid argument for FCST_SYR!  Received {}" \
              .format(sys.argv[1]))
        usage()
        sys.exit(1)

    # FCST_EYR
    try:
        FCST_EYR = int(sys.argv[2])
    except ValueError:
        print("[ERR] Invalid argument for FCST_EYR!  Received {}" \
              .format(sys.argv[2]))
        usage()
        sys.exit(1)
    if FCST_EYR < 0:
        print("[ERR] Invalid argument for FCST_EYR!  Received {}" \
              .format(sys.argv[2]))
        usage()
        sys.exit(1)

    # MONTH_ABBR
    MONTH_ABBR = str(sys.argv[3])

    # MONTH_NUM
    try:
        MONTH_NUM = int(sys.argv[4])
    except ValueError:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[4]))
        usage()
        sys.exit(1)
    if MONTH_NUM < 1:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[4]))
        usage()
        sys.exit(1)
    if MONTH_NUM > 12:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[4]))
        usage()
        sys.exit(1)

    # LAT1
    try:
        LAT1 = int(sys.argv[5])
    except ValueError:
        print("[ERR] Invalid argument for LAT1!  Received {}" \
              .format(sys.argv[5]))
        usage()
        sys.exit(1)

    # LAT2
    try:
        LAT2 = int(sys.argv[6])
    except ValueError:
        print("[ERR] Invalid argument for LAT2!  Received {}" \
              .format(sys.argv[6]))
        usage()
        sys.exit(1)

    # LON1
    try:
        LON1 = int(sys.argv[7])
    except ValueError:
        print("[ERR] Invalid argument for LON1!  Received {}" \
              .format(sys.argv[7]))
        usage()
        sys.exit(1)

    # LON2
    try:
        LON2 = int(sys.argv[8])
    except ValueError:
        print("[ERR] Invalid argument for LON2!  Received {}" \
              .format(sys.argv[8]))
        usage()
        sys.exit(1)

    # FCST_TYPE
    FCST_TYPE = str(sys.argv[9])

    # LEAD_MONTHS
    try:
        LEAD_MONTHS = int(sys.argv[10])
    except ValueError:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[10]))
        usage()
        sys.exit(1)
    if LEAD_MONTHS < 0:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[10]))
        usage()
        sys.exit(1)

    # ENS_NUM
    try:
        ENS_NUM = int(sys.argv[11])
    except ValueError:
        print("[ERR] Invalid argument for ENS_NUM!  Received {}" \
              .format(sys.argv[11]))
        usage()
        sys.exit(1)
    if ENS_NUM < 0:
        print("[ERR] Invalid argument for ENS_NUM!  Received {}" \
              .format(sys.argv[11]))
        usage()
        sys.exit(1)

    return FCST_SYR, FCST_EYR, MONTH_ABBR, MONTH_NUM, LAT1, LAT2, LON1, LON2, \
    FCST_TYPE, LEAD_MONTHS, ENS_NUM

def driver():
    """Main driver."""
    FCST_SYR, FCST_EYR, MONTH_ABBR, MONTH_NUM, LAT1, LAT2, LON1, LON2, \
    FCST_TYPE, LEAD_MONTHS, ENS_NUM = read_cmd_args()

    print("[INFO] Combining subdaily BC CFSv2 non-precip variables")
    for YEAR in range(FCST_SYR, (FCST_EYR + 1)):
        cmd = "sbatch"
        cmd += " {srcdir}/run_Combining.scr".format(srcdir=SRCDIR)
        cmd += " {srcdir}".format(srcdir=SRCDIR)
        cmd += " {year}".format(year=YEAR)
        cmd += " {month_num}".format(month_num=MONTH_NUM)
        cmd += " {ens_num}".format(ens_num=ENS_NUM)
        cmd += " {lead_months}".format(lead_months=LEAD_MONTHS)
        cmd += " {forcedir}".format(forcedir=FORCEDIR)
        cmd += " {fcst_type}".format(fcst_type=FCST_TYPE)
        cmd += " {logdir}".format(logdir=LOGDIR)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)

    print("[INFO] Completed CFSv2 combination for: {}".format(MONTH_ABBR))

#
# Main Method
#
if __name__ == "__main__":
    driver()

