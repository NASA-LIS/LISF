#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_11.py
#
# PURPOSE: Copy files to fill final "10-month" for writing out full average for
# initialized forecasts. Based on FORECAST_TASK_11.sh.
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
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

#
# Local constants.  FIXME:  Put in single location for whole system
#

# Path of the main project directory
PROJDIR="/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM"

# Path of the final nmme directory
NMME_DATA_DIR="{projdir}/data/forecast/NMME/final/6-Hourly"\
    .format(projdir=PROJDIR)

# Array of all NMME models
NMME_MODELS=["CFSv2", "GEOSv2", "CCSM4", "CCM4", "GNEMO", "GFDL"]

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} MONTH_ABBR MONTH_NUM CURRENT_YEAR"\
        .format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] MONTH_NUM: Integer number of the initialization month")
    print("[INFO] CURRENT_YEAR: Current year of forecast")
    print("[INFO] LEAD_MONTHS: Number of lead months")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 5:
        print("[ERR] Invalid number of command line arguments!")
        usage()
        sys.exit(1)

    # MONTH_ABBR
    MONTH_ABBR = str(sys.argv[1])

    # MONTH_NUM
    try:
        MONTH_NUM = int(sys.argv[2])
    except ValueError:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[2]))
        usage()
        sys.exit(1)
    if MONTH_NUM < 1:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[2]))
        usage()
        sys.exit(1)
    if MONTH_NUM > 12:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[2]))
        usage()
        sys.exit(1)

    # CURRENT_YEAR
    try:
        CURRENT_YEAR = int(sys.argv[3])
    except ValueError:
        print("[ERR] Invalid argument for CURRENT_YEAR!  Received {}" \
            .format(sys.argv[3]))
        usage()
        sys.exit(1)
    if CURRENT_YEAR < 0:
        print("[ERR] Invalid argument for CURRENT_YEAR!  Received {}" \
              .format(sys.argv[3]))
        usage()
        sys.exit(1)

    # LEAD_MONTHS
    try:
        LEAD_MONTHS = int(sys.argv[4])
    except ValueError:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[4]))
        usage()
        sys.exit(1)
    if LEAD_MONTHS < 0:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[4]))
        usage()
        sys.exit(1)

    return MONTH_ABBR, MONTH_NUM, CURRENT_YEAR, LEAD_MONTHS

def gather_ensemble_info(NMME_MODEL):
    """Gathers ensemble information based on NMME model."""

    # Number of ensembles in the forecast (ENS_NUM)
    if NMME_MODEL == "CFSv2":
        ENS_NUM=24
    elif NMME_MODEL == "GEOSv2":
        ENS_NUM=10
    elif NMME_MODEL == "CCM4":
        ENS_NUM=10
    elif NMME_MODEL == "GNEMO":
        ENS_NUM=10
    elif NMME_MODEL == "CCSM4":
        ENS_NUM=10
    elif NMME_MODEL == "GFDL":
        ENS_NUM=30
    else:
        print("[ERR] Invalid argument for NMME_MODEL!  Received {}" \
            .format(NMME_MODEL))
        sys.exit(1)

    return ENS_NUM

def gather_date_info(CURRENT_YEAR, MONTH_NUM, LEAD_MONTHS):
    """Gathers monthly date information based on fcst and lead months."""

    INIT_DATETIME = datetime(CURRENT_YEAR, MONTH_NUM, 1)
    SRC_DATETIME = INIT_DATETIME + relativedelta(months=(LEAD_MONTHS-1))
    DST_DATETIME = INIT_DATETIME + relativedelta(months=(LEAD_MONTHS))

    SRC_DATE = SRC_DATETIME.strftime("%Y%m")
    DST_DATE = DST_DATETIME.strftime("%Y%m")

    return SRC_DATE, DST_DATE

def driver():
    """Main driver."""
    MONTH_ABBR, MONTH_NUM, CURRENT_YEAR, LEAD_MONTHS = read_cmd_args()

    SRC_DATE, DST_DATE = gather_date_info(CURRENT_YEAR, MONTH_NUM, LEAD_MONTHS)

    for NMME_MODEL in NMME_MODELS:
        ENS_NUM = gather_ensemble_info(NMME_MODEL)

        for MEMBER in range(1,ENS_NUM+1):
            DIR="{}/{}/{}/{}01/ens{}".format(NMME_DATA_DIR, NMME_MODEL, \
                CURRENT_YEAR, MONTH_ABBR, MEMBER)

            cmd = "cp"
            cmd += " {}/PRECTOT.{}.nc4".format(DIR, SRC_DATE)
            cmd += " {}/PRECTOT.{}.nc4".format(DIR, DST_DATE)
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem calling copy subroutine!")
                sys.exit(1)

#
# Main Method
#
if __name__ == "__main__":
    driver()


