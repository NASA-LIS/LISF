#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_12.py
#
# PURPOSE: Prepares an all zero V10M variable for LIS preparation due to the
# USAF-LIS observational forcing only including average windspeed. Based on
# FORECAST_TASK_12.sh.
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
PROJDIR="/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM"

# Path for where forecast files are located:
FORCEDIR="{}/data/forecast/NMME/final/6-Hourly".format(PROJDIR)

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} MONTH_ABBR CURRENT_YEAR NMME_MODEL LEAD_MONTHS"\
        .format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] CURRENT_YEAR: Current year of forecast")
    print("[INFO] NMME_MODEL: NMME model name")
    print("[INFO] LEAD_MONTHS: Number of lead months")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 5:
        print("[ERR] Invalid number of command line arguments!")
        usage()
        sys.exit(1)

    # MONTH_ABBR
    MONTH_ABBR = str(sys.argv[1])

    # CURRENT_YEAR
    try:
        CURRENT_YEAR = int(sys.argv[2])
    except ValueError:
        print("[ERR] Invalid argument for CURRENT_YEAR!  Received {}" \
            .format(sys.argv[2]))
        usage()
        sys.exit(1)
    if CURRENT_YEAR < 0:
        print("[ERR] Invalid argument for CURRENT_YEAR!  Received {}" \
              .format(sys.argv[2]))
        usage()
        sys.exit(1)

    # NMME_MODEL
    NMME_MODEL = str(sys.argv[3])

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

    return MONTH_ABBR, CURRENT_YEAR, NMME_MODEL, LEAD_MONTHS

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

def driver():
    """Main driver."""

    MONTH_ABBR, CURRENT_YEAR, NMME_MODEL, LEAD_MONTHS = read_cmd_args()

    ENS_NUM = gather_ensemble_info(NMME_MODEL)

    FILE_COUNT=2*(LEAD_MONTHS + 1)*ENS_NUM
    FILE_SIZE=26*ENS_NUM

    FINAL_DIR="{}/{}/{}/{}01"\
        .format(FORCEDIR, NMME_MODEL, CURRENT_YEAR, MONTH_ABBR)

    print("{}: {} {}".format(NMME_MODEL, MONTH_ABBR, CURRENT_YEAR))
    print("Expected file count is: {}".format(FILE_COUNT))
    print("Actual file count is:")
    cmd = "ls {}/ens*/* | wc -l".format(FINAL_DIR)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print("[ERR] Problem calling file count subroutine!")
        sys.exit(1)

    print("Expected directory size is approximately: {} MB".format(FILE_SIZE))
    print("Actual directory size is:")
    cmd = "du -sm {} | cut -f1".format(FINAL_DIR)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print("[ERR] Problem calling file size subroutine!")
        sys.exit(1)
    print("")

#
# Main Method
#
if __name__ == "__main__":
    driver()
