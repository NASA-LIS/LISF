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
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

#
# Local constants.  FIXME:  Put in single location for whole system
#

# Path of the main project directory
PROJDIR="/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM"

# Path for where forecast files are located:
FORCEDIR="{}/data/forecast/CFSv2_25km/final/6-Hourly".format(PROJDIR)

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} MONTH_ABBR MONTH_NUM CURRENT_YEAR ENS_NUM"\
        .format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] MONTH_NUM: Integer number of the initialization month")
    print("[INFO] CURRENT_YEAR: Current year of forecast")
    print("[INFO] ENS_NUM: Integer number of ensembles")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 6:
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

    # ENS_NUM
    try:
        ENS_NUM = int(sys.argv[4])
    except ValueError:
        print("[ERR] Invalid argument for ENS_NUM!  Received {}" \
              .format(sys.argv[4]))
        usage()
        sys.exit(1)
    if ENS_NUM < 0:
        print("[ERR] Invalid argument for ENS_NUM!  Received {}" \
              .format(sys.argv[4]))
        usage()
        sys.exit(1)

    # LEAD_MONTHS
    try:
        LEAD_MONTHS = int(sys.argv[5])
    except ValueError:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[5]))
        usage()
        sys.exit(1)
    if LEAD_MONTHS < 0:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[5]))
        usage()
        sys.exit(1)

    return MONTH_ABBR, MONTH_NUM, CURRENT_YEAR, ENS_NUM, LEAD_MONTHS

def driver():
    """Main driver."""
    MONTH_ABBR, MONTH_NUM, CURRENT_YEAR, ENS_NUM, LEAD_MONTHS = read_cmd_args()

    for iens in range(1, (ENS_NUM + 1)):
        print("Ensemble {}/{}".format(iens, ENS_NUM))

        INDIR="{}/{}/{}01/ens{}"\
            .format(FORCEDIR, CURRENT_YEAR, MONTH_ABBR, iens)

        for ilead in range(LEAD_MONTHS):
            FCST_DATE=(datetime(CURRENT_YEAR, MONTH_NUM, 1) + \
                relativedelta(months=ilead)).strftime("%Y%m")

            SRCFILE="{}/CFSv2.{}.nc4".format(INDIR, FCST_DATE)
            DSTFILE=SRCFILE

            cmd = "ncap2 -O -s 'V10M=array(0.0,0.0,U10M)'"
            cmd += " {}".format(SRCFILE)
            cmd += " {}".format(DSTFILE)
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem creating V10M variable!")
                sys.exit(1)

#
# Main Method
#
if __name__ == "__main__":
    driver()
