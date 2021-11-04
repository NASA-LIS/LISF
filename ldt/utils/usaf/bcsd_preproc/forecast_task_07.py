#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_07.py
#
# PURPOSE: Combine all non-precip 6-hourly files into one file and copy BCSD
# precip files in to the same directory Based on FORECAST_TASK_07.sh.
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

# Number of precip ensembles needed
RANGE_ENS_FCST=list(range(1, 13)) + list(range(1,13)) + list(range(1,7))
RANGE_ENS_NMME=range(1,31)

#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} CURRENTYEAR MON".format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] CURRENTYEAR: Current year")
    print("[INFO] MONTH_ABBR: Current month")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 3:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # CURRENT_YEAR
    try:
        CURRENT_YEAR = int(sys.argv[1])
    except ValueError:
        print("[ERR] Invalid argument for CURRENT_YEAR!  Received {}" \
            .format(sys.argv[1]))
        _usage()
        sys.exit(1)
    if CURRENT_YEAR < 0:
        print("[ERR] Invalid argument for CURRENT_YEAR!  Received {}" \
              .format(sys.argv[1]))
        _usage()
        sys.exit(1)

    # MONTH_ABBR
    MONTH_ABBR = str(sys.argv[2])

    return CURRENT_YEAR, MONTH_ABBR

def _driver():
    """Main driver."""
    CURRENT_YEAR, MONTH_ABBR = _read_cmd_args()
    FCST_DATE = "{month_abbr}01".format(month_abbr=MONTH_ABBR)

    # Path for where forecast files are located:
    INDIR="{}/data/forecast/CFSv2_25km/raw/6-Hourly/{}/{}" \
        .format(PROJDIR, FCST_DATE, CURRENT_YEAR)

    # Path for where the linked precip files should be placed:
    OUTDIR="{}/data/forecast/NMME/linked_cfsv2_precip_files/{}/{}" \
        .format(PROJDIR, FCST_DATE, CURRENT_YEAR)

    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    for iens in range(len(RANGE_ENS_FCST)):
        SRC_FILE="{}/ens{}".format(INDIR, RANGE_ENS_FCST[iens])
        DST_FILE="{}/ens{}".format(OUTDIR, RANGE_ENS_NMME[iens])

        cmd = "ln -sfn {src_file} {dst_file}" \
            .format(src_file=SRC_FILE, dst_file=DST_FILE)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling creating symbolic links!")
            sys.exit(1)

    print("[INFO] Done creating symbolic links")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
