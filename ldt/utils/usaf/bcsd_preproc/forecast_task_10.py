#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_10.py
#
# PURPOSE: Copy BC NMME precip files in to the same directory. Based on
# FORECAST_TASK_10.sh.
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
#import calendar
from dateutil.relativedelta import relativedelta
#
# Local constants.  FIXME:  Put in single location for whole system
#

# Path of the main project directory
PROJDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'

# Path of the directory where all the BC codes are kept:
SRCDIR="{}/scripts/code_library".format(PROJDIR)

#  Log file output directory
LOGDIR="{}/scripts/log_files".format(PROJDIR)

# Path for the final 6-hourly forcing dataets:
FORCEDIR_FCST="{}/data/forecast/CFSv2_25km".format(PROJDIR)
FORCEDIR_NMME="{}/data/forecast/NMME".format(PROJDIR)

# Base model prefixes for forecast files
BASEMODNAME_SRC="CFSv2"
BASEMODNAME_DST="GEOS5"

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} MONTH_ABBR MONTH_NUM CURRENT_YEAR NMME_MODEL"\
        .format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] MONTH_NUM: Integer number of the initialization month")
    print("[INFO] CURRENT_YEAR: Current year of forecast")
    print("[INFO] NMME_MODEL: NMME model name")

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

    # NMME_MODEL
    NMME_MODEL = str(sys.argv[4])

    return MONTH_ABBR, MONTH_NUM, CURRENT_YEAR, NMME_MODEL

def gather_ensemble_info(NMME_MODEL):
    """Gathers ensemble information based on NMME model."""

    # Number of ensembles in the forecast (ENS_NUM)
    # Ensemble start index (ENS_START)
    # Ensemble end index (ENS_END)
    if NMME_MODEL == "CFSv2":
        ENS_NUM=24
        ENS_START=1
        ENS_END=24
        ENS_RANGE=list(range(1, 13)) + list(range(1,13))
    elif NMME_MODEL == "GEOSv2":
        ENS_NUM=10
        ENS_START=25
        ENS_END=34
        ENS_RANGE=range(1, 11)
    elif NMME_MODEL == "CCM4":
        ENS_NUM=10
        ENS_START=35
        ENS_END=44
        ENS_RANGE=range(1, 11)
    elif NMME_MODEL == "GNEMO":
        ENS_NUM=10
        ENS_START=45
        ENS_END=54
        ENS_RANGE=range(1, 11)
    elif NMME_MODEL == "CCSM4":
        ENS_NUM=10
        ENS_START=55
        ENS_END=64
        ENS_RANGE=range(1, 11)
    elif NMME_MODEL == "GFDL":
        ENS_NUM=30
        ENS_START=65
        ENS_END=94
        ENS_RANGE=list(range(1, 13)) + list(range(1,13)) + list(range(1,7))
    else:
        print("[ERR] Invalid argument for NMME_MODEL!  Received {}" \
            .format(NMME_MODEL))
        sys.exit(1)

    return ENS_NUM, ENS_START, ENS_END, ENS_RANGE

def copy_subdaily_precipitation(YEAR, MONTH_ABBR, ENS_NUM, INDIR_NMME, OUTDIR):
    """Copies the BC 6-Hourly precipition files to the final directory."""

    for iens in range(1, (ENS_NUM + 1)):
        print("[INFO] Copying precipitation files for: {}, ens {}/{}"\
            .format(YEAR, iens, ENS_NUM))

        INDIR_COMPLETE="{}/{}/ens{}".format(INDIR_NMME, YEAR, iens)

        OUTDIR_COMPLETE="{}/{}/{}01/ens{}"\
            .format(OUTDIR, YEAR, MONTH_ABBR, iens)

        if not os.path.exists(OUTDIR_COMPLETE):
            os.makedirs(OUTDIR_COMPLETE)

        cmd = "cp {}/PRECTOT* {}/".format(INDIR_COMPLETE, OUTDIR_COMPLETE)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling copy subroutine!")
            sys.exit(1)

def driver():
    """Main driver."""
    MONTH_ABBR, MONTH_NUM, CURRENT_YEAR, NMME_MODEL = read_cmd_args()

    init_datetime = datetime(CURRENT_YEAR, MONTH_NUM, 1)
    SRC_YYYYMM = [(init_datetime + relativedelta(months=0)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=1)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=2)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=3)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=4)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=5)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=6)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=7)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=8)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=8)).strftime("%Y%m")]

    DST_YYYYMM = [(init_datetime + relativedelta(months=0)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=1)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=2)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=3)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=4)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=5)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=6)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=7)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=8)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=9)).strftime("%Y%m")]

    ENS_NUM, ENS_START, ENS_END, ENS_RANGE = gather_ensemble_info(NMME_MODEL)

    INDIR_NMME="{}/bcsd/6-Hourly/{}01/{}"\
        .format(FORCEDIR_NMME, MONTH_ABBR, NMME_MODEL)
    OUTDIR="{}/final/6-Hourly/{}".format(FORCEDIR_NMME, NMME_MODEL)

    # Copy the precipitation files
    print("[INFO] NMME MODEL: {}".format(NMME_MODEL))
    copy_subdaily_precipitation(CURRENT_YEAR, MONTH_ABBR, ENS_NUM, INDIR_NMME,\
         OUTDIR)

    # Symbolically link the non-precip data
    print("[INFO] Creating symbolic links for non-precip data")
    for iens in range(len(ENS_RANGE)):
        ENS_NMME=iens + 1
        ENS_FCST=ENS_RANGE[iens]

        SRC_DIR="{}/final/6-Hourly/{}/{}01/ens{}"\
            .format(FORCEDIR_FCST, CURRENT_YEAR, MONTH_ABBR, ENS_FCST)
        DST_DIR="{}/{}/{}01/ens{}"\
            .format(OUTDIR, CURRENT_YEAR, MONTH_ABBR, ENS_NMME)

        for ilead in range(len(SRC_YYYYMM)):
            SRC_FILE="{}/{}.{}.nc4"\
                .format(SRC_DIR, BASEMODNAME_SRC, SRC_YYYYMM[ilead])
            DST_FILE="{}/{}.{}.nc4"\
                .format(DST_DIR, BASEMODNAME_DST, DST_YYYYMM[ilead])

            cmd = "ln -sfn {src_file} {dst_file}"\
                .format(src_file=SRC_FILE, dst_file=DST_FILE)
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem calling creating symbolic links!")
                sys.exit(1)

#
# Main Method
#
if __name__ == "__main__":
    driver()
