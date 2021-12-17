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

import configparser
import os
import subprocess
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} month_abbr month_num current_year "\
        "ens_num lead_months config_file"
    print(txt)
    print("[INFO] where")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] current_year: Current year of forecast")
    print("[INFO] ens_num: Integer number of ensembles")
    print("[INFO] lead_months: Number of lead months")
    print("[INFO] config_file: Config file that sets up environment")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 7:
        print("[ERR] Invalid number of command line arguments!")
        usage()
        sys.exit(1)

    # month_abbr
    month_abbr = str(sys.argv[1])

    # month_num
    try:
        month_num = int(sys.argv[2])
    except ValueError:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)
    if month_num < 1:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)
    if month_num > 12:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)

    # current_year
    try:
        current_year = int(sys.argv[3])
    except ValueError:
        print(f"[ERR] Invalid argument for current_year! Received {(sys.argv[3])}")
        usage()
        sys.exit(1)
    if current_year < 0:
        print(f"[ERR] Invalid argument for current_year! Received {(sys.argv[3])}")
        usage()
        sys.exit(1)

    # ens_num
    try:
        ens_num = int(sys.argv[4])
    except ValueError:
        print(f"[ERR] Invalid argument for ens_num! Received {(sys.argv[4])}")
        usage()
        sys.exit(1)
    if ens_num < 0:
        print(f"[ERR] Invalid argument for ens_num! Received {(sys.argv[4])}")
        usage()
        sys.exit(1)

    # lead_months
    try:
        lead_months = int(sys.argv[5])
    except ValueError:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[5])}")
        usage()
        sys.exit(1)
    if lead_months < 0:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[5])}")
        usage()
        sys.exit(1)

    # config_file
    config_file = sys.argv[6]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return month_abbr, month_num, current_year, ens_num, lead_months, \
        config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def driver():
    """Main driver."""
    month_abbr, month_num, current_year, ens_num, lead_months, \
        config_file = read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

    # Path for where forecast files are located:
    forcedir=f"{projdir}/data/forecast/CFSv2_25km/final/6-Hourly"

    for iens in range(1, (ens_num + 1)):
        print(f"Ensemble {iens}/{ens_num}")

        indir=f"{forcedir}/{current_year}/{month_abbr}01/ens{iens}"

        for ilead in range(lead_months):
            fcst_date=(datetime(current_year, month_num, 1) + \
                relativedelta(months=ilead)).strftime("%Y%m")

            srcfile=f"{indir}/CFSv2.{fcst_date}.nc4"
            dstfile=srcfile

            cmd = "ncap2 -O -s 'V10M=array(0.0,0.0,U10M)'"
            cmd += f" {srcfile}"
            cmd += f" {dstfile}"
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem creating V10M variable!")
                sys.exit(1)

#
# Main Method
#
if __name__ == "__main__":
    driver()
