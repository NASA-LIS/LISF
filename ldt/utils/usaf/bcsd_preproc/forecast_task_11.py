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
    "lead_months config_file"
    print(txt)
    print("[INFO] where")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] current_year: Current year of forecast")
    print("[INFO] lead_months: Number of lead months")
    print("[INFO] config_file: Config file that sets up environment")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 6:
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

    # lead_months
    try:
        lead_months = int(sys.argv[4])
    except ValueError:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[4])}")
        usage()
        sys.exit(1)
    if lead_months < 0:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[4])}")
        usage()
        sys.exit(1)

    # config_file
    config_file = sys.argv[5]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return month_abbr, month_num, current_year, lead_months, config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def gather_ensemble_info(nmme_model):
    """Gathers ensemble information based on NMME model."""

    # Number of ensembles in the forecast (ens_num)
    if nmme_model == "CFSv2":
        ens_num=24
    elif nmme_model == "GEOSv2":
        ens_num=10
    elif nmme_model == "CCM4":
        ens_num=10
    elif nmme_model == "GNEMO":
        ens_num=10
    elif nmme_model == "CCSM4":
        ens_num=10
    elif nmme_model == "GFDL":
        ens_num=30
    else:
        print(f"[ERR] Invalid argument for nmme_model!  Received {nmme_model}")
        sys.exit(1)

    return ens_num

def gather_date_info(current_year, month_num, lead_months):
    """Gathers monthly date information based on fcst and lead months."""

    init_datetime = datetime(current_year, month_num, 1)
    src_datetime = init_datetime + relativedelta(months=(lead_months-1))
    dst_datetime = init_datetime + relativedelta(months=(lead_months))

    src_date = src_datetime.strftime("%Y%m")
    dst_date = dst_datetime.strftime("%Y%m")

    return src_date, dst_date

def driver():
    """Main driver."""
    month_abbr, month_num, current_year, lead_months, \
        config_file = read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

    # Path of the final nmme directory
    nmme_data_dir=f"{projdir}/data/forecast/NMME/final/6-Hourly"

    # Array of all NMME models
    nmme_models=["CFSv2", "GEOSv2", "CCSM4", "CCM4", "GNEMO", "GFDL"]
    src_date, dst_date = gather_date_info(current_year, month_num, lead_months)

    for nmme_model in nmme_models:
        ens_num = gather_ensemble_info(nmme_model)

        for member in range(1,ens_num+1):
            outdir=f"{nmme_data_dir}/{nmme_model}/{current_year}/{month_abbr}01/ens{member}"

            cmd = "cp"
            cmd += f" {outdir}/PRECTOT.{src_date}.nc4"
            cmd += f" {outdir}/PRECTOT.{dst_date}.nc4"
            print(cmd)
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem calling copy subroutine!")
                sys.exit(1)

#
# Main Method
#
if __name__ == "__main__":
    driver()
