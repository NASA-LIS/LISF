#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

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

import subprocess
import sys
import argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
import yaml

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {(sys.argv[0])} -s current_year -m month_abbr \
                         -w cwd -n month_num -c config_file"
    print(txt)
    print("[INFO] current_year: Start year of forecast")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] config_file: Config file that sets up environment")
    print("[INFO] cwd: current working directory")

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
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--current_year', required=True, help='forecast start year')
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-m', '--month_abbr', required=True, help='month abbreviation')
    parser.add_argument('-n', '--month_num', required=True, help='month number')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')

    args = parser.parse_args()
    config_file = args.config_file
    current_year = int(args.current_year)
    month_abbr = args.month_abbr
    month_num = int(args.month_num)
    cwd = args.cwd

    # load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)
    
    lead_months = config['EXP']['lead_months']

    # Path of the main project directory
    projdir = cwd

    # Path of the final nmme directory
    nmme_data_dir = f"{projdir}/bcsd_fcst/NMME/final/6-Hourly"

    # Array of all NMME models
    nmme_models = config['EXP']['NMME_models']
    ensemble_sizes = config['EXP']['ensemble_sizes'][0]
    src_date, dst_date = gather_date_info(current_year, month_num, lead_months)

    for nmme_model in nmme_models:
        ens_num = ensemble_sizes[nmme_model]

        for member in range(1, ens_num+1):
            outdir = f"{nmme_data_dir}/{nmme_model}/{month_abbr}01/{current_year}/ens{member}"

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
