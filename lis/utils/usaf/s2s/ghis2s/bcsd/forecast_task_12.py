#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

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
    ens_num = config['BCSD']['nof_raw_ens']

    # Path of the main project directory
    projdir = cwd

    # Path for where forecast files are located:
    forcedir = f"{projdir}/bcsd_fcst/CFSv2_25km/final/6-Hourly"

    for iens in range(1, (ens_num + 1)):
        print(f"Ensemble {iens}/{ens_num}")

        indir = f"{forcedir}/{month_abbr}01/{current_year}/ens{iens}"

        for ilead in range(lead_months):
            fcst_date = (datetime(current_year, month_num, 1) + \
                relativedelta(months=ilead)).strftime("%Y%m")

            srcfile = f"{indir}/CFSv2.{fcst_date}.nc4"
            dstfile = srcfile

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
