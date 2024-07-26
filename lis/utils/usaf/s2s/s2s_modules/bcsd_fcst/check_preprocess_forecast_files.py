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
# SCRIPT: check_preprocess_forecast_files.py
#
# PURPOSE: Check the final file size and count of output. Based on
# Check_Preprocess_Finalfiles.F.sh.
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

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt =f"[INFO] Usage: {(sys.argv[0])} month_abbr current_year nmme_model lead_months config_file"
    print(txt)
    print("[INFO] where")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] current_year: Current year of forecast")
    print("[INFO] nmme_model: NMME model name")
    print("[INFO] lead_months: Number of lead months")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 6:
        print("[ERR] Invalid number of command line arguments!")
        usage()
        sys.exit(1)

    # month_abbr
    month_abbr = str(sys.argv[1])

    # current_year
    try:
        current_year = int(sys.argv[2])
    except ValueError:
        print(f"[ERR] Invalid argument for current_year! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)
    if current_year < 0:
        print(f"[ERR] Invalid argument for current_year! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)

    # nmme_model
    nmme_model = str(sys.argv[3])

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

    return month_abbr, current_year, nmme_model, lead_months, config_file

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
    elif nmme_model == "GNEMO5":
        ens_num=10
    elif nmme_model == "CCSM4":
        ens_num=10
    elif nmme_model == "GFDL":
        ens_num=30
    else:
        print(f"[ERR] Invalid argument for nmme_model!  Received {nmme_model}")
        sys.exit(1)

    return ens_num

def driver():
    """Main driver."""
    month_abbr, current_year, nmme_model, lead_months, config_file = read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

    # Path for where forecast files are located:
    forcedir=f"{projdir}/bcsd_fcst/forecast/NMME/final/6-Hourly"

    ens_num = gather_ensemble_info(nmme_model)

    file_count=2*(lead_months + 1)*ens_num
    file_size=26*ens_num

    final_dir=f"{forcedir}/{nmme_model}/{current_year}/{month_abbr}01"

    print(f"{nmme_model}: {month_abbr} {current_year}")
    print(f"Expected file count is: {file_count}")
    print("Actual file count is:")
    cmd = f"ls {final_dir}/ens*/* | wc -l"
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print("[ERR] Problem calling file count subroutine!")
        sys.exit(1)

    print(f"Expected directory size is approximately: {file_size} MB")
    print("Actual directory size is:")
    cmd = f"du -sm {final_dir} | cut -f1"
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
