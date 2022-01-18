#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_09.py
#
# PURPOSE: Combine all non-precip 6-hourly files into one file and copy BCSD
# precip files in to the same directory. Based on FORECAST_TASK_09.sh.
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
    txt = f"[INFO] Usage: {(sys.argv[0])} FCST_SYR FCST_EYR MONTH_ABBR "\
    	"MONTH_NUM FCST_TYPE LEAD_MONTHS ENS_NUM CONFIG_FILE"
    print(txt)
    print("[INFO] where")
    print("[INFO] FCST_SYR: Start year of forecast")
    print("[INFO] FCST_EYR: End year of forecast")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] MONTH_NUM: Integer number of the initialization month")
    print("[INFO] FCST_TYPE: Forecast type (should be nmme)")
    print("[INFO] LEAD_MONTHS: Number of lead months")
    print("[INFO] ENS_NUM: Integer number of ensembles")
    print("[INFO] CONFIG_FILE: Configfile that sets up environment")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 9:
        print("[ERR] Invalid number of command line arguments!")
        usage()
        sys.exit(1)

    # FCST_SYR
    try:
        fcst_syr = int(sys.argv[1])
    except ValueError:
        print(f"[ERR] Invalid argument for FCST_SYR! Received {(sys.argv[1])}")
        usage()
        sys.exit(1)
    if fcst_syr < 0:
        print(f"[ERR] Invalid argument for FCST_SYR! Received {(sys.argv[1])}")
        usage()
        sys.exit(1)

    # FCST_EYR
    try:
        fcst_eyr = int(sys.argv[2])
    except ValueError:
        print(f"[ERR] Invalid argument for FCST_EYR! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)
    if fcst_eyr < 0:
        print(f"[ERR] Invalid argument for FCST_EYR! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)

    # MONTH_ABBR
    month_abbr = str(sys.argv[3])

    # MONTH_NUM
    try:
        month_num = int(sys.argv[4])
    except ValueError:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[4])}")
        usage()
        sys.exit(1)
    if month_num < 1:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[4])}")
        usage()
        sys.exit(1)
    if month_num > 12:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[4])}")
        usage()
        sys.exit(1)

    # FCST_TYPE
    fcst_type = str(sys.argv[5])

    # LEAD_MONTHS
    try:
        lead_months = int(sys.argv[6])
    except ValueError:
        print(f"[ERR] Invalid argument for LEAD_MONTHS! Received {(sys.argv[6])}")
        usage()
        sys.exit(1)
    if lead_months < 0:
        print(f"[ERR] Invalid argument for LEAD_MONTHS! Received {(sys.argv[6])}")
        usage()
        sys.exit(1)

    # ENS_NUM
    try:
        ens_num = int(sys.argv[7])
    except ValueError:
        print(f"[ERR] Invalid argument for ENS_NUM!  Received {(sys.argv[7])}")
        usage()
        sys.exit(1)
    if ens_num < 0:
        print(f"[ERR] Invalid argument for ENS_NUM!  Received {(sys.argv[7])}")
        usage()
        sys.exit(1)

    # CONFIG_FILE
    config_file = sys.argv[8]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return fcst_syr, fcst_eyr, month_abbr, month_num, fcst_type, lead_months, \
    	ens_num, config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def driver():
    """Main driver."""
    fcst_syr, fcst_eyr, month_abbr, month_num, fcst_type, lead_months, \
    	ens_num, config_file = read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

    # Path of the directory where all the BC codes are kept:
    srcdir = config["bcsd_preproc"]["srcdir"]

    # Log file output directory
    logdir = config["bcsd_preproc"]["logdir"]

    # Path for the final 6-hourly forcing data:
    forcedir=f"{projdir}/data/forecast/CFSv2_25km"

    print("[INFO] Combining subdaily BC CFSv2 non-precip variables")
    for year in range(fcst_syr, (fcst_eyr + 1)):
        cmd = "sbatch"
        cmd += f" {srcdir}/run_Combining.scr"
        cmd += f" {srcdir}"
        cmd += f" {year}"
        cmd += f" {month_num}"
        cmd += f" {ens_num}"
        cmd += f" {lead_months}"
        cmd += f" {forcedir}"
        cmd += f" {fcst_type}"
        cmd += f" {logdir}"
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)

    print(f"[INFO] Completed CFSv2 combination for: {month_abbr}")

#
# Main Method
#
if __name__ == "__main__":
    driver()
