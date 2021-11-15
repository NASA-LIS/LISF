#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_03.py
#
# PURPOSE: Reorganizes the downloaded nmme data into a format for further
# processing. Based on FORECAST_TASK_03.sh.
#
# REVISION HISTORY:
# 24 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
import os
import subprocess
import sys

# Local methods
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} month_num current_year config_file"
    print(txt)
    print("[INFO] where")
    print("[INFO] month_num: Current integer month")
    print("[INFO] current_year: Current year")
    print("[INFO] config_file: Config file that sets up environment")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    try:
        month_num = int(sys.argv[1])
    except ValueError:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)
    if month_num < 0:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)

    try:
        current_year = int(sys.argv[2])
    except ValueError:
        print(f"[ERR] Invalid argument for current_year! Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)
    if current_year < 0:
        print(f"[ERR] Invalid argument for current_year! Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)

    # config_file
    config_file = sys.argv[3]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return month_num, current_year, config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def _driver():
    """Main driver."""
    month_num, current_year, config_file = _read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

    # Path of the directory where all the BC codes are kept
    srcdir = config["bcsd_preproc"]["srcdir"]

    # Path of the directory where supplementary files are kept
    supplementary_dir = config["bcsd_preproc"]["supplementary_dir"]

    # Path for where raw and bias corrected forecast files are located:
    nmme_download_dir = config["bcsd_preproc"]["nmme_download_dir"]
    forcedir = f"{projdir}/data/forecast/NMME"
    nmme_output_dir=f"{forcedir}/raw/Monthly"

    cmd = "python"
    cmd += f" {srcdir}/nmme_reorg_f.py"
    cmd += f" {month_num}"
    cmd += f" {current_year}"
    cmd += f" {nmme_download_dir}"
    cmd += f" {nmme_output_dir}"
    cmd += f" {supplementary_dir}"
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print(f"[ERR] Problem calling python script: {cmd}.")
        sys.exit(1)

if __name__ == "__main__":
    _driver()
