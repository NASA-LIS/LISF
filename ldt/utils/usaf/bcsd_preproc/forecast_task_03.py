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
    txt = f"[INFO] Usage: {(sys.argv[0])} MONTH_NUM CURRENT_YEAR CONFIG_FILE"
    print(txt)
    print("[INFO] where")
    print("[INFO] MONTH_NUM: Current integer month")
    print("[INFO] CURRENT_YEAR: Current year")
    print("[INFO] CONFIG_FILE: Config file that sets up environment")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    try:
        MONTH_NUM = int(sys.argv[1])
    except ValueError:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)
    if MONTH_NUM < 0:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)

    try:
        CURRENT_YEAR = int(sys.argv[2])
    except ValueError:
        print(f"[ERR] Invalid argument for CURRENT_YEAR! Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)
    if CURRENT_YEAR < 0:
        print(f"[ERR] Invalid argument for CURRENT_YEAR! Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)

    # CONFIG_FILE
    CONFIG_FILE = sys.argv[3]
    if not os.path.exists(CONFIG_FILE):
        print(f"[ERR] {CONFIG_FILE} does not exist!")
        sys.exit(1)

    return MONTH_NUM, CURRENT_YEAR, CONFIG_FILE

def read_config(CONFIG_FILE):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    return config

def _driver():
    """Main driver."""
    MONTH_NUM, CURRENT_YEAR, CONFIG_FILE = _read_cmd_args()

    # Setup local directories
    config = read_config(CONFIG_FILE)

    # Path of the main project directory
    PROJDIR = config["bcsd_preproc"]["projdir"]

    # Path of the directory where all the BC codes are kept
    SRCDIR = config["bcsd_preproc"]["srcdir"]

    # Path of the directory where supplementary files are kept
    SUPPLEMENTARY_DIR = config["bcsd_preproc"]["supplementary_dir"]

    # Path for where raw and bias corrected forecast files are located:
    NMME_DOWNLOAD_DIR = config["bcsd_preproc"]["nmme_download_dir"]
    FORCEDIR = f"{PROJDIR}/data/forecast/NMME"
    NMME_OUTPUT_DIR=f"{FORCEDIR}/raw/Monthly"

    cmd = "python"
    cmd += f" {SRCDIR}/nmme_reorg_f.py"
    cmd += f" {MONTH_NUM}"
    cmd += f" {CURRENT_YEAR}"
    cmd += f" {NMME_DOWNLOAD_DIR}"
    cmd += f" {NMME_OUTPUT_DIR}"
    cmd += f" {SUPPLEMENTARY_DIR}"
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print(f"[ERR] Problem calling python script: {cmd}.")
        sys.exit(1)

if __name__ == "__main__":
    _driver()

