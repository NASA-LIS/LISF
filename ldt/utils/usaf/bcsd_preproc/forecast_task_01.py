#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_01.py
#
# PURPOSE: Processes the CFSv2 forecast data and outputs in 6-hourly and
# monthly time resolutions.  Based on FORECAST_TASK_01.sh.
#
# REVISION HISTORY:
# 22 Oct 2021: Eric Kemp/SSAI, first version
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
    txt = f"[INFO] Usage: {(sys.argv[0])} fcst_syr fcst_eyr month_abbr "\
        "config_file"
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] fcst_eyr: End year of forecast")
    print("[INFO] month_abbr: Abbreviated month to start forecast")
    print("[INFO] config_file: Config file that sets up environment")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 5:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # fcst_syr
    try:
        fcst_syr = int(sys.argv[1])
    except ValueError:
        print(f"[ERR] Invalid argument for fcst_syr! Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)
    if fcst_syr < 0:
        print(f"[ERR] Invalid argument for fcst_syr! Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)

    # fcst_eyr
    try:
        fcst_eyr = int(sys.argv[2])
    except ValueError:
        print(f"[ERR] Invalid argument for fcst_eyr! Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)
    if fcst_eyr < 0:
        print(f"[ERR] Invalid argument for fcst_eyr! Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)

    # month_abbr
    month_abbr = sys.argv[3]

    # config_file
    config_file = sys.argv[4]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return fcst_syr, fcst_eyr, month_abbr, config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def calc_ic_dates(icmon):
    """Generates forecast initialization dates based on the initialization
    month."""

    # We'll store the dates in a dictionary, and then pull the appropriate
    # selection based on the initialization month code.
    ic_dates_all = {
        "jan01" : ['1217', '1222', '1227'],
        "feb01" : ['0121', '0126', '0131'],
        "mar01" : ['0215', '0220', '0225'],
        "apr01" : ['0317', '0322', '0327'],
        "may01" : ['0416', '0421', '0426'],
        "jun01" : ['0521', '0526', '0531'],
        "jul01" : ['0620', '0625', '0630'],
        "aug01" : ['0720', '0725', '0730'],
        "sep01" : ['0819', '0824', '0829'],
        "oct01" : ['0918', '0923', '0928'],
        "nov01" : ['1018', '1023', '1028'],
        "dec01" : ['1117', '1122', '1127'],
    }
    try:
        ic_dates = ic_dates_all[icmon]
    except KeyError:
        print(f"[ERR] Unknown initialization month {icmon}")
        sys.exit(1)
    return ic_dates

def _driver():
    """Main driver."""
    fcst_syr, fcst_eyr, month_abbr, config_file = _read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

    # Path of the directory where all the BC codes are kept
    srcdir = config["bcsd_preproc"]["srcdir"]

    # Path of the directory where patch files for missing data are kept
    patchdir = config["bcsd_preproc"]["patchdir"]

    # Path of the directory where supplementary files are kept
    supplementary_dir = config["bcsd_preproc"]["supplementary_dir"]

    # Log file output directory
    logdir = config["bcsd_preproc"]["logdir"]

    # Paths for the daily forecast data (input and output paths)
    forcedir = config["bcsd_preproc"]["fcst_download_dir"]
    outdir = f"{projdir}/data/forecast/CFSv2_25km/raw"
    griddesc = f"{supplementary_dir}/CFSv2_25km_AFRICOM_grid_description.txt"

    if not os.path.exists(logdir):
        os.makedirs(logdir)

    imon = f"{month_abbr}01"
    ic_dates = calc_ic_dates(imon)

    # Process 3-hrly CFSv2 forecasts and output in monthly and 6-hrly formats
    print("[INFO] Processing CFSv2 3-hrly forecast variables")
    for year in range(fcst_syr, (fcst_eyr + 1)):
        cmd = "sbatch"
        cmd += f" {srcdir}/run_process_forecast_data.scr"
        cmd += f" {year:04d}"
        cmd += f" {year:04d}"
        cmd += f" {imon}"
        cmd += f" {srcdir}"
        cmd += f" {outdir}"
        cmd += f" {forcedir}"
        cmd += f" {griddesc}"
        cmd += f" {patchdir}"
        for ic_date in ic_dates:
            cmd += f" {ic_date}"
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)
    print(f"[INFO] Jobs submitted to process CFSv2 forecast files for {imon}")

if __name__ == "__main__":
    _driver()
