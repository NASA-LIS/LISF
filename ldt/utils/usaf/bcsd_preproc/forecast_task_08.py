#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_08.py
#
# PURPOSE: Generate bias-corrected 6-hourly nmme forecasts using raw monthly
# forecasts, bias-corrected monthly forecasts and raw 6-hourly forecasts. Based
# on FORECAST_TASK_08.sh.
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

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} FCST_SYR FCST_EYR MONTH_ABBR "\
    "MONTH_NUM LAT1 LAT2 LON1 LON2 NMME_MODEL LEAD_MONTHS CONFIG_FILE"
    print(txt)
    print("[INFO] where")
    print("[INFO] FCST_SYR: Start year of forecast")
    print("[INFO] FCST_EYR: End year of forecast")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] MONTH_NUM: Integer number of the initialization month")
    print("[INFO] LAT1: Minimum latitudinal extent")
    print("[INFO] LAT2: Maximum latitudinal extent")
    print("[INFO] LON1: Minimum longitudinal extent")
    print("[INFO] LON2: Maximum longitudinal extent")
    print("[INFO] NMME_MODEL: NMME model name")
    print("[INFO] LEAD_MONTHS: Number of lead months")
    print("[INFO] CONFIG_FILE: Config file that sets up environment")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 12:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # FCST_SYR
    try:
        fcst_syr = int(sys.argv[1])
    except ValueError:
        print(f"[ERR] Invalid argument for FCST_SYR! Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)
    if fcst_syr < 0:
        print(f"[ERR] Invalid argument for FCST_SYR! Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)

    # FCST_EYR
    try:
        fcst_eyr = int(sys.argv[2])
    except ValueError:
        print(f"[ERR] Invalid argument for FCST_EYR! Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)
    if fcst_eyr < 0:
        print(f"[ERR] Invalid argument for FCST_EYR! Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)

    # MONTH_ABBR
    month_abbr = str(sys.argv[3])

    # MONTH_NUM
    try:
        month_num = int(sys.argv[4])
    except ValueError:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)
    if month_num < 1:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)
    if month_num > 12:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)

    # LAT1
    try:
        lat1 = int(sys.argv[5])
    except ValueError:
        print(f"[ERR] Invalid argument for LAT1! Received {(sys.argv[5])}")
        _usage()
        sys.exit(1)

    # LAT2
    try:
        lat2 = int(sys.argv[6])
    except ValueError:
        print(f"[ERR] Invalid argument for LAT2! Received {(sys.argv[6])}")
        _usage()
        sys.exit(1)

    # LON1
    try:
        lon1 = int(sys.argv[7])
    except ValueError:
        print(f"[ERR] Invalid argument for LON1! Received {(sys.argv[7])}")
        _usage()
        sys.exit(1)

    # LON2
    try:
        lon2 = int(sys.argv[8])
    except ValueError:
        print(f"[ERR] Invalid argument for LON2! Received {(sys.argv[8])}")
        _usage()
        sys.exit(1)

    # NMME_MODEL
    nmme_model = str(sys.argv[9])

    # LEAD_MONTHS
    try:
        lead_months = int(sys.argv[10])
    except ValueError:
        print(f"[ERR] Invalid argument for LEAD_MONTHS! Received {(sys.argv[10])}")
        _usage()
        sys.exit(1)
    if lead_months < 0:
        print(f"[ERR] Invalid argument for LEAD_MONTHS! Received {(sys.argv[10])}")
        _usage()
        sys.exit(1)

    # CONFIG_FILE
    CONFIG_FILE = sys.argv[11]
    if not os.path.exists(CONFIG_FILE):
        print(f"[ERR] {CONFIG_FILE} does not exist!")
        sys.exit(1)

    return fcst_syr, fcst_eyr, month_abbr, month_num, lat1, lat2, lon1, lon2, \
        nmme_model, lead_months, CONFIG_FILE

def read_config(CONFIG_FILE):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    return config

def gather_ensemble_info(nmme_model):
    """Gathers ensemble information based on NMME model."""

    # Number of ensembles in the forecast (ens_num)
    # Ensemble start index (ens_start)
    # Ensemble end index (ens_end)
    if nmme_model == "CFSv2":
        ens_num=24
        ens_start=1
        ens_end=24
    elif nmme_model == "GEOSv2":
        ens_num=10
        ens_start=25
        ens_end=34
    elif nmme_model == "CCM4":
        ens_num=10
        ens_start=35
        ens_end=44
    elif nmme_model == "GNEMO":
        ens_num=10
        ens_start=45
        ens_end=54
    elif nmme_model == "CCSM4":
        ens_num=10
        ens_start=55
        ens_end=64
    elif nmme_model == "GFDL":
        ens_num=30
        ens_start=65
        ens_end=94
    else:
        print(f"[ERR] Invalid argument for NMME_MODEL! Received {nmme_model}")
        sys.exit(1)

    return ens_num, ens_start, ens_end

def _driver():
    """Main driver."""
    fcst_syr, fcst_eyr, month_abbr, month_num, lat1, lat2, lon1, lon2, \
    nmme_model, lead_months, CONFIG_FILE = read_cmd_args()

    # Setup local directories
    config = read_config(CONFIG_FILE)

    # Path of the main project directory
    PROJDIR = config["bcsd_preproc"]["projdir"]

    # Path of the directory where all the BC codes are kept
    SRCDIR = config["bcsd_preproc"]["srcdir"]

    # Log file output directory
    LOGDIR = config["bcsd_preproc"]["logdir"]

    # Path of the directory where supplementary files are kept
    SUPPLEMENTARY_DIR = config["bcsd_preproc"]["supplementary_dir"]

    # Path for where forecast files are located:
    FORCEDIR=f"{PROJDIR}/data/forecast/NMME"

    # Mask file
    MASK_FILE_PRECIP=f"{SUPPLEMENTARY_DIR}/Mask_nafpa.nc"
    MASK_FILE_NONPRECIP=f"{SUPPLEMENTARY_DIR}/Mask_nafpa.nc"

    #  Calculate bias correction for different variables separately:
    OBS_VAR="PRECTOT"
    FCST_VAR="PRECTOT"
    UNIT="kg/m^2/s"
    VAR_TYPE='PRCP'

    # Path for where forecast and bias corrected files are located:
    subdaily_raw_fcst_dir=f"{FORCEDIR}/linked_cfsv2_precip_files/{month_abbr}01"
    monthly_raw_fcst_dir=f"{FORCEDIR}/raw/Monthly/{month_abbr}01"
    monthly_bc_fcst_dir="{FORCEDIR}/bcsd/Monthly/{month_abbr}01"

    outdir="{FORCEDIR}/bcsd/6-Hourly/{month_abbr}01/{nmme_model}"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ens_num, ens_start, ens_end = gather_ensemble_info(nmme_model)
    print("[INFO] Processing temporal disaggregation of CFSv2 variables")
    for year in range(fcst_syr, (fcst_eyr + 1)):
        cmd = "sbatch"
        cmd += f" {SRCDIR}/run_NMME_Temporal_disagg.scr"
        cmd += f" {SRCDIR}"
        cmd += f" {OBS_VAR}"
        cmd += f" {FCST_VAR}"
        cmd += f" {month_num}"
        cmd += f" {VAR_TYPE}"
        cmd += f" {UNIT}"
        cmd += f" {lat1}"
        cmd += f" {lat2}"
        cmd += f" {lon1}"
        cmd += f" {lon2}"
        cmd += f" {nmme_model}"
        cmd += f" {ens_num}"
        cmd += f" {lead_months}"
        cmd += f" {year}"
        cmd += f" {year}"
        cmd += f" {MASK_FILE_PRECIP}"
        cmd += f" {MASK_FILE_NONPRECIP}"
        cmd += f" {monthly_bc_fcst_dir}"
        cmd += f" {monthly_raw_fcst_dir}"
        cmd += f" {subdaily_raw_fcst_dir}"
        cmd += f" {outdir}"
        cmd += f" {LOGDIR}"
        cmd += f" {ens_start}"
        cmd += f" {ens_end}"
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)

    print(f"[INFO] Completed NMME temporal disaggregation for: {month_abbr}")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
