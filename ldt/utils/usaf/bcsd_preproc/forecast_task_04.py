#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_04.py
#
# PURPOSE: Computes the bias correction for the forecast (CFSv2) dataset. Based
# on FORECAST_TASK_04.sh.
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
    txt = f"[INFO] Usage: {(sys.argv[0])} FCST_SYR FCST_EYR CLIM_SYR CLIM_EYR "\
        "month_abbr iMonNo lat1 lat2 lon1 lon2 lead_months ens_numc CONFIG_FILE"
    print(txt)
    print("[INFO] where")
    print("[INFO] FCST_SYR: Start year of forecast")
    print("[INFO] FCST_EYR: End year of forecast")
    print("[INFO] CLIM_SYR: Start year of the climatological period")
    print("[INFO] CLIM_EYR: End year of the climatological period")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] lat1: Minimum latitudinal extent")
    print("[INFO] lat2: Maximum latitudinal extent")
    print("[INFO] lon1: Minimum longitudinal extent")
    print("[INFO] lon2: Maximum longitudinal extent")
    print("[INFO] lead_months: Number of lead months")
    print("[INFO] ens_num: Number of ensembles")
    print("[INFO] CONFIG_FILE: Config file that sets up environment")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 14:
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

    # CLIM_SYR
    try:
        clim_syr = int(sys.argv[3])
    except ValueError:
        print(f"[ERR] Invalid argument for CLIM_SYR! Received {(sys.argv[3])}")
        _usage()
        sys.exit(1)
    if clim_syr < 0:
        print(f"[ERR] Invalid argument for CLIM_SYR! Received {(sys.argv[3])}")
        _usage()
        sys.exit(1)

    # CLIM_EYR
    try:
        clim_eyr = int(sys.argv[4])
    except ValueError:
        print(f"[ERR] Invalid argument for CLIM_EYR! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)
    if clim_eyr < 0:
        print(f"[ERR] Invalid argument for CLIM_EYR! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)

    # MONTH_ABBR
    month_abbr = str(sys.argv[5])

    # MONTH_NUM
    try:
        month_num = int(sys.argv[6])
    except ValueError:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[6])}")
        _usage()
        sys.exit(1)
    if month_num < 1:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[6])}")
        _usage()
        sys.exit(1)
    if month_num > 12:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[6])}")
        _usage()
        sys.exit(1)

    # LAT1
    try:
        lat1 = int(sys.argv[7])
    except ValueError:
        print(f"[ERR] Invalid argument for LAT1! Received {(sys.argv[7])}")
        _usage()
        sys.exit(1)

    # LAT2
    try:
        lat2 = int(sys.argv[8])
    except ValueError:
        print(f"[ERR] Invalid argument for LAT2! Received {(sys.argv[8])}")
        _usage()
        sys.exit(1)

    # LON1
    try:
        lon1 = int(sys.argv[9])
    except ValueError:
        print(f"[ERR] Invalid argument for LON1! Received {(sys.argv[9])}")
        _usage()
        sys.exit(1)

    # LON2
    try:
        lon2 = int(sys.argv[10])
    except ValueError:
        print(f"[ERR] Invalid argument for LON2! Received {(sys.argv[10])}")
        _usage()
        sys.exit(1)

    # LEAD_MONTHS
    try:
        lead_months = int(sys.argv[11])
    except ValueError:
        print(f"[ERR] Invalid argument for LEAD_MONTHS! Received {(sys.argv[11])}")
        _usage()
        sys.exit(1)
    if lead_months < 0:
        print(f"[ERR] Invalid argument for LEAD_MONTHS! Received {(sys.argv[11])}")
        _usage()
        sys.exit(1)

    # ENS_NUM
    try:
        ens_num = int(sys.argv[12])
    except ValueError:
        print(f"[ERR] Invalid argument for ENS_NUM! Received {(sys.argv[12])}")
        _usage()
        sys.exit(1)
    if ens_num < 0:
        print(f"[ERR] Invalid argument for ENS_NUM! Received {(sys.argv[12])}")
        _usage()
        sys.exit(1)

    # CONFIG_FILE
    CONFIG_FILE = sys.argv[13]
    if not os.path.exists(CONFIG_FILE):
        print(f"[ERR] {CONFIG_FILE} does not exist!")
        sys.exit(1)

    return fcst_syr, fcst_eyr, clim_syr, clim_eyr, month_abbr, month_num,\
    lat1, lat2, lon1, lon2, lead_months, ens_num, CONFIG_FILE

def read_config(CONFIG_FILE):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    return config

def _driver():
    """Main driver."""
    fcst_syr, fcst_eyr, clim_syr, clim_eyr, month_abbr, month_num, lat1, lat2, \
        lon1, lon2, lead_months, ens_num, CONFIG_FILE = _read_cmd_args()

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

    # Path for where observational & forecast files are located:
    FORCEDIR=f"{PROJDIR}/data/forecast"
    OBS_INDIR=f"{FORCEDIR}/USAF-LIS7.3rc8_25km"
    FCST_INDIR=f"{FORCEDIR}/CFSv2_25km"

    # Mask file
    MASK_FILE=f"{SUPPLEMENTARY_DIR}/Mask_nafpa.nc"

    #  Calculate bias correction for different variables separately:
    OBS_VAR_LIST=["Rainf_f_tavg", "LWdown_f_tavg", "SWdown_f_tavg", \
        "Psurf_f_tavg", "Qair_f_tavg", "Tair_f_tavg", "Wind_f_tavg"]
    FCST_VAR_LIST=["PRECTOT", "LWS", "SLRSF", "PS", "Q2M", "T2M", "WIND10M"]
    UNIT_LIST=["kg/m^2/s", "W/m^2", "W/m^2", "Pa", "kg/kg", "K", "m/s"]

    # BC output directory for FCST:
    outdir=f"{FCST_INDIR}/bcsd/Monthly/{month_abbr}01"

    print("[INFO] Processing forecast bias correction of CFSv2 variables")

#    for var_num in range(len(OBS_VAR_LIST)):
    for var_num,var_value in enumerate(OBS_VAR_LIST):
        if var_num in [0, 2]:
            var_type="PRCP"
        else:
            var_type="TEMP"

        obs_var=OBS_VAR_LIST[var_num]
        fcst_var=FCST_VAR_LIST[var_num]
        unit=UNIT_LIST[var_num]
        print(f"{var_num} {fcst_var}")

        cmd = "sbatch"
        cmd += f" {SRCDIR}/run_BCSD_calctest.scr"
        cmd += f" {SRCDIR}"
        cmd += f" {obs_var}"
        cmd += f" {fcst_var}"
        cmd += f" {month_num}"
        cmd += f" {var_type}"
        cmd += f" {unit}"
        cmd += f" {lat1}"
        cmd += f" {lat2}"
        cmd += f" {lon1}"
        cmd += f" {lon2}"
        cmd += f" {ens_num}"
        cmd += f" {lead_months}"
        cmd += f" {fcst_syr}"
        cmd += f" {fcst_eyr}"
        cmd += f" {clim_syr}"
        cmd += f" {clim_eyr}"
        cmd += f" {MASK_FILE}"
        cmd += f" {OBS_INDIR}"
        cmd += f" {FCST_INDIR}"
        cmd += f" {outdir}"
        cmd += f" {LOGDIR}"
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)

    print(f"[INFO] Completed processing forecast bias correction for: {month_abbr}")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
