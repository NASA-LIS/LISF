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
    txt = f"[INFO] Usage: {(sys.argv[0])} fcst_syr fcst_eyr clim_syr clim_eyr "\
    "month_abbr month_num lat1 lat2 lon1 lon2 lead_months ens_num config_file"
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] fcst_eyr: End year of forecast")
    print("[INFO] clim_syr: Start year of the climatological period")
    print("[INFO] clim_eyr: End year of the climatological period")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] lat1: Minimum latitudinal extent")
    print("[INFO] lat2: Maximum latitudinal extent")
    print("[INFO] lon1: Minimum longitudinal extent")
    print("[INFO] lon2: Maximum longitudinal extent")
    print("[INFO] lead_months: Number of lead months")
    print("[INFO] ens_num: Number of ensembles")
    print("[INFO] config_file: Config file that sets up environment")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 14:
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

    # clim_syr
    try:
        clim_syr = int(sys.argv[3])
    except ValueError:
        print(f"[ERR] Invalid argument for clim_syr! Received {(sys.argv[3])}")
        _usage()
        sys.exit(1)
    if clim_syr < 0:
        print(f"[ERR] Invalid argument for clim_syr! Received {(sys.argv[3])}")
        _usage()
        sys.exit(1)

    # clim_eyr
    try:
        clim_eyr = int(sys.argv[4])
    except ValueError:
        print(f"[ERR] Invalid argument for clim_eyr! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)
    if clim_eyr < 0:
        print(f"[ERR] Invalid argument for clim_eyr! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)

    # month_abbr
    month_abbr = str(sys.argv[5])

    # month_num
    try:
        month_num = int(sys.argv[6])
    except ValueError:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[6])}")
        _usage()
        sys.exit(1)
    if month_num < 1:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[6])}")
        _usage()
        sys.exit(1)
    if month_num > 12:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[6])}")
        _usage()
        sys.exit(1)

    # lat1
    try:
        lat1 = int(sys.argv[7])
    except ValueError:
        print(f"[ERR] Invalid argument for lat1! Received {(sys.argv[7])}")
        _usage()
        sys.exit(1)

    # lat2
    try:
        lat2 = int(sys.argv[8])
    except ValueError:
        print(f"[ERR] Invalid argument for lat2! Received {(sys.argv[8])}")
        _usage()
        sys.exit(1)

    # lon1
    try:
        lon1 = int(sys.argv[9])
    except ValueError:
        print(f"[ERR] Invalid argument for lon1! Received {(sys.argv[9])}")
        _usage()
        sys.exit(1)

    # lon2
    try:
        lon2 = int(sys.argv[10])
    except ValueError:
        print(f"[ERR] Invalid argument for lon2! Received {(sys.argv[10])}")
        _usage()
        sys.exit(1)

    # lead_months
    try:
        lead_months = int(sys.argv[11])
    except ValueError:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[11])}")
        _usage()
        sys.exit(1)
    if lead_months < 0:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[11])}")
        _usage()
        sys.exit(1)

    # ens_num
    try:
        ens_num = int(sys.argv[12])
    except ValueError:
        print(f"[ERR] Invalid argument for ens_num! Received {(sys.argv[12])}")
        _usage()
        sys.exit(1)
    if ens_num < 0:
        print(f"[ERR] Invalid argument for ens_num! Received {(sys.argv[12])}")
        _usage()
        sys.exit(1)

    # config_file
    config_file = sys.argv[13]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return fcst_syr, fcst_eyr, clim_syr, clim_eyr, month_abbr, month_num,\
    lat1, lat2, lon1, lon2, lead_months, ens_num, config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def _driver():
    """Main driver."""
    fcst_syr, fcst_eyr, clim_syr, clim_eyr, month_abbr, month_num, lat1, lat2, \
        lon1, lon2, lead_months, ens_num, config_file = _read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

    # Path of the directory where all the BC codes are kept
    srcdir = config["bcsd_preproc"]["srcdir"]

    # Log file output directory
    logdir = config["bcsd_preproc"]["logdir"]

    # Path of the directory where supplementary files are kept
    supplementary_dir = config["bcsd_preproc"]["supplementary_dir"]

    # Path for where observational & forecast files are located:
    forcedir=f"{projdir}/data/forecast"
    obs_indir=f"{forcedir}/USAF-LIS7.3rc8_25km"
    fcst_indir=f"{forcedir}/CFSv2_25km"

    # Mask file
    mask_file=f"{supplementary_dir}/Mask_nafpa.nc"

    #  Calculate bias correction for different variables separately:
    obs_var_list=["Rainf_f_tavg", "LWdown_f_tavg", "SWdown_f_tavg", \
        "Psurf_f_tavg", "Qair_f_tavg", "Tair_f_tavg", "Wind_f_tavg"]
    fcst_var_list=["PRECTOT", "LWS", "SLRSF", "PS", "Q2M", "T2M", "WIND10M"]
    unit_list=["kg/m^2/s", "W/m^2", "W/m^2", "Pa", "kg/kg", "K", "m/s"]

    # BC output directory for FCST:
    outdir=f"{fcst_indir}/bcsd/Monthly/{month_abbr}01"

    print("[INFO] Processing forecast bias correction of CFSv2 variables")

#    for var_num in range(len(obs_var_list)):
    for var_num,var_value in enumerate(obs_var_list):
        if var_num in [0, 2]:
            var_type="PRCP"
        else:
            var_type="TEMP"

        obs_var=var_value
        fcst_var=fcst_var_list[var_num]
        unit=unit_list[var_num]
        print(f"{var_num} {fcst_var}")

        cmd = "sbatch"
        cmd += f" {srcdir}/run_BCSD_calctest.scr"
        cmd += f" {srcdir}"
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
        cmd += f" {mask_file}"
        cmd += f" {obs_indir}"
        cmd += f" {fcst_indir}"
        cmd += f" {outdir}"
        cmd += f" {logdir}"
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
