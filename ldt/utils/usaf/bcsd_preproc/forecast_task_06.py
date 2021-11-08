#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_06.py
#
# PURPOSE: Generate bias-corrected 6-hourly forecasts using raw monthly
# forecasts, bias-corrected monthly forecasts and raw 6-hourly forecasts. Based
# on FORECAST_TASK_06.sh.
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
    txt = f"[INFO] Usage: {(sys.argv[0])} fcst_syr fcst_eyr month_abbr "\
    "month_num lat1 lat2 lon1 lon2 model_name lead_months ens_num config_file"
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] fcst_eyr: End year of forecast")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] lat1: Minimum latitudinal extent")
    print("[INFO] lat2: Maximum latitudinal extent")
    print("[INFO] lon1: Minimum longitudinal extent")
    print("[INFO] lon2: Maximum longitudinal extent")
    print("[INFO] model_name: Model name (should be CFSv2)")
    print("[INFO] lead_months: Number of lead months")
    print("[INFO] ens_num: Integer number of ensembles")
    print("[INFO] config_file: Config file that sets up environment")

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 13:
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
        print(f"[ERR] Invalid argument for fcst_eyr!  Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)
    if fcst_eyr < 0:
        print(f"[ERR] Invalid argument for fcst_eyr!  Received {(sys.argv[2])}")
        _usage()
        sys.exit(1)

    # month_abbr
    month_abbr = str(sys.argv[3])

    # month_num
    try:
        month_num = int(sys.argv[4])
    except ValueError:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)
    if month_num < 1:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)
    if month_num > 12:
        print(f"[ERR] Invalid argument for month_num! Received {(sys.argv[4])}")
        _usage()
        sys.exit(1)

    # lat1
    try:
        lat1 = int(sys.argv[5])
    except ValueError:
        print(f"[ERR] Invalid argument for lat1! Received {(sys.argv[5])}")
        _usage()
        sys.exit(1)

    # lat2
    try:
        lat2 = int(sys.argv[6])
    except ValueError:
        print(f"[ERR] Invalid argument for lat2! Received {(sys.argv[6])}")
        _usage()
        sys.exit(1)

    # lon1
    try:
        lon1 = int(sys.argv[7])
    except ValueError:
        print(f"[ERR] Invalid argument for lon1! Received {(sys.argv[7])}")
        _usage()
        sys.exit(1)

    # lon2
    try:
        lon2 = int(sys.argv[8])
    except ValueError:
        print(f"[ERR] Invalid argument for lon2! Received {(sys.argv[8])}")
        _usage()
        sys.exit(1)

    # model_name
    model_name = str(sys.argv[9])

    # lead_months
    try:
        lead_months = int(sys.argv[10])
    except ValueError:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[10])}")
        _usage()
        sys.exit(1)
    if lead_months < 0:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[10])}")
        _usage()
        sys.exit(1)

    # ens_num
    try:
        ens_num = int(sys.argv[11])
    except ValueError:
        print(f"[ERR] Invalid argument for ens_num! Received {(sys.argv[11])}")
        _usage()
        sys.exit(1)
    if ens_num < 0:
        print(f"[ERR] Invalid argument for ens_num! Received {(sys.argv[11])}")
        _usage()
        sys.exit(1)

    # config_file
    config_file = sys.argv[12]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return fcst_syr, fcst_eyr, month_abbr, month_num, lat1, lat2, lon1, lon2, \
    model_name, lead_months, ens_num, config_file

def _driver():
    """Main driver."""
    fcst_syr, fcst_eyr, month_abbr, month_num, lat1, lat2, lon1, lon2, \
    model_name, lead_months, ens_num, config_file = _read_cmd_args()

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

    # Path for where forecast files are located:
    forcedir=f"{projdir}/data/forecast/CFSv2_25km"

    # Mask file
    mask_file_precip=f"{supplementary_dir}/Mask_nafpa.nc"
    mask_file_nonprecip=f"{supplementary_dir}/Mask_nafpa.nc"

    #  Calculate bias correction for different variables separately:
    obs_var_list=["LWGAB", "SWGDN", "PS", "QV2M", "T2M", "U10M"]
    fcst_var_list=["LWS", "SLRSF", "PS", "Q2M", "T2M", "WIND10M"]
    unit_list=["W/m^2", "W/m^2", "Pa", "kg/kg", "K", "m/s"]

    # Path for where forecast and bias corrected files are located:
    subdaily_raw_fcst_dir=f"{forcedir}/raw/6-Hourly/{month_abbr}01"
    monthly_raw_fcst_dir=f"{forcedir}/raw/Monthly/{month_abbr}01"
    monthly_bc_fcst_dir=f"{forcedir}/bcsd/Monthly/{month_abbr}01"

    outdir=f"{forcedir}/bcsd/6-Hourly/{month_abbr}01"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("[INFO] Processing temporal disaggregation of CFSv2 variables")
    for year in range(fcst_syr, (fcst_eyr + 1)):
        for var_num, var_value in enumerate(obs_var_list):
            if var_num == 1:
                var_type="PRCP"
            else:
                var_type="TEMP"

            obs_var=var_value
            fcst_var=fcst_var_list[var_num]
            unit=unit_list[var_num]

            cmd = "sbatch"
            cmd += f" {srcdir}/run_Temporal_disagg.scr"
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
            cmd += f" {model_name}"
            cmd += f" {ens_num}"
            cmd += f" {lead_months}"
            cmd += f" {year}"
            cmd += f" {year}"
            cmd += f" {mask_file_precip}"
            cmd += f" {mask_file_nonprecip}"
            cmd += f" {monthly_bc_fcst_dir}"
            cmd += f" {monthly_raw_fcst_dir}"
            cmd += f" {subdaily_raw_fcst_dir}"
            cmd += f" {outdir}"
            cmd += f" {logdir}"
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem calling sbatch!")
                sys.exit(1)

    print(f"[INFO] Completed CFSv2 temporal disaggregation for: {(month_abbr)}")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
