#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_05.py
#
# PURPOSE: Computes the bias correction for the NMME dataset. Based on
# FORECAST_TASK_03.sh.
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
    txt = "[INFO] Usage: {(sys.argv[0])} fcst_syr fcst_eyr clim_syr clim_eyr "\
    	"month_abbr month_num lat1 lat2 lon1 lon2 nmme_model lead_months "\
    	"config_file"
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
    print("[INFO] nmme_model: NMME model name")
    print("[INFO] lead_months: Number of lead months")
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

    # nmme_model
    nmme_model = str(sys.argv[11])

    # lead_months
    try:
        lead_months = int(sys.argv[12])
    except ValueError:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[12])}")
        _usage()
        sys.exit(1)
    if lead_months < 0:
        print(f"[ERR] Invalid argument for lead_months! Received {(sys.argv[12])}")
        _usage()
        sys.exit(1)

    # config_file
    config_file = sys.argv[13]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return fcst_syr, fcst_eyr, clim_syr, clim_eyr, month_abbr, month_num,\
    lat1, lat2, lon1, lon2, nmme_model, lead_months, config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def _gather_ensemble_info(nmme_model):
    """Gathers ensemble information based on NMME model."""

    # Number of ensembles in the forecast (ENS_NUMF)
    # Number of ensembles in the climatology (ENS_NUMC)
    # Ensemble start index (ENS_START)
    # Ensemble end index (ENS_END)
    if nmme_model == "CFSv2":
        ens_numf=24
        ens_numc=12
        ens_start=1
        ens_end=24
    elif nmme_model == "GEOSv2":
        ens_numf=10
        ens_numc=4
        ens_start=25
        ens_end=34
    elif nmme_model == "CCM4":
        ens_numf=10
        ens_numc=10
        ens_start=35
        ens_end=44
    elif nmme_model == "GNEMO":
        ens_numf=10
        ens_numc=10
        ens_start=45
        ens_end=54
    elif nmme_model == "CCSM4":
        ens_numf=10
        ens_numc=10
        ens_start=55
        ens_end=64
    elif nmme_model == "GFDL":
        ens_numf=30
        ens_numc=15
        ens_start=65
        ens_end=94
    else:
        print(f"[ERR] Invalid argument for nmme_model! Received {(nmme_model)}")
        sys.exit(1)

    return ens_numf, ens_numc, ens_start, ens_end

def _driver():
    """Main driver."""
    fcst_syr, fcst_eyr, clim_syr, clim_eyr, month_abbr, month_num, lat1, lat2,\
        lon1, lon2, nmme_model, lead_months, config_file = _read_cmd_args()

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

	# Path for where observational files are located:
    forcedir=f"{projdir}/data/forecast"
    obs_clim_indir=f"{forcedir}/USAF-LIS7.3rc8_25km/raw/Climatology"

	# Mask file
    mask_file=f"{supplementary_dir}/Mask_nafpa.nc"

	#  Calculate bias correction for different variables separately:
    obs_var="Rainf_f_tavg"
    fcst_var="PRECTOT"
    unit="kg/m^2/s"
    var_type="PRCP"

    # Path for where nmme forecast files are located:
    fcst_clim_indir=f"{forcedir}/NMME/raw/Climatology/{month_abbr}01"
    fcst_indir=f"{forcedir}/NMME/raw/Monthly/{month_abbr}01"

    # Path for where output BC forecast file are located:
    outdir=f"{forcedir}/NMME/bcsd/Monthly/{month_abbr}01"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print(f"[INFO] Processing forecast bias correction of NMME-{nmme_model} precip")

    ens_numf, ens_numc, ens_start, ens_end = _gather_ensemble_info(nmme_model)

    for year in range(fcst_syr, (fcst_eyr + 1)):
        cmd = "sbatch"
        cmd += f" {srcdir}/run_NMME_BCSD_calctest.scr"
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
        cmd += f" {ens_numc}"
        cmd += f" {ens_numf}"
        cmd += f" {nmme_model}"
        cmd += f" {lead_months}"
        cmd += f" {year}"
        cmd += f" {year}"
        cmd += f" {clim_syr}"
        cmd += f" {clim_eyr}"
        cmd += f" {mask_file}"
        cmd += f" {fcst_clim_indir}"
        cmd += f" {obs_clim_indir}"
        cmd += f" {fcst_indir}"
        cmd += f" {outdir}"
        cmd += f" {logdir}"
        cmd += f" {ens_start}"
        cmd += f" {ens_end}"
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)

    print(f"[INFO] Completed processing NMME bias correction for: {(month_abbr)}")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
