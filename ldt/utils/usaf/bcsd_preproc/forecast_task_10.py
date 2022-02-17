#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_10.py
#
# PURPOSE: Copy BC NMME precip files in to the same directory. Based on
# FORECAST_TASK_10.sh.
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
from datetime import datetime
from dateutil.relativedelta import relativedelta

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {(sys.argv[0])} MONTH_ABBR MONTH_NUM CURRENT_YEAR NMME_MODEL CONFIG_FILE"
    print(txt)
    print("[INFO] where")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] MONTH_NUM: Integer number of the initialization month")
    print("[INFO] CURRENT_YEAR: Current year of forecast")
    print("[INFO] NMME_MODEL: NMME model name")
    print("[INFO] CONFIG_FILE: Config file that sets up environment")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 6:
        print("[ERR] Invalid number of command line arguments!")
        usage()
        sys.exit(1)

    # MONTH_ABBR
    month_abbr = str(sys.argv[1])

    # MONTH_NUM
    try:
        month_num = int(sys.argv[2])
    except ValueError:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)
    if month_num < 1:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)
    if month_num > 12:
        print(f"[ERR] Invalid argument for MONTH_NUM! Received {(sys.argv[2])}")
        usage()
        sys.exit(1)

    # CURRENT_YEAR
    try:
        current_year = int(sys.argv[3])
    except ValueError:
        print(f"[ERR] Invalid argument for CURRENT_YEAR! Received {(sys.argv[3])}")
        usage()
        sys.exit(1)
    if current_year < 0:
        print(f"[ERR] Invalid argument for CURRENT_YEAR! Received {(sys.argv[3])}")
        usage()
        sys.exit(1)

    # NMME_MODEL
    nmme_model = str(sys.argv[4])

    # CONFIG_FILE
    config_file = sys.argv[5]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return month_abbr, month_num, current_year, nmme_model, config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def gather_ensemble_info(nmme_model):
    """Gathers ensemble information based on NMME model."""

    # Number of ensembles in the forecast (ENS_NUM)
    # Ensemble start index (ENS_START)
    # Ensemble end index (ENS_END)
    if nmme_model == "CFSv2":
        ens_num=24
        ens_range=list(range(1, 13)) + list(range(1,13))
    elif nmme_model == "GEOSv2":
        ens_num=10
        ens_range=range(1, 11)
    elif nmme_model == "CCM4":
        ens_num=10
        ens_range=range(1, 11)
    elif nmme_model == "GNEMO":
        ens_num=10
        ens_range=range(1, 11)
    elif nmme_model == "CCSM4":
        ens_num=10
        ens_range=range(1, 11)
    elif nmme_model == "GFDL":
        ens_num=30
        ens_range=list(range(1, 13)) + list(range(1,13)) + list(range(1,7))
    else:
        print(f"[ERR] Invalid argument for NMME_MODEL! Received {nmme_model}")
        sys.exit(1)

    return ens_num, ens_range

def copy_subdaily_precipitation(year, month_abbr, ens_num, indir_nmme, outdir):
    """Copies the BC 6-Hourly precipition files to the final directory."""

    for iens in range(1, (ens_num + 1)):
        print(f"[INFO] Copying precipitation files for: {year}, ens {iens}/{ens_num}")

        indir_complete=f"{indir_nmme}/{year}/ens{iens}"

        outdir_complete=f"{outdir}/{year}/{month_abbr}01/ens{iens}"

        if not os.path.exists(outdir_complete):
            os.makedirs(outdir_complete)

        cmd = f"cp {indir_complete}/PRECTOT* {outdir_complete}/"
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling copy subroutine!")
            sys.exit(1)

def driver():
    """Main driver."""
    month_abbr, month_num, current_year, nmme_model, config_file = \
    	read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

    # Path for the final 6-hourly forcing dataets:
    forcedir_fcst = f"{projdir}/data/forecast/CFSv2_25km"
    forcedir_nmme = f"{projdir}/data/forecast/NMME"

    # Base model prefixes for forecast files
    basemodname_src="CFSv2"
    basemodname_dst="CFSv2"

    init_datetime = datetime(current_year, month_num, 1)
    src_yyyymm = [(init_datetime + relativedelta(months=0)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=1)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=2)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=3)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=4)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=5)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=6)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=7)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=8)).strftime("%Y%m"),
                  (init_datetime + relativedelta(months=8)).strftime("%Y%m")]

    dst_yyyymm = [(init_datetime + relativedelta(months=0)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=1)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=2)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=3)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=4)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=5)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=6)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=7)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=8)).strftime("%Y%m"),
              (init_datetime + relativedelta(months=9)).strftime("%Y%m")]

    ens_num, ens_range = gather_ensemble_info(nmme_model)

    indir_nmme=f"{forcedir_nmme}/bcsd/6-Hourly/{month_abbr}01/{nmme_model}"
    outdir=f"{forcedir_nmme}/final/6-Hourly/{nmme_model}"

    # Copy the precipitation files
    print(f"[INFO] NMME MODEL: {nmme_model}")
    copy_subdaily_precipitation(current_year, month_abbr, ens_num, indir_nmme,\
         outdir)

    # Symbolically link the non-precip data
    print("[INFO] Creating symbolic links for non-precip data")
#    for iens in range(len(ens_range)):
    for iens,ens_value in enumerate(ens_range):
        ens_nmme=iens + 1
        ens_fcst=ens_value

        src_dir=f"{forcedir_fcst}/final/6-Hourly/{current_year}/{month_abbr}01/ens{ens_fcst}"
        dst_dir=f"{outdir}/{current_year}/{month_abbr}01/ens{ens_nmme}"

#        for ilead in range(len(src_yyyymm)):
        for ilead,src_lead in enumerate(src_yyyymm):
            src_file=f"{src_dir}/{basemodname_src}.{src_lead}.nc4"
            dst_file=f"{dst_dir}/{basemodname_dst}.{dst_yyyymm[ilead]}.nc4"

            cmd = f"ln -sfn {src_file} {dst_file}"
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem calling creating symbolic links!")
                sys.exit(1)

#
# Main Method
#
if __name__ == "__main__":
    driver()
