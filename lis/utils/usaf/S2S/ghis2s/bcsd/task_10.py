#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

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

import os
import subprocess
import sys
import argparse
import yaml
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {(sys.argv[0])} -s current_year -m month_abbr \
                       -w cwd -n month_num -M nmme_model"
    print(txt)
    print("[INFO] current_year: Start year of forecast")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] nmme_model: NMME_MODEL")
    print("[INFO] cwd: current working directory")

def copy_subdaily_precipitation(year, month_abbr, ens_num, indir_nmme, outdir):
    """Copies the BC 6-Hourly precipition files to the final directory."""

    for iens in range(1, (ens_num + 1)):
        print(f"[INFO] Copying precipitation files for: {year}, ens {iens}/{ens_num}")

        indir_complete = f"{indir_nmme}/{year}/ens{iens}"

        outdir_complete = f"{outdir}/{month_abbr}01/{year}/ens{iens}"

        if not os.path.exists(outdir_complete):
            os.makedirs(outdir_complete)

        cmd = f"cp {indir_complete}/PRECTOT* {outdir_complete}/"
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling copy subroutine!")
            sys.exit(1)

def driver():
    """Main driver."""

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-s', '--current_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--month_abbr', required=True, help='month abbreviation')
    parser.add_argument('-n', '--month_num', required=True, help='month number')
    parser.add_argument('-M', '--nmme_model', required=True, help='NMME Model')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')

    args = parser.parse_args()
    config_file = args.config_file
    current_year = int(args.current_year)
    month_abbr = args.month_abbr
    month_num = int(args.month_num)
    cwd = args.cwd
    nmme_model = args.nmme_model

    # load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)
   
    nmme_models = config['EXP']['NMME_models']
    ensemble_sizes = config['EXP']['ensemble_sizes'][0]
    nof_raw_ens = config['BCSD']['nof_raw_ens']

    # Path of the main project directory
    projdir = cwd

    # Path for the final 6-hourly forcing dataets:
    forcedir_fcst = f"{projdir}/bcsd_fcst/CFSv2_25km"
    forcedir_nmme = f"{projdir}/bcsd_fcst/NMME"

    # Base model prefixes for forecast files
    basemodname_src = "CFSv2"
    basemodname_dst = "CFSv2"

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

    ens_num = ensemble_sizes[nmme_model]
    nreps = int(np.ceil(ens_num/nof_raw_ens))
    ens_range = np.tile(list(range(1,nof_raw_ens+1)), nreps)[0:ens_num]
    
    indir_nmme = f"{forcedir_nmme}/bcsd/6-Hourly/{month_abbr}01/{nmme_model}"
    outdir = f"{forcedir_nmme}/final/6-Hourly/{nmme_model}"

    # Copy the precipitation files
    print(f"[INFO] NMME MODEL: {nmme_model}")
    copy_subdaily_precipitation(current_year, month_abbr, ens_num, indir_nmme,\
         outdir)

    # Symbolically link the non-precip data
    print("[INFO] Creating symbolic links for non-precip data")
    for iens, ens_value in enumerate(ens_range):
        ens_nmme = iens + 1
        ens_fcst = ens_value

        src_dir = f"{forcedir_fcst}/final/6-Hourly/{month_abbr}01/{current_year}/ens{ens_fcst}"
        dst_dir = f"{outdir}/{month_abbr}01/{current_year}/ens{ens_nmme}"

        for ilead, src_lead in enumerate(src_yyyymm):
            src_file = f"{src_dir}/{basemodname_src}.{src_lead}.nc4"
            dst_file = f"{dst_dir}/{basemodname_dst}.{dst_yyyymm[ilead]}.nc4"

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
