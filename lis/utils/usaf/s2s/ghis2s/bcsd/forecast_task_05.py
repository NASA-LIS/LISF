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

import os
import sys
import argparse
import yaml

#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {(sys.argv[0])}  -s fcst_syr -e fcst_eyr -m month_abbr \
                          -w cwd -n month_num -c config_file -j job_name \
                          -t ntasks -H hours -M NMME_MODEL"
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] fcst_eyr: End year of forecast")
    print("[INFO] clim_syr: Start year of the climatological period")
    print("[INFO] clim_eyr: End year of the climatological period")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] nmme_model: NMME model name")
    print("[INFO] lead_months: Number of lead months")
    print("[INFO] config_file: Config file that sets up environment")
    print("[INFO] cwd: current working directory")

def _driver():
    """Main driver."""

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--fcst_syr', required=True, help='forecast start year')
    parser.add_argument('-e', '--fcst_eyr', required=True, help='forecast end year')
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-m', '--month_abbr', required=True, help='month abbreviation')
    parser.add_argument('-n', '--month_num', required=True, help='month number')
    parser.add_argument('-j', '--job_name', required=True, help='job_name')
    parser.add_argument('-t', '--ntasks', required=True, help='ntasks')
    parser.add_argument('-H', '--hours', required=True, help='hours')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-M', '--nmme_model', required=True, help='NMME Model')

    args = parser.parse_args()
    config_file = args.config_file
    fcst_syr = args.fcst_syr
    fcst_eyr = args.fcst_eyr
    month_abbr = args.month_abbr
    month_num = args.month_num
    job_name = args.job_name
    ntasks = args.ntasks
    hours = args.hours
    cwd = args.cwd
    nmme_model = args.nmme_model

    # load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    # import local module
    sys.path.append(config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.shared import utils

    # Path of the main project directory
    projdir = cwd

    # Path of the directory where all the BC codes are kept
    srcdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/bcsd_fcst/bcsd_library/'

    # Log file output directory
    logdir = cwd + '/log_files'

    # Path of the directory where supplementary files are kept
    supplementary_dir = config['SETUP']['supplementarydir'] + '/bcsd_fcst/'

    lead_months = config['EXP']['lead_months']
    clim_syr = config['BCSD']['clim_start_year']
    clim_eyr = config['BCSD']['clim_end_year']
    datatype = config['SETUP']['DATATYPE']

    # Path for where observational files are located:
    forcedir = f"{projdir}/bcsd_fcst"
    obs_clim_indir = f"{forcedir}/USAF-LIS7.3rc8_25km/raw/Climatology"

    #  Calculate bias correction for different variables separately:
    #obs_var = "Rainf_f_tavg"
    obs_var = "PRECTOT"
    fcst_var = "PRECTOT"
    unit = "kg/m^2/s"
    var_type = "PRCP"

    # Path for where nmme forecast files are located:
    fcst_clim_indir = f"{forcedir}/NMME/raw/Climatology/{month_abbr}01"
    fcst_indir = f"{forcedir}/NMME/raw/Monthly/{month_abbr}01"

    # Path for where output BC forecast file are located:
    outdir = f"{forcedir}/NMME/bcsd/Monthly/{month_abbr}01"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print(f"[INFO] Processing forecast bias correction of NMME-{nmme_model} precip")

    ensemble_sizes = config['EXP']['ensemble_sizes'][0]
    ens_num = ensemble_sizes[nmme_model]

    for year in range(int(fcst_syr), (int(fcst_eyr) + 1)):
        cmd = "python"
        cmd += f" {srcdir}/bias_correction_nmme_modulefast.py"
        cmd += f" {obs_var}"
        cmd += f" {fcst_var}"
        cmd += f" {var_type}"
        cmd += f" {unit}"
        cmd += f" {month_num}"
        cmd += f" {nmme_model}"
        cmd += f" {lead_months}"
        cmd += f" {ens_num}"
        cmd += f" {year}"
        cmd += f" {year}"
        cmd += f" {clim_syr}"
        cmd += f" {clim_eyr}"     
        cmd += f" {fcst_clim_indir}"
        cmd += f" {obs_clim_indir}"
        cmd += f" {fcst_indir}"
        cmd += f" {config_file}"
        cmd += f" {outdir}"
        #cmd += f" {logdir}"
        jobfile = job_name + '_' + nmme_model + '_run.j'
        jobname = job_name + '_' + nmme_model + '_'
        utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, in_command=cmd)

    print(f"[INFO] Completed writing NMME bias correction scripts for: {(month_abbr)}")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
