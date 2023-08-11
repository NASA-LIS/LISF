#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

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
import os
import sys
import argparse
import yaml
#pylint: disable=import-outside-toplevel, too-many-locals
# Local methods
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} -s fcst_syr -e fcst_eyr -m month_abbr \
            -w cwd - c config_file -j job_name -t ntasks -H hours "
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] fcst_eyr: End year of forecast")
    print("[INFO] month_abbr: Abbreviated month to start forecast")
    print("[INFO] cwd: current working directory")
    print("[INFO] config_file: Config file that sets up environment")
    print("[INFO] job_name: SLURM job_name")
    print("[INFO] ntasks: SLURM ntasks")
    print("[INFO] hours: SLURM time hours")

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

    # Parse command arguements
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--fcst_syr', required=True, help='forecast start year')
    parser.add_argument('-e', '--fcst_eyr', required=False, help='forecast end year')
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-m', '--month_abbr', required=True, help='month abbreviation')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-j', '--job_name', required=True, help='job_name')
    parser.add_argument('-t', '--ntasks', required=True, help='ntasks')
    parser.add_argument('-H', '--hours', required=True, help='hours')

    args = parser.parse_args()
    config_file = args.config_file
    syear = int(args.fcst_syr)
    if args.fcst_eyr is None:
        eyear = None
    else:
        eyear = int(args.fcst_eyr)
    month_abbr = args.month_abbr
    cwd = args.cwd
    job_name = args.job_name
    ntasks = args.ntasks
    hours = args.hours

    # Load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    # Import local module
    sys.path.append(config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.shared import utils

    # Path of the main project directory
    projdir = cwd

    # Path of the directory where all the BC codes are kept
    srcdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/bcsd_fcst/bcsd_library/'
    # Log file output directory
    logdir = cwd + '/log_files'

    # Paths for the daily forecast data (input and output paths)
    outdir = f"{projdir}/bcsd_fcst/CFSv2_25km/raw"

    if not os.path.exists(logdir):
        os.makedirs(logdir)

    imon = f"{month_abbr}01"
    ic_dates = calc_ic_dates(imon)

    # Process 6-hrly CFSv2 forecasts and output in monthly and 6-hrly formats
    print("[INFO] Processing CFSv2 6-hrly forecast variables")
    nof_raw_ens = config['BCSD']['nof_raw_ens']
    for ens_num in range(1, nof_raw_ens + 1):
        if eyear is not None:
            cmd_list = []
            for cyear in range(syear,eyear+1):
                cmd = "python"
                cmd += f" {srcdir}/process_forecast_data.py"
                cmd += f" {cyear:04d}"
                cmd += f" {ens_num:02d}"
                cmd += f" {imon}"
                cmd += f" {outdir}"
                cmd += f" {config_file}"
                for ic_date in ic_dates:
                    cmd += f" {ic_date}"
                cmd_list.append(cmd)
            jobfile = job_name + '_' + str(ens_num).zfill(2) + '_run.j'
            jobname = job_name + '_' + str(ens_num).zfill(2) + '_'
            utils.job_script(config_file, jobfile, jobname, ntasks,
                             hours, cwd, command_list=cmd_list)

        else:
            cmd = "python"
            cmd += f" {srcdir}/process_forecast_data.py"
            cmd += f" {syear:04d}"
            cmd += f" {ens_num:02d}"
            cmd += f" {imon}"
            cmd += f" {outdir}"
            cmd += f" {config_file}"
            for ic_date in ic_dates:
                cmd += f" {ic_date}"
            jobfile = job_name + '_' + str(ens_num).zfill(2) + '_run.j'
            jobname = job_name + '_' + str(ens_num).zfill(2) + '_'
            utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, in_command=cmd)

    print(f"[INFO] Write command to process CFSv2 files for {imon}")

if __name__ == "__main__":
    _driver()
