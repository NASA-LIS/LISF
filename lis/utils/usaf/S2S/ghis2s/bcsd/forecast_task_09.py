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
# SCRIPT: forecast_task_09.py
#
# PURPOSE: Combine all non-precip 6-hourly files into one file and copy BCSD
# precip files in to the same directory. Based on FORECAST_TASK_09.sh.
#
# REVISION HISTORY:
# 24 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""

#
# Standard modules
#

import sys
import argparse
import yaml

#
# Local methods
#

def usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} -s fcst_syr -e fcst_eyr -m month_abbr \
            -w cwd -n month_num -c config_file -j job_name -t ntasks -H hours -M fcst_type"
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] fcst_eyr: End year of forecast")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] config_file: Config file that sets up environment")
    print("[INFO] fcst_type: CFSv2")
    print("[INFO] job_name: SLURM job_name")
    print("[INFO] ntasks: SLURM ntasks")
    print("[INFO] hours: SLURM time hours")
    print("[INFO] cwd: current working directory")

def driver():
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
    parser.add_argument('-M', '--fcst_type', required=True, help='NMME Model')
    parser.add_argument('-p', '--project_directory', required=True, help='Project (E2ES) directory')

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
    fcst_type = args.fcst_type

    # load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    # import local module
    sys.path.append(config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.shared import utils

    lead_months = config['EXP']['lead_months']
    ens_num = config['BCSD']['nof_raw_ens']

    # Path of the main project directory
    projdir = args.project_directory

    # Path of the directory where all the BC codes are kept:
    srcdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/bcsd_fcst/bcsd_library/'
    srcdir2 = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/bcsd_fcst/'

    # Path for the final 6-hourly forcing data:
    forcedir = f"{projdir}/bcsd_fcst/CFSv2_25km"

    print("[INFO] Combining subdaily BC CFSv2 non-precip variables")
    for year in range(int(fcst_syr), (int(fcst_eyr) + 1)):
        cmd = "python"
        cmd += f" {srcdir}/combine_sub_daily_downscaled_forcings.py"
        cmd += f" {year}"
        cmd += f" {month_num}"
        cmd += f" {fcst_type}"
        cmd += f" {ens_num}"
        cmd += f" {lead_months}"
        cmd += f" {forcedir}"
        jobfile = job_name + '_run.j'
        jobname = job_name + '_'
        utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, in_command=cmd)

    print(f"[INFO] Wrote  CFSv2 combination script for: {month_abbr}")

    # Now write task 10 scripts

    for nmme_model in  config['EXP']['NMME_models']:
        cmd = "python"
        cmd += f" {srcdir2}/forecast_task_10.py"
        cmd += f" -c {config_file}"
        cmd += f" -s {year}"
        cmd += f" -m {month_abbr}"
        cmd += f" -w {projdir}"
        cmd += f" -n {month_num}"
        cmd += f" -M {nmme_model}"
        jobfile = 'bcsd10_' + nmme_model + '_run.j'
        jobname = 'bcsd10_' + nmme_model + '_'
        utils.job_script(config_file, jobfile, jobname, ntasks, '1', cwd, in_command=cmd)

    # Now write task 11 scripts

    cmd = "python"
    cmd += f" {srcdir2}/forecast_task_11.py"
    cmd += f" -s {year}"
    cmd += f" -m {month_abbr}"
    cmd += f" -c {config_file}"
    cmd += f" -w {cwd}"
    cmd += f" -n {month_num}"
    jobfile = 'bcsd11_run.j'
    jobname = 'bcsd11_'
    utils.job_script(config_file, jobfile, jobname, ntasks, '1', cwd, in_command=cmd)

    # Now write task 12 scripts

    cmd = "python"
    cmd += f" {srcdir2}/forecast_task_12.py"
    cmd += f" -s {year}"
    cmd += f" -m {month_abbr}"
    cmd += f" -c {config_file}"
    cmd += f" -w {cwd}"
    cmd += f" -n {month_num}"
    jobfile = 'bcsd12_run.j'
    jobname = 'bcsd12_'
    utils.job_script(config_file, jobfile, jobname, ntasks, '3', cwd, in_command=cmd)

#
# Main Method
#
if __name__ == "__main__":
    driver()
