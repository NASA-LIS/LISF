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
# SCRIPT: forecast_task_03.py
#
# PURPOSE: Reorganizes the downloaded NMME data into a format for further
# processing. Based on FORECAST_TASK_03.sh.
#
# REVISION HISTORY:
# 24 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""


# Standard modules
import sys
import argparse
import yaml

# Local methods
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} -m month_num -s current_year -c config_file \
                          -w cwd -j job_name -t ntasks -H hours"
    print(txt)
    print("[INFO] where")
    print("[INFO] month_num: Current integer month")
    print("[INFO] current_year: Current year")
    print("[INFO] config_file: Config file that sets up environment")
    print("[INFO] cwd: current working directory")
    print("[INFO] job_name: SLURM job_name")
    print("[INFO] ntasks: SLURM ntasks")
    print("[INFO] hours: SLURM time hours")

def _driver():
    """Main driver."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--current_year', required=True, help='forecast start year')
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-m', '--month_num', required=True, type=int, help='month number')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-j', '--job_name', required=True, help='job_name')
    parser.add_argument('-t', '--ntasks', required=True, help='ntasks')
    parser.add_argument('-H', '--hours', required=True, help='hours')

    args = parser.parse_args()
    config_file = args.config_file
    current_year = args.current_year
    month_num = args.month_num
    job_name = args.job_name
    ntasks = args.ntasks
    hours = args.hours
    cwd = args.cwd

    # load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    # import local module
    sys.path.append(config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
#pylint: disable=import-outside-toplevel
    from s2s_modules.shared import utils

    # Path of the main project directory
    projdir = cwd

    # Path of the directory where all the BC codes are kept
    srcdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/bcsd_fcst/bcsd_library/'

    # Path of the directory where supplementary files are kept
    supplementary_dir = config['SETUP']['supplementarydir'] + '/bcsd_fcst/'

    # List of NMME models and ensemble sizes to use
    ensemble_sizes = config['EXP']['ensemble_sizes'][0]

    # Path for where raw and bias corrected forecast files are located:
    nmme_download_dir = config['BCSD']['nmme_download_dir']
    forcedir = f"{projdir}/bcsd_fcst/NMME"
    nmme_output_dir = f"{forcedir}/raw/Monthly"

    for nmme_model in  config['EXP']['NMME_models']:
        ensemble_size = ensemble_sizes[nmme_model]
        cmd = "python"
        cmd += f" {srcdir}/nmme_reorg_f.py"
        cmd += f" {month_num}"
        cmd += f" {current_year}"
        cmd += f" {nmme_download_dir}"
        cmd += f" {nmme_output_dir}"
        cmd += f" {supplementary_dir}"
        cmd += f" {nmme_model}"
        cmd += f" {ensemble_size}"
        cmd += f" {config_file}"
        jobfile = job_name + '_' + nmme_model + '_run.j'
        jobname = job_name + '_' + nmme_model + '_'
        utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, in_command=cmd)

if __name__ == "__main__":
    _driver()
