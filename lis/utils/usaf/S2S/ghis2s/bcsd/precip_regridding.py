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
# SCRIPT: precip_regridding.py
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
from ghis2s.shared import utils

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

def main(config_file,  month_num, job_name, ntasks, hours, cwd, current_year=None, py_call=False):
    """Main driver."""    

    # load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    # Path of the main project directory
    projdir = cwd

    # Path of the directory where all the BC codes are kept
    srcdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/S2S/ghis2s/bcsd/bcsd_library/'

    # Path for where raw and bias corrected forecast files are located:
    forcedir = f"{projdir}/bcsd_fcst/NMME"
    nmme_output_dir = f"{forcedir}/raw/Monthly"
    slurm_commands = []
    for nmme_model in  config['EXP']['NMME_models']:
        cmd = "python"
        cmd += f" {srcdir}/nmme_module.py"
        cmd += f" {month_num}"
        if current_year is not None:
            cmd += f" {current_year}"
        cmd += f" {nmme_output_dir}"
        cmd += f" {nmme_model}"
        cmd += f" {config_file}"
        jobfile = job_name + '_' + nmme_model + '_run.j'
        jobname = job_name + '_' + nmme_model + '_'
        if py_call:
            slurm_commands.append(cmd)
        else:
            utils.job_script(config_file, jobfile, jobname, ntasks, hours,
                             cwd, None, in_command=cmd)

    if py_call:
        return slurm_commands

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--current_year', required=True, help='forecast start year')
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-m', '--month_num', required=True, type=int, help='month number')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-j', '--job_name', required=True, help='job_name')
    parser.add_argument('-t', '--ntasks', required=True, help='ntasks')
    parser.add_argument('-H', '--hours', required=True, help='hours')
    args = parser.parse_args()

    main(args.config_file, args.month_num, args.job_name,
         args.ntasks, args.hours, args.cwd, current_year=args.current_year)
