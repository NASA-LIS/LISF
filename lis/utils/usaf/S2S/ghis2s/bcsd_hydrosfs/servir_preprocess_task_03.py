#!/usr/bin/env python3

"""
#-------------------------------------------------------------------------
#
# SCRIPT: servir_preprocess_task_03.py
#
# PURPOSE: Creates forecast climatology data
#
# REVISION HISTORY:
# 22 Jul 2024: Ryan Zamora, first version
#
#-------------------------------------------------------------------------
"""

#
# Modules
#
import sys
import argparse
import yaml

#
# Custom Modules
#
from dict_variables import get_all_hydrosfs_names
import utils

#
# Functions
#
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} -c config_file -w cwd -t ntasks -H hours -j job_name"
    print(txt)
    print('[INFO] where')
    print('[INFO] config_file: Config file that sets up environment')
    print('[INFO] cwd: current working directory')
    print('[INFO] ntasks: SLURM number of tasks')
    print('[INFO] hours: SLURM time hours')
    print('[INFO] job_name: SLURM job name')

def _driver():
    """Main driver."""

    # Parse command arguements
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-m', '--fcst_init_month', required=True, help='initial forecast month (integer)')    
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-t', '--ntasks', required=True, help='ntasks')
    parser.add_argument('-H', '--hours', required=True, help='hours')
    parser.add_argument('-j', '--job_name', required=True, help='job_name')

    # Create local variables for arguments (optional)
    args = parser.parse_args()
    config_file = args.config_file
    fcst_init_month = args.fcst_init_month
    cwd = args.cwd
    ntasks = args.ntasks
    hours = args.hours
    job_name = args.job_name

    # Load the S2S config file
    with open(config_file, 'r', encoding = 'utf-8') as file:
        config = yaml.safe_load(file)

    # Path of the directory where all of the bias-correction code is kept
    dir_src = config['SETUP']['DIR_SCRIPTS'] + '/bcsd'

    # Process hourly observations and output in monthly format
    print('[INFO] Generating job scripts')


    for var_name in get_all_hydrosfs_names():
        jobname = f"{job_name}_{var_name}_"
        jobfile = f"{jobname}run.j"
        
        cmd = 'python'
        cmd += f" {dir_src}/create_hindcast_climatology.py"
        cmd += f" {config_file}"
        cmd += f" {fcst_init_month}"
        cmd += f" {var_name}"
    
        # Write job script
        utils.job_script(
            config_file, jobfile, jobname, ntasks, hours, cwd, in_command = cmd, mem = 480)

    print('[INFO] Completed generating monthly forecast data job scripts')

#
# Main Method
#

if __name__ == '__main__':
    _driver()
