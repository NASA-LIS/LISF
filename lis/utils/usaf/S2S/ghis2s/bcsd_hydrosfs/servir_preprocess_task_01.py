#!/usr/bin/env python3
"""
#-------------------------------------------------------------------------
#
# SCRIPT: servir_preprocess_task_01.py
#
# PURPOSE: Processes the observational data and outputs in monthly time
# resolution.
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
import utils

#
# Functions
#
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} -c config_file -y year -m month \
            -w cwd -t ntasks -H hours -j job_name"
    print(txt)
    print('[INFO] where')
    print('[INFO] config_file: Config file that sets up environment')
    print('[INFO] year: Year (integer)')
    print('[INFO] month: Month (integer)')
    print('[INFO] cwd: current working directory')
    print('[INFO] ntasks: SLURM number of tasks')
    print('[INFO] hours: SLURM time hours')
    print('[INFO] job_name: SLURM job name')

def _driver():
    """Main driver."""

    # Parse command arguements
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-y', '--year', required=True, help='year (integer)')
    parser.add_argument('-m', '--month', required=True, help='month (integer)')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-t', '--ntasks', required=True, help='ntasks')
    parser.add_argument('-H', '--hours', required=True, help='hours')
    parser.add_argument('-j', '--job_name', required=True, help='job_name')

    # Create local variables for arguments (optional)
    args = parser.parse_args()
    config_file = args.config_file
    year = int(args.year)
    month = int(args.month)
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
    print('[INFO] Generating job script')
    jobfile = f"{job_name}_{year:04}{month:02}_run.j"
    jobname = f"{job_name}_{year:04}{month:02}_"

    cmd = 'python'
    cmd += f" {dir_src}/create_monthly_observation_data.py"
    cmd += f" {config_file}"
    cmd += f" {year}"
    cmd += f" {month}"
    utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, in_command = cmd, mem = 20, packable = True)

    print('[INFO] Completed generating monthly observation data job scripts')

#
# Main Method
#

if __name__ == '__main__':
    _driver()
