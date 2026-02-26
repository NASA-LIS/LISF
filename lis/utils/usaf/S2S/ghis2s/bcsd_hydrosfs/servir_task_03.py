#!/usr/bin/env python3
""" Wrapper function to run the first task of Bias Correction in the SERVIR-S2S project

Computes the temporal disaggregation for the forecast dataset at 6-hourly time resolution

Example:
    To run this script from the command line:

    $ python servir_task_02.py -c $BWD/$CFILE -y $YYYY -m $mon -w ${CWD} -t 1 -H 9 -j ${jobname}
    
Author: Ryan Zamora
"""

#
# Modules
#
import sys
import argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
import yaml

#
# Custom Modules
#
import utils
from dict_variables import (
    get_all_key_variable_names,
    get_hydrosfs_name,
    get_disaggregation_type,
    get_units
)

#
# Functions
#
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} -c config_file -y fcst_init_year -m fcst_init_month \
                          -w cwd -t ntasks -H hours -j job_name"
    print(txt)
    print("[INFO] where")
    print("[INFO] config_file: Config file that sets up environment")
    print('[INFO] fcst_init_year: Initial forecast year (integer)')
    print('[INFO] fcst_init_month: Initial forecast month (integer)')
    print('[INFO] cwd: current working directory')
    print('[INFO] ntasks: SLURM number of tasks')
    print('[INFO] hours: SLURM time hours')
    print('[INFO] job_name: SLURM job name')

def _driver():
    """Main driver."""

    # Setup local directories
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-y', '--fcst_init_year', required=True, help='initial forecast year')
    parser.add_argument('-m', '--fcst_init_month', required=True, help='initial forecast month')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-t', '--ntasks', required=True, help='ntasks')
    parser.add_argument('-H', '--hours', required=True, help='hours')
    parser.add_argument('-j', '--job_name', required=True, help='job_name')

    # Create local variables for arguments (optional)
    args = parser.parse_args()
    config_file = args.config_file
    fcst_init_year = int(args.fcst_init_year)
    fcst_init_month = int(args.fcst_init_month)
    cwd = args.cwd
    ntasks = args.ntasks
    hours = args.hours
    job_name = args.job_name

    # load config file
    with open(config_file, 'r', encoding = 'utf-8') as file:
        config = yaml.safe_load(file)

    # Path of the directory where all of the bias-correction code is kept
    dir_src = config['SETUP']['DIR_SCRIPTS'] + '/bcsd/'

    # Number of forecast ensemble members
    nof_ensemble_members = config['BCSD']['nof_raw_ens']

    # Number of forecast lead months
    nof_lead_months = config['EXP']['lead_months']

    print('[INFO] Generating a separate job script for each variable')
    for ensemble_number in range(1, nof_ensemble_members + 1):
        for var_name in get_all_key_variable_names():
            fcst_var_name = get_hydrosfs_name(var_name)
            var_type = get_disaggregation_type(var_name)
            var_unit = get_units(var_name)
            
            jobname = f"{job_name}_ens{ensemble_number:02}_{fcst_var_name}_"
            jobfile = f"{jobname}run.j"

            cmd_list = []
            for lead_month in range(nof_lead_months):
    
                # Calculate lead days for given lead month
                init_date = datetime(fcst_init_year, fcst_init_month, 1)
                start_date = datetime(fcst_init_year, fcst_init_month, 1) + relativedelta(months = lead_month)
                end_date = start_date + relativedelta(months = 1)
    
                start_lead_day = (start_date - init_date).days
                end_lead_day   = (end_date - init_date).days
    
                for lead_day in range(start_lead_day, end_lead_day):
                    cmd = 'python'
                    cmd += f" {dir_src}/create_temporal_disaggregated_6hourly_data.py"
                    cmd += f" {config_file}"
                    cmd += f" {fcst_init_year}"
                    cmd += f" {fcst_init_month}"
                    cmd += f" {ensemble_number}"
                    cmd += f" {lead_month}"
                    cmd += f" {lead_day}"
                    cmd += f" {fcst_var_name}"
                    cmd += f" {var_type}"
                    cmd += f" {var_unit}"
                    cmd_list.append(cmd)

            # Write job script
            utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, command_list = cmd_list, mem = 20, packable = True)

    print('[INFO] Completed generating forecast temporal disaggregation job scripts')

#
# Main
#
if __name__ == "__main__":
    _driver()
