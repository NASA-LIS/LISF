#!/usr/bin/env python3
""" Wrapper function to run the first task of Bias Correction in the SERVIR-S2S project

Processes the forecast data and outputs in 6-hourly and monthly time resolutions.

Example:
    To run this script from the command line:

    $ python servir_task_01.py -c $BWD/$CFILE -y $YYYY -m $mmm -w ${CWD} -t 1 -H 7 -j ${jobname}_${YYYY}
    
Author: Ryan Zamora
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
    txt = f"[INFO] Usage: {(sys.argv[0])} -c config_file -y fcst_init_year -m fcst_init_month \
            -w cwd -t ntasks -H hours -j job_name"
    print(txt)
    print('[INFO] where')
    print('[INFO] config_file: Config file that sets up environment')
    print('[INFO] fcst_init_year: Initial forecast year (integer)')
    print('[INFO] fcst_init_month: Initial forecast month (integer)')
    print('[INFO] cwd: current working directory')
    print('[INFO] ntasks: SLURM number of tasks')
    print('[INFO] hours: SLURM time hours')
    print('[INFO] job_name: SLURM job name')

def _driver():
    """Main driver."""

    # Parse command arguements
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

    # Load the S2S config file
    with open(config_file, 'r', encoding = 'utf-8') as file:
        config = yaml.safe_load(file)

    # Path of the directory where all of the bias-correction code is kept
    dir_src = config['SETUP']['DIR_SCRIPTS'] + '/bcsd/'

    # Number of forecast ensemble members
    nof_ensemble_members = config['BCSD']['nof_raw_ens']

    # Number of forecast lead months
    nof_lead_months = config['EXP']['lead_months']

    # Determine rescale_reorg function to run
    if config['BCSD']['fcst_data_type'] == 'CFSv2':
        rescale_reorg_function = 'rescale_reorg_cfsv2_data.py'
    elif config['BCSD']['fcst_data_type'] == 'GEOSv3':
        rescale_reorg_function = 'rescale_reorg_geosv3_data.py'
    else:
        raise ValueError('forecast_data_type in config must either be: CFSv2 or GEOSv3')

    print('[INFO] Generating a separate job script for each ensemble member and lead month')
    for ensemble_number in range(1, nof_ensemble_members + 1):
        for lead_month in range(nof_lead_months):
            jobname = f"{job_name}_ens{ensemble_number:02}_lead_{lead_month}_"
            jobfile = f"{jobname}run.j"

            # Generate command list to submit multiple python commands
            cmd_list = []

            # Method: rescale and reorganize forecast data
            cmd = 'python'
            cmd += f" {dir_src}/{rescale_reorg_function}"
            cmd += f" {config_file}"
            cmd += f" {fcst_init_year}"
            cmd += f" {fcst_init_month}"
            cmd += f" {ensemble_number}"
            cmd += f" {lead_month}"
            cmd_list.append(cmd)

            # Method: create monthly averaged forecast data
            cmd = 'python'
            cmd += f" {dir_src}/create_monthly_forecast_data.py"
            cmd += f" {config_file}"
            cmd += f" {fcst_init_year}"
            cmd += f" {fcst_init_month}"
            cmd += f" {ensemble_number}"
            cmd += f" {lead_month}"
            cmd_list.append(cmd)

            # Write job script
            utils.job_script(config_file, jobfile, jobname, ntasks,
                             hours, cwd, command_list = cmd_list, mem = 20, packable = True)

    print('[INFO] Completed generating forecast rescale and reorganization job scripts')

#
# Main Method
#

if __name__ == '__main__':
    _driver()
