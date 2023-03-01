#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_08.py
#
# PURPOSE: Generate bias-corrected 6-hourly nmme forecasts using raw monthly
# forecasts, bias-corrected monthly forecasts and raw 6-hourly forecasts. Based
# on FORECAST_TASK_08.sh.
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
    txt = f"[INFO] Usage: {(sys.argv[0])} -s fcst_syr -e fcst_eyr -m month_abbr \
             -w cwd -n month_num -c config_file -j job_name -t ntasks -H hours -M NMME_MODEL"
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] fcst_eyr: End year of forecast")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] nmme_model: NMME model name")
    print("[INFO] config_file: Config file that sets up environment")
    print("[INFO] job_name: SLURM job_name")
    print("[INFO] ntasks: SLURM ntasks")
    print("[INFO] hours: SLURM time hours")

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
    #logdir = cwd + '/log_files'

    # Path of the directory where supplementary files are kept
    supplementary_dir = config['SETUP']['supplementarydir'] + '/bcsd_fcst/'

    # domain
    domain = config['EXP']['DOMAIN']

    lead_months = config['EXP']['lead_months']
    datatype = config['SETUP']['DATATYPE']
    ensemble_sizes = config['EXP']['ensemble_sizes'][0]
    ens_num = ensemble_sizes[nmme_model]

    # Path for where forecast files are located:
    forcedir = f"{projdir}/bcsd_fcst/NMME"

    #  Calculate bias correction for different variables separately:
    obs_var = "PRECTOT"
    fcst_var = "PRECTOT"
    unit = "kg/m^2/s"
    var_type = 'PRCP'

    # Path for where forecast and bias corrected files are located:
    subdaily_raw_fcst_dir = f"{forcedir}/linked_cfsv2_precip_files/{month_abbr}01"
    monthly_raw_fcst_dir = f"{forcedir}/raw/Monthly/{month_abbr}01"
    monthly_bc_fcst_dir = f"{forcedir}/bcsd/Monthly/{month_abbr}01"

    outdir = f"{forcedir}/bcsd/6-Hourly/{month_abbr}01/{nmme_model}"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("[INFO] Processing temporal disaggregation of CFSv2 variables")
    for year in range(int(fcst_syr), (int(fcst_eyr) + 1)):
        cmd = "python"
        cmd += f" {srcdir}/temporal_disaggregation_nmme_6hourly_module.py"
        cmd += f" {obs_var}"
        cmd += f" {fcst_var}"
        cmd += f" {year}"
        cmd += f" {month_num}"
        cmd += f" {var_type}"
        cmd += f" {unit}"
        cmd += f" {nmme_model}"
        cmd += f" {ens_num}"
        cmd += f" {lead_months}"
        cmd += f" {year}"
        cmd += f" {year}"
        cmd += f" {config_file}"
        cmd += f" {monthly_bc_fcst_dir}"
        cmd += f" {subdaily_raw_fcst_dir}"
        cmd += f" {outdir}"
        cmd += f" {domain}"
        jobfile = job_name + '_' + nmme_model + '_run.j'
        jobname = job_name + '_' + nmme_model + '_'
        utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, in_command=cmd)

    print(f"[INFO] Wrote NMME temporal disaggregation script for: {month_abbr}")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
