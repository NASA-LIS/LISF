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
# SCRIPT: forecast_task_04.py
#
# PURPOSE: Computes the bias correction for the forecast (CFSv2) dataset. Based
# on FORECAST_TASK_04.sh.
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

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])}  -s fcst_syr -e fcst_eyr -m month_abbr -n month_num \
                           -w cwd - c config_file -j job_name -t ntasks -H hours"
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
    ens_num = config['BCSD']['nof_raw_ens']
    clim_syr = config['BCSD']['clim_start_year']
    clim_eyr = config['BCSD']['clim_end_year']

    # Path for where observational & forecast files are located:
    forcedir = f"{projdir}/bcsd_fcst"
    obs_indir = f"{forcedir}/USAF-LIS7.3rc8_25km"
    fcst_indir = f"{forcedir}/CFSv2_25km"

    #  Calculate bias correction for different variables separately:
    #obs_var_list = ["Rainf_f_tavg", "LWdown_f_tavg", "SWdown_f_tavg", \
    #    "Psurf_f_tavg", "Qair_f_tavg", "Tair_f_tavg", "Wind_f_tavg"]
    obs_var_list = ["PRECTOT", "LWGAB", "SWGDN", \
        "PS", "QV2M", "T2M", "U10M"]
    fcst_var_list = ["PRECTOT", "LWS", "SLRSF", "PS", "Q2M", "T2M", "WIND10M"]
    unit_list = ["kg/m^2/s", "W/m^2", "W/m^2", "Pa", "kg/kg", "K", "m/s"]

    # BC output directory for FCST:
    outdir = f"{fcst_indir}/bcsd/Monthly/{month_abbr}01"

    print("[INFO] Processing forecast bias correction of CFSv2 variables")

#    for var_num in range(len(obs_var_list)):
    for var_num, var_value in enumerate(obs_var_list):
        if var_num in [0, 2]:
            var_type = "PRCP"
        else:
            var_type = "TEMP"

        obs_var = var_value
        fcst_var = fcst_var_list[var_num]
        unit = unit_list[var_num]
        #print(f"{var_num} {fcst_var}")
        cmd = "python"
        cmd += f" {srcdir}/bias_correction_modulefast.py"
        cmd += f" {obs_var}"
        cmd += f" {fcst_var}"
        cmd += f" {var_type}"
        cmd += f" {unit}"
        cmd += f" {month_num}"
        cmd += f" {lead_months}"
        cmd += f" {ens_num}"
        cmd += f" {fcst_syr}"
        cmd += f" {fcst_eyr}"
        cmd += f" {clim_syr}"
        cmd += f" {clim_eyr}"
        cmd += f" {obs_indir}"
        cmd += f" {fcst_indir}"
        cmd += f" {config_file}"
        cmd += f" {outdir}"
        cmd += f" {logdir}"
        jobfile = job_name + '_' + obs_var + '_run.j'
        jobname = job_name + '_' + obs_var + '_'
        utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, in_command=cmd)

    print(f"[INFO] Completed writing forecast bias correction scripts for: {month_abbr}")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
