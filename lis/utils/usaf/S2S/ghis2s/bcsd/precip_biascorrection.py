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
# SCRIPT: precip_biascorrection.py
#
# PURPOSE: Computes the bias correction for the NMME dataset. Based on
# FORECAST_TASK_03.sh.
#
# REVISION HISTORY:
# 24 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""

#
# Standard modules
#

import argparse
import yaml
from ghis2s.shared import utils
from ghis2s.bcsd.bcsd_library.nmme_module import NMMEParams
#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {(sys.argv[0])}  -s fcst_syr -m month_abbr \
                          -w cwd -n month_num -c config_file -j job_name \
                          -t ntasks -H hours -M NMME_MODEL"
    print(txt)
    print("[INFO] where")
    print("[INFO] fcst_syr: Start year of forecast")
    print("[INFO] clim_syr: Start year of the climatological period")
    print("[INFO] clim_eyr: End year of the climatological period")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] nmme_model: NMME model name")
    print("[INFO] lead_months: Number of lead months")
    print("[INFO] config_file: Config file that sets up environment")
    print("[INFO] cwd: current working directory")

def main(config_file, fcst_syr, month_abbr, month_num, job_name,
         ntasks, hours, cwd, nmme_model, py_call=False):
    """Main driver."""
    # load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    lats, _ = utils.get_domain_info(config_file, coord=True)
    resol = f'{round((lats[1] - lats[0])*100)}km'

    # Path of the main project directory
    projdir = cwd

    # Path of the directory where all the BC codes are kept
    srcdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/S2S/ghis2s/bcsd/bcsd_library/'

    lead_months = config['EXP']['lead_months']
    clim_syr = config['BCSD']['clim_start_year']
    clim_eyr = config['BCSD']['clim_end_year']

    # Path for where observational files are located:
    forcedir = f"{projdir}/bcsd_fcst"
    obs_clim_indir = f"{forcedir}/USAF-LIS7.3rc8_{resol}/raw/Climatology"

    #  Calculate bias correction for different variables separately:
    obs_var = "PRECTOT"
    fcst_var = "PRECTOT"
    unit = "kg/m^2/s"
    var_type = "PRCP"

    # Path for where nmme forecast files are located:
    fcst_indir = f"{forcedir}/NMME"

    print(f"[INFO] Processing forecast bias correction of NMME-{nmme_model} precip")

    ens_num = NMMEParams(nmme_model).ens_num
    slurm_commands = []
    cmd = "python"
    cmd += f" {srcdir}/bias_correction_module.py"
    cmd += f" {obs_var}"
    cmd += f" {fcst_var}"
    cmd += f" {var_type}"
    cmd += f" {unit}"
    cmd += f" {month_num}"
    cmd += f" {fcst_syr}"
    cmd += f" {nmme_model}"
    cmd += f" {ens_num}"
    cmd += f" {clim_syr}"
    cmd += f" {clim_eyr}"
    cmd += f" {config_file}"
    cmd += f" {obs_clim_indir}"
    cmd += f" {fcst_indir}"

    jobfile = job_name + '_' + nmme_model + '_run.j'
    jobname = job_name + '_' + nmme_model + '_'

    if py_call:
        for lead_num in range(lead_months):
            slurm_commands.append(cmd+f" {lead_num}")
    else:
        utils.job_script(config_file, jobfile, jobname, ntasks, hours,
                         cwd, None, in_command=cmd)

    print(f"[INFO] Completed writing NMME bias correction scripts for: {(month_abbr)}")
    if py_call:
        return slurm_commands
#
# Main Method
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--fcst_syr', required=True, help='forecast start year')
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-m', '--month_abbr', required=True, help='month abbreviation')
    parser.add_argument('-n', '--month_num', required=True, help='month number')
    parser.add_argument('-j', '--job_name', required=True, help='job_name')
    parser.add_argument('-t', '--ntasks', required=True, help='ntasks')
    parser.add_argument('-H', '--hours', required=True, help='hours')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-M', '--nmme_model', required=True, help='NMME Model')

    args = parser.parse_args()

    main(args.config_file, args.fcst_syr, args.month_abbr, args.month_num,
         args.job_name, args.ntasks, args.hours, args.cwd, args.nmme_model)
