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
# SCRIPT: combine_forcings.py
#
# PURPOSE: Combine all non-precip 6-hourly files into one file and copy BCSD
# precip files in to the same directory. 
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
from ghis2s.shared import utils
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

def main(config_file, fcst_syr, fcst_eyr, month_abbr, month_num, job_name,
         ntasks, hours, cwd, projdir, fcst_type, py_call=False):
    """Main driver."""
    # load config file
    with open(config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    # Base forecast model
    fcst_model = config['BCSD']['fcst_data_type']
    
    # get resolution
    lats, lons = utils.get_domain_info(config_file, coord=True)
    resol = f'{round((lats[1] - lats[0])*100)}km'
        
    lead_months = config['EXP']['lead_months']
    ens_num = config['BCSD']['nof_raw_ens']

    # Path of the directory where all the BC codes are kept:
    srcdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/S2S/ghis2s/bcsd/bcsd_library/'
    srcdir2 = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/S2S/ghis2s/bcsd/'

    # Path for the final 6-hourly forcing data:
    forcedir = f"{projdir}/bcsd_fcst/{fcst_model}_{resol}"

    print("[INFO] Combining subdaily BC CFSv2 non-precip variables")
    slurm_9_10 = []
    for year in range(int(fcst_syr), (int(fcst_eyr) + 1)):
        cmd = "python"
        cmd += f" {srcdir}/combine_sub_daily_downscaled_forcings.py"
        cmd += f" {year}"
        cmd += f" {month_num}"
        cmd += f" {fcst_type}"
        cmd += f" {ens_num}"
        cmd += f" {lead_months}"
        cmd += f" {forcedir}"
        cmd += f" {config_file}"
        jobfile = job_name + '_run.j'
        jobname = job_name + '_'

        if py_call:
            slurm_9_10.append(cmd)
        else:
            utils.job_script(config_file, jobfile, jobname, ntasks, hours, cwd, in_command=cmd)

    print(f"[INFO] Wrote  CFSv2 combination script for: {month_abbr}")

    if py_call:
        return slurm_9_10
#
# Main Method
#
if __name__ == "__main__":
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
    
    main(args.config_file, args.fcst_syr, args.fcst_eyr, args.month_abbr, args.month_num, args.job_name,
         args.ntasks, args.hours, args.cwd, args.project_directory, args.fcst_type)
