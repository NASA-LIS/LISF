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
# SCRIPT: forecast_task_07.py
#
# PURPOSE: Combine all non-precip 6-hourly files into one file and copy BCSD
# precip files in to the same directory Based on FORECAST_TASK_07.sh.
#
# REVISION HISTORY:
# 24 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""

#
# Standard modules
#

import subprocess
import os
import sys
import argparse

#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} -s current_year -m month_abbr -w cwd"
    print(txt)
    print("[INFO] where")
    print("[INFO] current_year: Current year")
    print("[INFO] month_abbr: Current month")
    print("[INFO] cwd: current working directory")

def _driver():
    """Main driver."""

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--current_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--month_abbr', required=True, help='month abbreviation')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')

    args = parser.parse_args()
    current_year = args.current_year
    month_abbr = args.month_abbr
    cwd = args.cwd

    # Path of the main project directory
    projdir = cwd

    # Number of precip ensembles needed
    range_ens_fcst = list(range(1, 13)) + list(range(1, 13)) + list(range(1, 7))
    range_ens_nmme = range(1, 31)

    fcst_date = f"{month_abbr}01"

    # Path for where forecast files are located:
    indir = f"{projdir}/bcsd_fcst/CFSv2_25km/raw/6-Hourly/{fcst_date}/{current_year}"

    # Path for where the linked precip files should be placed:
    outdir = f"{projdir}/bcsd_fcst/NMME/linked_cfsv2_precip_files/{fcst_date}/{current_year}"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for iens, ens_value in enumerate(range_ens_fcst):
        src_file = f"{indir}/ens{ens_value}"
        dst_file = f"{outdir}/ens{range_ens_nmme[iens]}"

        cmd = f"ln -sfn {src_file} {dst_file}"
        print(cmd)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling creating symbolic links!")
            sys.exit(1)

    print("[INFO] Done creating symbolic links")

#
# Main Method
#
if __name__ == "__main__":
    _driver()
