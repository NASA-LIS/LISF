#!/usr/bin/env python3
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

import configparser
import os
import subprocess
import sys

#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} current_year month_abbr config_file"
    print(txt)
    print("[INFO] where")
    print("[INFO] current_year: Current year")
    print("[INFO] month_abbr: Current month")
    print("[INFO] config_file: Config file that sets up environment")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # current_year
    try:
        current_year = int(sys.argv[1])
    except ValueError:
        print(f"[ERR] Invalid argument for current_year!  Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)
    if current_year < 0:
        print(f"[ERR] Invalid argument for current_year!  Received {(sys.argv[1])}")
        _usage()
        sys.exit(1)

    # month_abbr
    month_abbr = str(sys.argv[2])

    # config_file
    config_file = sys.argv[3]
    if not os.path.exists(config_file):
        print(f"[ERR] {config_file} does not exist!")
        sys.exit(1)

    return current_year, month_abbr, config_file

def read_config(config_file):
    """Read from bcsd_preproc config file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def _driver():
    """Main driver."""
    current_year, month_abbr, config_file = _read_cmd_args()

    # Setup local directories
    config = read_config(config_file)

    # Path of the main project directory
    projdir = config["bcsd_preproc"]["projdir"]

	# Number of precip ensembles needed
    range_ens_fcst=list(range(1, 13)) + list(range(1,13)) + list(range(1,7))
    range_ens_nmme=range(1,31)

    fcst_date = f"{month_abbr}01"

    # Path for where forecast files are located:
    indir=f"{projdir}/data/forecast/CFSv2_25km/raw/6-Hourly/{fcst_date}/{current_year}"

    # Path for where the linked precip files should be placed:
    outdir=f"{projdir}/data/forecast/NMME/linked_cfsv2_precip_files/{fcst_date}/{current_year}"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for iens, ens_value in enumerate(range_ens_fcst):
        src_file=f"{indir}/ens{ens_value}"
        dst_file=f"{outdir}/ens{range_ens_nmme[iens]}"

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
