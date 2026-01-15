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
# SCRIPT: customize_lis_config.py
#
# PURPOSE: Customizes lis.config file for 1-day run of LIS in AGRMET Ops mode,
# single ensemble member, over global domain.  This will extend the AGRMET
# forcing used by the LIS S2S system; based on the NAFPA forcing generation
# (Kemp et al., 2022).
#
# REQUIREMENTS as of 21 Nov 2022:
# * Python 3.9 or higher.
#
# REVISION HISTORY:
# 21 Sep 2021: Eric Kemp (SSAI), first version.
# 02 Nov 2021: Eric Kemp/SSAI, tweaks to satisfy pylint.
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import sys

def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} lis_config_template restart_dir YYYYMMDD"
    print(txt)
    print("[INFO] where: ")
    print("[INFO]  lis_config_template: Path to sample lis.config file.")
    print("[INFO]  restart_dir: Path to LIS restart files.")
    print("[INFO]  YYYYMMDD: year/month/day of start of next LIS run.")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        print(len(sys.argv))
        print(sys.argv[:])
        _usage()
        sys.exit(1)

    # Check if lis.config template exists.
    lis_config_template = sys.argv[1]
    if not os.path.exists(lis_config_template):
        print(f"[ERR] {lis_config_template} does not exist!")
        sys.exit(1)

    # Check if directory for restart files exists. Actual restart file
    # shall be checked later.
    restart_dir = sys.argv[2]
    if not os.path.exists(restart_dir):
        print(f"[ERR] Directory {restart_dir} does not exist!")
        sys.exit(1)

    # Get start date of new LIS run.
    yyyymmdd = sys.argv[3]
    if len(yyyymmdd) != 8:
        print("[ERR] Invalid length for YYYYMMDD, must be 8 characters!")
        sys.exit(1)
    year = int(yyyymmdd[0:4])
    month = int(yyyymmdd[4:6])
    day = int(yyyymmdd[6:8])
    try:
        startdate = datetime.date(year, month, day)
    except ValueError:
        print("[ERR] Invalid YYYYMMDD passed to script!")
        sys.exit(1)

    return lis_config_template, restart_dir, startdate

def _build_restart_filename(restart_dir, startdate):
    """Construct restart filename."""
    restart_file = f"{restart_dir}/LIS_RST_NOAH39_" + \
        f"{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}" + \
        "0000.d01.nc"
    return restart_file

def _check_restart_file(restart_file):
    """Check if restart file exists"""
    if not os.path.exists(restart_file):
        print(f"[ERR] Restart file {restart_file} does not exist!")
        sys.exit(1)

def _customize_lis_config(lis_config_template, restart_dir, startdate):
    """Build new lis.config file customized for current startdate."""

    restart_file = _build_restart_filename(restart_dir, startdate)
    _check_restart_file(restart_file)

    enddate = startdate + datetime.timedelta(days=1)

    # Build dictionary storing replacements for certain lines
    linedict = {
        "Start mode:" : "Start mode: restart\n",
        "Starting year:" : f"Starting year: {startdate.year}\n",
        "Starting month:" : f"Starting month: {startdate.month}\n",
        "Starting day:" : f"Starting day: {startdate.day}\n",
        "Starting hour:" : "Starting hour: 0\n",
        "Starting minute:" : "Starting minute: 0\n",
        "Starting second:" : "Starting second: 0\n",
        "Ending year:" : f"Ending year: {enddate.year}\n",
        "Ending month:" : f"Ending month: {enddate.month}\n",
        "Ending day:" : f"Ending day: {enddate.day}\n",
        "Ending hour:" : "Ending hour: 0\n",
        "Ending minute:" : "Ending minute: 5\n",
        "Ending second:" : "Ending second: 0\n",
        "Noah.3.9 restart file:" : \
            f"Noah.3.9 restart file: {restart_file}\n",
    }

    newlines = [] # List of lines for new file

    with open(lis_config_template, "r", encoding="ascii") as fobj:
        lines = fobj.readlines()
        for line in lines:
            for key, value in linedict.items():
                if key in line:
                    line = value
            newlines.append(line)

    newfile = "lis.config"
    print(f"[INFO] Writing {newfile} customized for new LIS run")
    with open(newfile, "w", encoding="ascii") as fobj:
        for line in newlines:
            fobj.write(line)
        fobj.close()

def _driver():
    """Main driver."""
    lis_config_template, restart_dir, startdate = _read_cmd_args()
    _customize_lis_config(lis_config_template, restart_dir, startdate)

# Invoke driver
if __name__ == "__main__":
    _driver()
