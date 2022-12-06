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
# SCRIPT: make_lvt_config_usafsipost.py
#
# PURPOSE: Reads template lvt.config file for running LVT in USAFSIpost mode,
# and customizes the start and end date/time based on user-provided values.
#
# REQUIREMENTS as of 15 Jul 2021:
# * Python 3.8
#
# REVISION HISTORY:
# 15 Jul 2021: Eric Kemp (SSAI), first version.
# 06 Dec 2022: Eric Kemp (SSAI), changes to boost pylint score.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import sys

def _usage():
    """Prints usage statement for script."""
    print(f"USAGE: {sys.argv[0]} template YYYYMMDDHH")
    print("  where:")
    print("    template is path to lvt.config template file")
    print("    YYYYMMDDHH is valid date/time of USAFSI analysis in UTC")

def _main():
    """Main driver"""
    # Check command-line arguments
    if (len(sys.argv)) != 3:
        _usage()
        sys.exit(1)

    template = sys.argv[1]
    if not os.path.exists(template):
        print(f"[ERR] {template} does not exist!")
        sys.exit(1)

    yyyymmddhh = sys.argv[2]
    year = int(yyyymmddhh[0:4])
    month = int(yyyymmddhh[4:6])
    day = int(yyyymmddhh[6:8])
    hour = int(yyyymmddhh[8:10])
    rundt = datetime.datetime(year=year, month=month, day=day, hour=hour)

    # Open the template lvt.config file, read in the contents, and update as
    # needed.
    with open(template, 'r', encoding="ascii") as file:
        lines = file.readlines()
    newlines = []
    config_dict = {
        "Starting year:" : f"Starting year: {rundt.year}\n",
        "Starting month:" : f"Starting month: {rundt.month}\n",
        "Starting day:" : f"Starting day: {rundt.day}\n",
        "Starting hour:" : f"Starting hour: {rundt.hour}\n",
        "Starting minute:" : "Starting minute: 00\n",
        "Starting second:" : "Starting second: 00\n",
        "Ending year:" : f"Ending year: {rundt.year}\n",
        "Ending month:" : f"Ending month: {rundt.month}\n",
        "Ending day:" : f"Ending day: {rundt.day}\n",
        "Ending minute:" : "Ending minute: 00\n",
        "Ending second:" : "Ending second: 00\n"
    }
    for line in lines:
        newline = line
        for _, key in enumerate(config_dict):
            if key in line:
                newline = config_dict[key]
                break
        newlines.append(newline)

    # Create the new, customized lvt.config file
    newline = "lvt.config.usafsipost"
    with open(newline, "w", encoding="ascii") as file:
        for line in newlines:
            file.write(line)

# Main driver
if __name__ == "__main__":
    _main()
