#!/usr/bin/env python3
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
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import sys

def _usage():
    """Prints usage statement for script."""
    print("USAGE: %s template YYYYMMDDHH" %(sys.argv[0]))
    print("  where:")
    print("    template is path to lvt.config template file")
    print("    YYYYMMDDHH is valid date/time of USAFSI analysis in UTC")

# Main driver
if __name__ == "__main__":

    # Check command-line arguments
    if (len(sys.argv)) != 3:
        _usage()
        sys.exit(1)

    template = sys.argv[1]
    if not os.path.exists(template):
        print("[ERR] %s does not exist!" %(template))
        sys.exit(1)

    yyyymmddhh = sys.argv[2]
    year = int(yyyymmddhh[0:4])
    month = int(yyyymmddhh[4:6])
    day = int(yyyymmddhh[6:8])
    hour = int(yyyymmddhh[8:10])
    rundt = datetime.datetime(year=year, month=month, day=day, hour=hour)

    # Open the template lvt.config file, read in the contents, and update as
    # needed.
    lines = open(template, 'r').readlines()
    newlines = []
    for line in lines:
        # NOTE:  pylint insists that NEWLINE is a constant and should have
        # UPPER_CASE naming style.  This is likely a bug in pylint, but we
        # will humor it.
        if "Starting year:" in line:
            NEWLINE = "Starting year: %s\n" %(rundt.year)
        elif "Starting month:" in line:
            NEWLINE = "Starting month: %s\n" %(rundt.month)
        elif "Starting day:" in line:
            NEWLINE = "Starting day: %s\n" %(rundt.day)
        elif "Starting hour:" in line:
            NEWLINE = "Starting hour: %s\n" %(rundt.hour)
        elif "Starting minute:" in line:
            NEWLINE = "Starting minute: 00\n"
        elif "Starting second:" in line:
            NEWLINE = "Starting second: 00\n"
        # For simplicity, we set the end date/time same as the start
        elif "Ending year:" in line:
            NEWLINE = "Ending year: %s\n" %(rundt.year)
        elif "Ending month:" in line:
            NEWLINE = "Ending month: %s\n" %(rundt.month)
        elif "Ending day:" in line:
            NEWLINE = "Ending day: %s\n" %(rundt.day)
        elif "Ending hour:" in line:
            NEWLINE = "Ending hour: %s\n" %(rundt.hour)
        elif "Ending minute:" in line:
            NEWLINE = "Ending minute: 00\n"
        elif "Ending second:" in line:
            NEWLINE = "Ending second: 00\n"
        else:
            NEWLINE = line
        newlines.append(NEWLINE)

    # Create the new, customized lvt.config file
    NEWFILE = "lvt.config.usafsipost"
    f = open(NEWFILE, "w")
    for line in newlines:
        f.write(line)
    f.close()
