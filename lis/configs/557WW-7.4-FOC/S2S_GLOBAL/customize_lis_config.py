#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: customize_lis_config_nafpa.py
#
# PURPOSE: Customizes lis.config file for 1-day run of LIS in AGRMET Ops mode,
# single ensemble member, over global domain.  This will extend the AGRMET
# forcing used by the LIS S2S system.
#
# REQUIREMENTS as of 21 Sep 2021:
# * Python 3.8 or higher.
#
# REVISION HISTORY:
# 21 Sep 2021: Eric Kemp (SSAI), first version.
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import sys

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s lis_config_template restart_dir YYYYMMDD" \
        %(sys.argv[0])
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
        print("[ERR] %s does not exist!" %(lis_config_template))
        sys.exit(1)

    # Check if directory for restart files exists. Actual restart file
    # shall be checked later.
    restart_dir = sys.argv[2]
    if not os.path.exists(restart_dir):
        print("[ERR] Directory %s does not exist!" %(restart_dir))
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
    restart_file = "%s/LIS_RST_NOAH39_%4.4d%2.2d%2.2d0000.d01.nc" \
        %(restart_dir, startdate.year, startdate.month, startdate.day)
    return restart_file

def _check_restart_file(restart_file):
    """Check if restart file exists"""
    if not os.path.exists(restart_file):
        print("[ERR] Restart file %s does not exist!" %(restart_file))
        sys.exit(1)

def _customize_lis_config(lis_config_template, restart_dir, startdate):
    """Build new lis.config file customized for current startdate."""

    restart_file = _build_restart_filename(restart_dir, startdate)
    _check_restart_file(restart_file)

    enddate = startdate + datetime.timedelta(days=1)

    # Build dictionary storing replacements for certain lines
    linedict = {
        "Start mode:" : "Start mode: restart\n",
        "Starting year:" : "Starting year: %s\n" %(startdate.year),
        "Starting month:" : "Starting month: %s\n" %(startdate.month),
        "Starting day:" : "Starting day: %s\n" %(startdate.day),
        "Starting hour:" : "Starting hour: 0\n",
        "Starting minute:" : "Starting minute: 0\n",
        "Starting second:" : "Starting second: 0\n",
        "Ending year:" : "Ending year: %s\n" %(enddate.year),
        "Ending month:" : "Ending month: %s\n" %(enddate.month),
        "Ending day:" : "Ending day: %s\n" %(enddate.day),
        "Ending hour:" : "Ending hour: 0\n",
        "Ending minute:" : "Ending minute: 5\n",
        "Ending second:" : "Ending second: 0\n",
        "Noah.3.9 restart file:" : "Noah.3.9 restart file: %s\n" \
                                   %(restart_file),
    }

    newlines = [] # List of lines for new file
    lines = open(lis_config_template, "r").readlines()
    for line in lines:
        for key in linedict:
            if key in line:
                line = linedict[key]
        newlines.append(line)

    newfile = "lis.config"
    print("[INFO] Writing %s customized for new LIS run" %(newfile))
    fobj = open(newfile, "w")
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
