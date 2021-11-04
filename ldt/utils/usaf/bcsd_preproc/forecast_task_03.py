#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_03.py
#
# PURPOSE: Reorganizes the downloaded nmme data into a format for further
# processing. Based on FORECAST_TASK_03.sh.
#
# REVISION HISTORY:
# 24 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""

# Standard modules
import os
import subprocess
import sys

# Local constants.  FIXME:  Put in single location for whole system

# Path of the main project directory:
_PROJDIR = '/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'

# Path of the directory where all the download codes are kept:
_SRCDIR = "{}/scripts/code_library".format(_PROJDIR)
#
# Path for where raw and bias corrected forecast files are located:
#
_FORCEDIR="{}/data/forecast/NMME".format(_PROJDIR)
NMME_DOWNLOAD_DIR="{}/raw/download".format(_FORCEDIR)
NMME_OUTPUT_DIR="{}/raw/Monthly".format(_FORCEDIR)
SUPPLEMENTARY_DIR="{}/supplementary_files".format(_SRCDIR)

# Local methods
def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} currentmon1 currentyear".format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] currentmon1: Current integer month")
    print("[INFO] currentyear: Current year")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 3:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    try:
        currentmon1 = int(sys.argv[1])
    except ValueError:
        print("[ERR] Invalid argument for currentmon1!  Received {}" \
            .format(sys.argv[1]))
        _usage()
        sys.exit(1)
    if currentmon1 < 0:
        print("[ERR] Invalid argument for currentmon1!  Received {}" \
              .format(sys.argv[1]))
        _usage()
        sys.exit(1)

    try:
        currentyear = int(sys.argv[2])
    except ValueError:
        print("[ERR] Invalid argument for currentyear!  Received %s" \
              %(sys.argv[2]))
        _usage()
        sys.exit(1)
    if currentyear < 0:
        print("[ERR] Invalid argument for currentyear!  Received %s" \
              %(sys.argv[2]))
        _usage()
        sys.exit(1)

    return currentmon1, currentyear

def _driver():
    """Main driver."""
    currentmon1, currentyear = _read_cmd_args()
    cmd = "python"
    cmd += " {srcdir}/nmme_reorg_f.py".format(srcdir=_SRCDIR)
    cmd += " {currentmon1}".format(currentmon1=currentmon1)
    cmd += " {currentyear}".format(currentyear=currentyear)
    cmd += " {NMME_DOWNLOAD_DIR}".format(NMME_DOWNLOAD_DIR=NMME_DOWNLOAD_DIR)
    cmd += " {NMME_OUTPUT_DIR}".format(NMME_OUTPUT_DIR=NMME_OUTPUT_DIR)
    cmd += " {SUPPLEMENTARY_DIR}".format(SUPPLEMENTARY_DIR=SUPPLEMENTARY_DIR)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print("[ERR] Problem calling python script: {}".format(cmd))
        sys.exit(1)

if __name__ == "__main__":
    _driver()
