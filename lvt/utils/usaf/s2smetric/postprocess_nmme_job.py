#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: postprocess_nmme_job.py
#
# PURPOSE: Main script for generating metrics for NMME-forced LIS forecasts.
# Based on Postprocess_NMME_job.sh and job_run_convert_Dyn_FCST_to_postproc.scr
#
# REVISION HISTORY:
# 22 Oct 2021: Eric Kemp/SSAI, first version
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import subprocess
import sys

# Local constants
_CSYR = 2008
_CEYR = 2020

# Temporary.  FIXME Put in common location for all scripts
_RUNDIR = "/discover/nobackup/projects/lis_aist17/emkemp/AFWA"
_RUNDIR += "/lis74_s2s_patches/LISF/lvt/utils/usaf/s2smetric"

#_NMME_MODELS = ["CCM4", "CCSM4", "CFSv2", "GEOSv2", "GFDL", "GNEMO"]
_NMME_MODELS = ["GEOSv2"] # For testing

_LSM_MODEL = "NOAHMP"

_LEAD = 9 # 9 month forecasts

_DOMAIN = "AFRICOM"

_BATCH_SCRIPT = "/discover/nobackup/projects/lis_aist17/emkemp/AFWA"
_BATCH_SCRIPT += "/lis74_s2s_patches/LISF/lvt/utils/usaf/s2smetric"
_BATCH_SCRIPT += "/run_Convert_Dyn_FCST_postproc.scr"

# Local methods
def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s" %(sys.argv[0])
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""
    if len(sys.argv) != 1:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

def _handle_dates():
    """Collect and return data information."""
    currentdate = datetime.date.today()
    print("[INFO] Current year / month: %4.4d / %2.2d" \
          %(currentdate.year, currentdate.month))
    return currentdate

def _submit_metric_batch_jobs(currentdate, model):
    """Submit batch jobs for calculationg metrics."""

    for py_script in ["convert_dyn_fcst_to_anom.py",
                      "convert_dyn_fcst_to_sanom.py"]:
        cmd = "sbatch"
        cmd += " %s" %(_BATCH_SCRIPT)
        cmd += " %s" %(py_script)
        cmd += " %2.2d" %(currentdate.month)
        cmd += " %s" %(_LSM_MODEL)
        cmd += " %d" %(_LEAD)
        cmd += " %s" %(_DOMAIN)
        cmd += " %4.4d" %(currentdate.year)
        cmd += " %s" %(model)
        cmd += " %4.4d" %(_CSYR)
        cmd += " %4.4d" %(_CEYR)
        print(cmd)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem submitting anomaly batch script!")
            sys.exit(1)

def _driver():
    """Main driver"""
    _read_cmd_args()
    currentdate = _handle_dates()
    if not os.path.exists(_RUNDIR):
        print("[ERR] %s does not exist!" %(_RUNDIR))
        sys.exit(1)
    os.chdir(_RUNDIR)
    for model in _NMME_MODELS:
        _submit_metric_batch_jobs(currentdate, model)

if __name__ == "__main__":
    _driver()
