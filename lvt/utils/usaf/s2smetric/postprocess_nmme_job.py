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
import time

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
#_BATCH_SCRIPT += "/run_Convert_Dyn_FCST_postproc.scr"
_BATCH_SCRIPT += "/run_generate_metrics.sh"

_PYLIBDIR = "/discover/nobackup/projects/lis_aist17/emkemp/AFWA"
_PYLIBDIR += "/lis74_s2s_patches/LISF/lvt/utils/usaf/s2smetric"
_PYLIBDIR += "/lib_bcsd_metrics"

_BASEOUTDIR = "/discover/nobackup/projects/lis_aist17/emkemp/AFWA"
_BASEOUTDIR += "/lis74_s2s_patches/work/POST/forecasts"

_METRIC_VARS = ["RootZone-SM", "Streamflow", "Surface-SM"]

_METRICS = ["RootZone_SM_ANOM", "RootZone_SM_SANOM",
            "Streamflow_ANOM",  "Streamflow_SANOM",
            "Surface_SM_ANOM",  "Surface_SM_SANOM"]

# Local methods
def _usage():
    """Print command line usage."""
    argv = sys.argv[0]
    txt = f"[INFO] Usage: {argv}"
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
    year = currentdate.year
    month = currentdate.month
    print(f"[INFO] Current year / month: {year:04d} / {month:02d}")
    return currentdate

def _submit_metric_batch_jobs(currentdate, model):
    """Submit batch jobs for calculationg metrics."""

    for py_script in ["convert_dyn_fcst_to_anom.py",
                      "convert_dyn_fcst_to_sanom.py"]:
        cmd = "sbatch"
        cmd += f" {_BATCH_SCRIPT}"
        cmd += f" {py_script}"
        cmd += f" {currentdate.month:02d}"
        cmd += f" {_LSM_MODEL}"
        cmd += f" {_LEAD}"
        cmd += f" {_DOMAIN}"
        cmd += f" {currentdate.year:04d}"
        cmd += f" {model}"
        cmd += f" {_CSYR:04d}"
        cmd += f" {_CEYR:04d}"
        cmd += f" {_PYLIBDIR}"
        print(cmd)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem submitting anomaly batch script!")
            sys.exit(1)

def _run_convert_s2s_anom_cf(currentdate):
    """Automate convert_s2s_anom_cf.py for each NMME run."""
    cfoutdir = f"{_BASEOUTDIR}/metrics_cf/{_DOMAIN}/{_LSM_MODEL}"
    if not os.path.exists(cfoutdir):
        os.makedirs(cfoutdir)
    year = currentdate.year
    month = currentdate.month
    for nmme_model in _NMME_MODELS:
        for anom_type in ("anom", "sanom"):
            touchfile = f"{_BASEOUTDIR}/DYN_"
            touchfile += f"{anom_type.upper()}"
            touchfile += f"/{_DOMAIN}"
            touchfile += f"/{_LSM_MODEL}"
            touchfile += f"/{anom_type}.{_LSM_MODEL}.{nmme_model}.done"
            while not os.path.exists(touchfile):
                print(f"[INFO] Waiting for {touchfile}... " + time.asctime() )
                time.sleep(30)
            for metric_var in _METRIC_VARS:
                metricfile = os.path.dirname(touchfile)
                metricfile += f"/{nmme_model}_{metric_var}"
                metricfile += f"_{anom_type.upper()}_init_monthly_"
                metricfile += f"{month:02d}_{year:04d}.nc"
                cmd = f"{_RUNDIR}/convert_s2s_anom_cf.py"
                cmd += f" {metricfile} {cfoutdir}"
                print(cmd)
                if subprocess.call(cmd, shell=True) != 0:
                    print("[ERR] Problem running convert_s2s_anom_cf.py!")
                    sys.exit(1)

def _calc_enddate(startdate):
    """Calculates end date based on number of forecast months"""
    count = 0
    year = startdate.year
    month = startdate.month
    while count < _LEAD:
        count += 1
        month += 1
        if month > 12:
            month = 1
            year += 1
    enddate = datetime.date(year=year, month=month, day=1)
    return enddate

def _run_merge_s2s_anom_cf(currentdate):
    """Automate merge_s2s_anom_cf.py"""
    input_dir = f"{_BASEOUTDIR}/metrics_cf/{_DOMAIN}/{_LSM_MODEL}"
    output_dir = input_dir
    startdate = datetime.date(year=currentdate.year,
                              month=currentdate.month,
                              day=1)
    enddate = _calc_enddate(startdate)
    for nmme_model in _NMME_MODELS:
        cmd =  f"{_RUNDIR}/merge_s2s_anom_cf.py"
        cmd += f" {input_dir} {output_dir}"
        cmd += f" {startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
        cmd += f" {enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
        cmd += f" {nmme_model}"
        print(cmd)
        if subprocess.call(cmd, shell=True) != 0:
            print("[ERR] Problem calling merge_s2s_anom_cf.py!")
            sys.exit(1)

def _run_make_s2s_median_metric_geotiff():
    """Automate make_s2s_median_metric_geotiff.py"""
    input_dir = f"{_BASEOUTDIR}/metrics_cf/{_DOMAIN}/{_LSM_MODEL}"
    for metric in _METRICS:
        cmd = f"{os.path.dirname(_BATCH_SCRIPT)}"
        cmd += "/make_s2s_median_metric_geotiff.py"
        cmd += f" {input_dir} {metric}"
        print(cmd)
        if subprocess.call(cmd, shell=True) != 0:
            print("[ERR] Problem running make_s2s_median_metric_geotiff.py")
            sys.exit(1)

def _driver():
    """Main driver"""
    _read_cmd_args()
    currentdate = _handle_dates()
    if not os.path.exists(_RUNDIR):
        print(f"[ERR] {_RUNDIR} does not exist!")
        sys.exit(1)
    os.chdir(_RUNDIR)
    for model in _NMME_MODELS:
        _submit_metric_batch_jobs(currentdate, model)
        time.sleep(1) # Don't overwhelm SLURM.
    _run_convert_s2s_anom_cf(currentdate)
    _run_merge_s2s_anom_cf(currentdate)
    _run_make_s2s_median_metric_geotiff()

if __name__ == "__main__":
    _driver()
