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
# 29 Oct 2021: Eric Kemp/SSAI, add config file
#
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
import datetime
import os
import subprocess
import sys
import time

# Local methods
def _usage():
    """Print command line usage."""
    print(f"[INFO] Usage: {sys.argv[0]} configfile")
    print("[INFO]  where:")
    print("[INFO]    configfile: path to s2smetric config file")

def _read_cmd_args():
    """Read command line arguments."""
    if len(sys.argv) != 2:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    configfile = sys.argv[1]
    if not os.path.exists(configfile):
        print(f"[ERR] Config file {configfile} does not exist!")
        sys.exit(1)

    return configfile

def _read_config(configfile):
    """Read s2smetric config file."""
    config = configparser.ConfigParser()
    config.read(configfile)
    return config

def _handle_dates():
    """Collect and return data information."""
    currentdate = datetime.date.today()
    year = currentdate.year
    month = currentdate.month
    print(f"[INFO] Current year / month: {year:04d} / {month:02d}")
    return currentdate

def _submit_metric_batch_jobs(config, currentdate, nmme_model, configfile):
    """Submit batch jobs for calculating metrics."""

    batch_script = config["s2smetric"]["batch_script"]
    pylibdir = config["s2smetric"]["pylibdir"]
    for py_script in ["convert_dyn_fcst_to_anom.py",
                      "convert_dyn_fcst_to_sanom.py"]:
        cmd = "sbatch"
        cmd += f" {batch_script}"
        cmd += f" {py_script}"
        cmd += f" {currentdate.month:02d}"
        cmd += f" {currentdate.year:04d}"
        cmd += f" {nmme_model}"
        cmd += f" {configfile}"
        cmd += f" {pylibdir}"
        print(cmd)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem submitting anomaly batch script!")
            sys.exit(1)
        time.sleep(1) # Don't overwhelm SLURM

def _run_convert_s2s_anom_cf(config, currentdate, configfile):
    """Automate convert_s2s_anom_cf.py for each NMME run."""
    cfoutdir = \
        f"{config['s2smetric']['baseoutdir']}" + \
        f"/metrics_cf/{config['s2smetric']['domain']}" + \
        f"/{config['s2smetric']['lsm_model']}"
    if not os.path.exists(cfoutdir):
        os.makedirs(cfoutdir)
    year = currentdate.year
    month = currentdate.month
    nmme_models = config["s2smetric"]["nmme_models"].split()
    metric_vars = config["s2smetric"]["metric_vars"].split()
    for nmme_model in nmme_models:
        for anom_type in ("anom", "sanom"):
            touchfile = f"{config['s2smetric']['baseoutdir']}/DYN_"
            touchfile += f"{anom_type.upper()}"
            touchfile += f"/{config['s2smetric']['domain']}"
            touchfile += f"/{config['s2smetric']['lsm_model']}"
            touchfile += f"/{anom_type}.{config['s2smetric']['lsm_model']}"
            touchfile += f".{nmme_model}.done"
            while not os.path.exists(touchfile):
                print(f"[INFO] Waiting for {touchfile}... " + time.asctime() )
                time.sleep(30)
            for metric_var in metric_vars:
                metricfile = os.path.dirname(touchfile)
                metricfile += f"/{nmme_model}_{metric_var}"
                metricfile += f"_{anom_type.upper()}_init_monthly_"
                metricfile += f"{month:02d}_{year:04d}.nc"
                cmd = f"{config['s2smetric']['rundir']}/convert_s2s_anom_cf.py"
                cmd += f" {metricfile} {cfoutdir} {configfile}"
                print(cmd)
                if subprocess.call(cmd, shell=True) != 0:
                    print("[ERR] Problem running convert_s2s_anom_cf.py!")
                    sys.exit(1)

def _calc_enddate(config, startdate):
    """Calculates end date based on number of forecast months"""
    count = 0
    year = startdate.year
    month = startdate.month
    lead = int(config["s2smetric"]["lead"])
    while count < lead:
        count += 1
        month += 1
        if month > 12:
            month = 1
            year += 1
    enddate = datetime.date(year=year, month=month, day=1)
    return enddate

def _run_merge_s2s_anom_cf(config, currentdate, configfile):
    """Automate merge_s2s_anom_cf.py"""
    baseoutdir = config["s2smetric"]["baseoutdir"]
    domain = config["s2smetric"]["domain"]
    lsm_model = config["s2smetric"]["lsm_model"]
    input_dir = f"{baseoutdir}/metrics_cf/{domain}/{lsm_model}"
    output_dir = input_dir
    startdate = datetime.date(year=currentdate.year,
                              month=currentdate.month,
                              day=1)
    enddate = _calc_enddate(config, startdate)
    nmme_models = config["s2smetric"]["nmme_models"].split()
    rundir = config["s2smetric"]["rundir"]
    for nmme_model in nmme_models:
        cmd =  f"{rundir}/merge_s2s_anom_cf.py"
        cmd += f" {input_dir} {output_dir}"
        cmd += f" {startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
        cmd += f" {enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
        cmd += f" {nmme_model} {configfile}"
        print(cmd)
        if subprocess.call(cmd, shell=True) != 0:
            print("[ERR] Problem calling merge_s2s_anom_cf.py!")
            sys.exit(1)

def _run_make_s2s_median_metric_geotiff(config, configfile):
    """Automate make_s2s_median_metric_geotiff.py"""
    baseoutdir = config["s2smetric"]["baseoutdir"]
    domain = config["s2smetric"]["domain"]
    lsm_model = config["s2smetric"]["lsm_model"]
    input_dir = f"{baseoutdir}/metrics_cf/{domain}/{lsm_model}"
    metrics = config["s2smetric"]["metrics"].split()
    batch_script = config["s2smetric"]["batch_script"]
    for metric in metrics:
        cmd = f"{os.path.dirname(batch_script)}"
        cmd += "/make_s2s_median_metric_geotiff.py"
        cmd += f" {input_dir} {metric} {configfile}"
        print(cmd)
        if subprocess.call(cmd, shell=True) != 0:
            print("[ERR] Problem running make_s2s_median_metric_geotiff.py")
            sys.exit(1)

def _driver():
    """Main driver"""
    configfile = _read_cmd_args()
    config = _read_config(configfile)
    currentdate = _handle_dates()
    rundir = config["s2smetric"]["rundir"]
    if not os.path.exists(rundir):
        print(f"[ERR] {rundir} does not exist!")
        sys.exit(1)
    os.chdir(rundir)
    nmme_models = config["s2smetric"]["nmme_models"].split()
    for model in nmme_models:
        _submit_metric_batch_jobs(config, currentdate, model, configfile)
    _run_convert_s2s_anom_cf(config, currentdate, configfile)
    _run_merge_s2s_anom_cf(config, currentdate, configfile)
    _run_make_s2s_median_metric_geotiff(config, configfile)

if __name__ == "__main__":
    _driver()
