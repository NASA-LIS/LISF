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
# 02 Jun 2023: K. Arsenault + S. Mahanama, updated 557 WW filename conventions.
#
#------------------------------------------------------------------------------
"""




# Standard modules
import os
import sys
import argparse
import datetime
import subprocess
import yaml
# pylint: disable=consider-using-f-string, too-many-locals, import-outside-toplevel
# Local methods
def _usage():
    """Print command line usage."""
    print(f"[INFO] Usage: {sys.argv[0]} configfile")
    print("[INFO]  where:")
    print("[INFO]    configfile: path to s2smetric config file")

def _handle_dates(year, month):
    """Collect and return data information."""
    currentdate = datetime.date(year, month, 1)
    year = currentdate.year
    month = currentdate.month
    print(f"[INFO] Current year / month: {year:04d} / {month:02d}")
    return currentdate

def _run_convert_s2s_anom_cf(config, currentdate, baseoutdir):
    """Automate convert_s2s_anom_cf.py for each NMME run."""
    cfoutdir = \
        f"{baseoutdir}" + \
        f"/metrics_cf/{config['EXP']['lsmdir']}"
    if not os.path.exists(cfoutdir):
        os.makedirs(cfoutdir)

    year = currentdate.year
    month = currentdate.month
    metric_vars = config["POST"]["metric_vars"]
    nmme_models = config["EXP"]["NMME_models"]
    for nmme_model in nmme_models:
        for anom_type in ("anom", "sanom"):
            touchfile = f"{baseoutdir}/DYN_"
            touchfile += f"{anom_type.upper()}"
            touchfile += f"/{config['EXP']['lsmdir']}"
            touchfile += f"/{anom_type}.{config['EXP']['lsmdir']}"
            touchfile += f".{nmme_model}.done"
           #while not os.path.exists(touchfile):
           #    print(f"[INFO] Waiting for {touchfile}... " + time.asctime() )
           #    time.sleep(30)
            rundir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/s2smetric/'
            for metric_var in metric_vars:
                metricfile = os.path.dirname(touchfile)
                metricfile += f"/{nmme_model}_{metric_var}"
                metricfile += f"_{anom_type.upper()}_init_monthly_"
                metricfile += f"{month:02d}_{year:04d}.nc"
                cmd = f"python {rundir}/convert_s2s_anom_cf.py"
                cmd += f" {metricfile} {cfoutdir}"
                print(cmd)
                if subprocess.call(cmd, shell=True) != 0:
                    print("[ERR] Problem running convert_s2s_anom_cf.py!")
                    sys.exit(1)

def _calc_enddate(config, startdate):
    """Calculates end date based on number of forecast months"""
    count = 0
    year = startdate.year
    month = startdate.month
    lead = int(config["EXP"]["lead_months"])
    while count < lead:
        count += 1
        month += 1
        if month > 12:
            month = 1
            year += 1
    enddate = datetime.date(year=year, month=month, day=1)
    return enddate

def _run_merge_s2s_anom_cf(config, currentdate, configfile, baseoutdir):
    """Automate merge_s2s_anom_cf.py"""
    lsm_model = config["EXP"]["lsmdir"]
    input_dir = f"{baseoutdir}/metrics_cf/{lsm_model}"
    output_dir = input_dir
    startdate = datetime.date(year=currentdate.year,
                              month=currentdate.month,
                              day=1)
    enddate = _calc_enddate(config, startdate)
    rundir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/s2smetric/'
    nmme_models = config["EXP"]["NMME_models"]
    for nmme_model in nmme_models:
        cmd = "python"
        cmd += f" {rundir}/merge_s2s_anom_cf.py"
        cmd += f" {input_dir} {output_dir}"
        cmd += f" {startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
        cmd += f" {enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
        cmd += f" {nmme_model} {configfile}"
        print(cmd)
        if subprocess.call(cmd, shell=True) != 0:
            print("[ERR] Problem calling merge_s2s_anom_cf.py!")
            sys.exit(1)

def _run_make_s2s_median_metric_geotiff(config, configfile, baseoutdir):
    """Automate make_s2s_median_metric_geotiff.py"""
    lsm_model = config["EXP"]["lsmdir"]
    input_dir = f"{baseoutdir}/metrics_cf/{lsm_model}"
    rundir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/s2smetric/'

    metrics1 = config["POST"]["metric_vars"]
    metrics = ["{}{}".format(i, '_ANOM') for i in metrics1]
    metrics.extend(["{}{}".format(i, '_SANOM') for i in metrics1])

    for metric1 in metrics:
        metric = metric1.replace('-', '_')
        cmd = "python"
        cmd += f" {rundir}"
        cmd += "/make_s2s_median_metric_geotiff.py"
        cmd += f" {input_dir} {metric} {configfile}"
        print(cmd)
        if subprocess.call(cmd, shell=True) != 0:
            print("[ERR] Problem running make_s2s_median_metric_geotiff.py")
            sys.exit(1)

def _driver():
    """Main driver"""

    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--fcst_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--fcst_mon', required=True, help='forecast end year')
    parser.add_argument('-c', '--configfile', required=True, help='config file name')
    parser.add_argument('-j', '--jobname', help='job_name')
    parser.add_argument('-t', '--ntasks', help='ntasks')
    parser.add_argument('-H', '--hours', help='hours')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-M', '--nmme_model', required=False, help='NMME Model')

    args = parser.parse_args()
    configfile = args.configfile
    fcst_year = args.fcst_year
    fcst_mon = args.fcst_mon
    cwd = args.cwd
    baseoutdir = cwd + '/s2smetric/' + fcst_year  + fcst_mon
    currentdate = _handle_dates(int(fcst_year), int(fcst_mon))
    nmme_model = args.nmme_model
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if nmme_model is None:
        _run_convert_s2s_anom_cf(config, currentdate, baseoutdir)
        _run_merge_s2s_anom_cf(config, currentdate, configfile, baseoutdir)
        _run_make_s2s_median_metric_geotiff(config, configfile, baseoutdir)
    else:
        jobname = args.jobname
        ntasks = args.ntasks
        hours = args.hours
        sys.path.append(config['SETUP']['LISFDIR']  + '/lis/utils/usaf/s2s/')
        from s2s_modules.shared import utils
        pylibdir = config['SETUP']['LISFDIR'] + \
            '/lis/utils/usaf/s2s/s2s_modules/s2smetric/metrics_library/'
        for anom_type in ["anom", "sanom"]:
            py_script = "convert_dyn_fcst_to_" + anom_type + ".py"
            cmd = "python"
            cmd += f" {pylibdir}{py_script}"
            cmd += f" {currentdate.month:02d}"
            cmd += f" {currentdate.year:04d}"
            cmd += f" {nmme_model}"
            cmd += f" {configfile}"
            cmd += f" {baseoutdir}"
            jobfile = jobname + '_' + nmme_model + '_' + anom_type + '_run.j'
            job_name = jobname + '_' + nmme_model + '_' + anom_type + '_'
            print(cmd)
            utils.job_script(configfile, jobfile, job_name, ntasks, hours, cwd, in_command=cmd)

if __name__ == "__main__":
    _driver()
