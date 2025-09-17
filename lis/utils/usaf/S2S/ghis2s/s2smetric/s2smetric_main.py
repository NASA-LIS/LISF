#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: s2smetric_main.py
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
from concurrent.futures import ProcessPoolExecutor
import yaml
from ghis2s.shared import utils
from ghis2s.shared.logging_utils import TaskLogger
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

def driver(configfile, fcst_year, fcst_mon, cwd, nmme_model=None, jobname=None,
         ntasks=None, hours=None, py_call=False, weekly=False):
    """Main driver"""

    baseoutdir = cwd + '/s2smetric/{:04d}{:02d}'.format(fcst_year, fcst_mon)
    currentdate = _handle_dates(fcst_year, fcst_mon)
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    pylibdir = config['SETUP']['LISFDIR'] + \
        '/lis/utils/usaf/S2S/ghis2s/s2smetric/'
    slurm_commands = []
    if not weekly:
        py_script = "convert_dyn_fcst_to_anom.py"
        cmd = "python"
        cmd += f" {pylibdir}{py_script}"
        cmd += f" {currentdate.month:02d}"
        cmd += f" {currentdate.year:04d}"
        cmd += f" {nmme_model}"
        cmd += f" {configfile}"
        cmd += f" {baseoutdir}"
        jobfile = jobname + '_' + nmme_model + '_anom_run.j'
        job_name = jobname + '_' + nmme_model + '_anom_'
        print(cmd)
        if py_call:
            slurm_commands.append(cmd)
        else:
            utils.job_script(configfile, jobfile, job_name, ntasks, hours, cwd, in_command=cmd)
    else:
        py_script = "compute_weekly_anom.py"
        cmd = "python"
        cmd += f" {pylibdir}{py_script}"
        cmd += f" {currentdate.month:02d}"
        cmd += f" {currentdate.year:04d}"
        cmd += f" {nmme_model}"
        cmd += f" {configfile}"
        cmd += f" {baseoutdir}"
        jobfile = jobname + '_' + nmme_model + '_anom_run.j'
        job_name = jobname + '_' + nmme_model + '_anom_'
        print(cmd)
        if py_call:
            slurm_commands.append(cmd)
        else:
            utils.job_script(configfile, jobfile, job_name, ntasks, hours, cwd, in_command=cmd)

    if py_call:
        return slurm_commands

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--fcst_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--fcst_mon', required=True, help='forecast end year')
    parser.add_argument('-c', '--configfile', required=True, help='config file name')
    parser.add_argument('-j', '--jobname', help='job_name')
    parser.add_argument('-t', '--ntasks', help='ntasks')
    parser.add_argument('-H', '--hours', help='hours')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-M', '--nmme_model', required=False, help='NMME Model')
    parser.add_argument('-W', '--weekly', action='store_true', help='weekly metrics (default: False)?')

    args = parser.parse_args()
    if args.nmme_model is None:
        task_name = os.environ.get('SCRIPT_NAME')
        if args.weekly:
            logger = TaskLogger(task_name,
                                os.getcwd(),
                                f's2smetric/postprocess_nmme_job.py-> make_s2s_median_metric_geotiff.py to generate weekly TIF files.')
        else:
            logger = TaskLogger(task_name,
                                os.getcwd(),
                                f's2smetric/postprocess_nmme_job.py-> make_s2s_median_metric_geotiff.py to generate monthly TIF files.')

        def run_tiff_py(cmd):
            if subprocess.call(cmd, shell=True) != 0:
                print("[ERR] Problem running make_s2s_median_metric_geotiff.py")
                sys.exit(1)
                
        baseoutdir = args.cwd + '/s2smetric/{:04d}{:02d}'.format(int(args.fcst_year), int(args.fcst_mon))
        with open(args.configfile, 'r', encoding="utf-8") as file:
            config = yaml.safe_load(file)

        input_dir = f"{baseoutdir}/"
        rundir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/S2S/ghis2s/s2smetric/'

        if args.weekly:
            metrics1 = config["POST"]["weekly_vars"]
        else:
            metrics1 = config["POST"]["metric_vars"]
        num_vars = 2*len(metrics1)
        num_workers = int(os.environ.get('NUM_WORKERS', num_vars))
        metrics = ["{}{}".format(i, '_ANOM') for i in metrics1]
        metrics.extend(["{}{}".format(i, '_SANOM') for i in metrics1])

        # ProcessPoolExecutor parallel processing
        logger.info("Starting parallel processing of variables")
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = []
            for metric1 in metrics:
                metric = metric1.replace('-', '_')
                logger.info(f"Submitting TIF generating job for {metric}", subtask=metric)
                cmd = "python"
                cmd += f" {rundir}"
                cmd += "/make_s2s_median_metric_geotiff.py"
                cmd += f" {input_dir} {metric} {args.configfile}"
                if args.weekly:
                    cmd += f" weekly"
                future = executor.submit(run_tiff_py, cmd)
                futures.append(future)

            for future in futures:
                try:
                    result = future.result()
                except Exception as e:
                    logger.error(f"Failed writing TIF files: {str(e)}", subtask=metric)

        logger.info(f"Writing TIF files completed successfully")
                
    else:
        driver(args.configfile, int(args.fcst_year), int(args.fcst_mon), args.cwd,
             nmme_model=args.nmme_model, jobname=args.jobname, ntasks=args.ntasks,
             hours=args.hours)
