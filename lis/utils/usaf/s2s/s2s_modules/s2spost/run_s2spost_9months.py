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
# SCRIPT: run_s2spost_9months.py
#
# PURPOSE: Loops through multiple months of a completed LIS S2S forecast, and
# submits batch jobs (one per month) to generate CF-convention netCDF files.
#
# REVISION HISTORY:
# 24 Sep 2021: Eric Kemp (SSAI), first version.
# 27 Oct 2021: Eric Kemp/SSAI, address pylint string objections.
# 29 Oct 2021: Eric Kemp/SSAI, add config file.
# 02 Jun 2023: K. Arsenault, updated the s2spost filenaming conventions
#
#------------------------------------------------------------------------------
"""


# Standard modules
import sys
import argparse
import datetime
import yaml
# pylint: disable=too-many-locals, import-outside-toplevel
def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} -c configfile -y fcst_year -m fcst_mon \
                        -w cwd  -j job_name -t ntasks -H hours -M NMME_MODEL"
    txt += " model_forcing [--collect_output]"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]   configfile: s2s config file")
    print("[INFO]   fcst_mon is month of start of LIS forecast")
    print("[INFO]   fcst_year is year of start of LIS forecast")
    print("[INFO]   collect_output] is option to collect files into")
    print("[INFO] job_name: SLURM job_name")
    print("[INFO] ntasks: SLURM ntasks")
    print("[INFO] hours: SLURM time hours")
    print("[INFO] nmme_model: NMME model name")

def _advance_date_by_month(curdate):
    """Calculate new date one month in advance."""
    if curdate.month == 12:
        newdate = datetime.date(year=(curdate.year + 1),
                                month=1,
                                day=1)
    else:
        newdate = datetime.date(year=curdate.year,
                                month=(curdate.month + 1),
                                day=1)
    return newdate

def _submit_batch_jobs(args):
    """Submit batch jobs for processing LIS forecast."""
    configfile = args.configfile
    fcst_year = args.fcst_year
    fcst_mon = args.fcst_mon
    job_name = args.jobname
    ntasks = args.ntasks
    hours = args.hours
    cwd = args.cwd
    model_forcing = args.nmme_model

    startdate = datetime.datetime(int(fcst_year), int(fcst_mon), day=1)
    topdatadir = cwd + '/' + model_forcing + '/'

    # load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    total_months = int(config["EXP"]["lead_months"])
    scriptdir = config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/s2s_modules/s2spost/'

    sys.path.append(config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.shared import utils

    # One batch job per month (call to run_s2spost_1month.py):
    fcstdate = startdate
    curdate = startdate
    for _ in range(0, total_months):
        txt = "[INFO] Submitting batch job for"
        txt += f" cf_{model_forcing}_{curdate.year:04d}{curdate.month:02d}"
        print(txt)
        cmd = "python"
        cmd += f" {scriptdir}/run_s2spost_1month.py"
        cmd += f" {configfile} {topdatadir}"
        cmd += f" {fcstdate.year:04d}{fcstdate.month:02d}"
        cmd += f" {curdate.year:04d}{curdate.month:02d} {model_forcing}"
        jobfile = job_name + '_' + model_forcing +  '_' + curdate.strftime("%Y") \
            + curdate.strftime("%m")+  '_run.j'
        jobname = job_name + '_' + model_forcing +  '_' + curdate.strftime("%Y") \
            + curdate.strftime("%m")+ '_'
        utils.job_script(configfile, jobfile, jobname, ntasks, hours, cwd, in_command=cmd)

        newdate = _advance_date_by_month(curdate)
        curdate = newdate

def _driver():
    """Main driver."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--fcst_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--fcst_mon', required=True, help='initial forecast month')
    parser.add_argument('-c', '--configfile', required=True, help='config file name')
    parser.add_argument('-j', '--jobname', required=True, help='SLURM job name')
    parser.add_argument('-t', '--ntasks', required=True, help='number of SLURM tasks')
    parser.add_argument('-H', '--hours', required=True, help='SLURM job hours')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-M', '--nmme_model', required=True, help='NMME model')

    args = parser.parse_args()
    _submit_batch_jobs(args)

# Invoke driver
if __name__ == "__main__":
    _driver()
