#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

import sys
import argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
import yaml
from shared import utils
from bcsd.bcsd_library import convert_forecast_data_to_netcdf as cfdn
#pylint: disable=wrong-import-position
#pylint: disable=import-error

'''
This script:
 (1) writes Python/LIS job files
 (2) writes JOB_SCHEDULE table
 (3) produces production status report
 (4) runs CFSv2 file checker
'''

def print_status_report(CWD, YYYYMMDIR):
    utils.print_status_report(CWD, YYYYMMDIR)
    sys.exit()

def cfsv2_file_checker(CFSV2_FILE, YYYYMMDIR):
    ds = cfdn.wgrib2_to_netcdf(CFSV2_FILE)
    date_str = np.datetime_as_string(ds['valid_time'][-1].values, unit='s')
    grib_lastdate = datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S')
    fcst_lastdate = datetime.strptime(YYYYMMDIR, "%Y%m") + relativedelta(months=9)
    if grib_lastdate >= fcst_lastdate:
        sys.exit(0)
    else:
        sys.exit(1)

def lis_job_file(CONFIGFILE, JOBFILE, JOBNAME, CWD, HOURS):
    utils.job_script_lis(CONFIGFILE, JOBFILE, JOBNAME, CWD, hours = HOURS)
    sys.exit()

def python_job_file(CONFIGFILE, JOBFILE, JOBNAME, NTASKS, HOURS, CWD, group_jobs):
    if group_jobs is None:
        utils.job_script(CONFIGFILE, JOBFILE, JOBNAME, NTASKS, HOURS, CWD)
    else:
        with open(group_jobs, 'r') as file:
            commands = [line.strip() for line in file if line.strip()]
        NTASKS = len(commands)
        utils.job_script(CONFIGFILE, JOBFILE, JOBNAME, NTASKS, HOURS, CWD, group_jobs=commands)
    sys.exit()
    
def update_job_schedule(SCHEDULE_FILE, MYID, JOBFILE, AFTERID):
    utils.update_job_schedule(SCHEDULE_FILE, MYID, JOBFILE, AFTERID)
    sys.exit()

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-f', '--JOBFILE', required=False, help='job file name')
    PARSER.add_argument('-t', '--NTASKS', required=False, help='NTASKS')
    PARSER.add_argument('-c', '--CONFIGFILE', required=False, help='config file name')
    PARSER.add_argument('-C', '--group_jobs', required=False, help='list of commands for group jobs')
    PARSER.add_argument('-H', '--HOURS', required=False, help='time HOURS')
    PARSER.add_argument('-j', '--JOBNAME', required=False, help='job-name')
    PARSER.add_argument('-w', '--CWD', required=False, help='current working directory')
    PARSER.add_argument('-m', '--MYID', required=False, help='my job id')
    PARSER.add_argument('-a', '--AFTERID', required=False, help='after id')
    PARSER.add_argument('-s', '--SCHEDULE_FILE', required=False, help='schedule file')
    PARSER.add_argument('-r', '--REPORT', required=False, help='print report')
    PARSER.add_argument('-d', '--YYYYMMDIR', required=False, help='yyyymm directory')
    PARSER.add_argument('-L', '--RUN_LIS', required=False, help='running LISF executables')
    PARSER.add_argument('-i', '--CFSV2_FILE', required=False, help='run CFSv2 file-checker')

    ARGS = PARSER.parse_args()

    if ARGS.REPORT is not None:
        print_status_report(ARGS.CWD, ARGS.YYYYMMDIR)

    if ARGS.CFSV2_FILE is not None:
        cfsv2_file_checker(ARGS.CFSV2_FILE, ARGS.YYYYMMDIR)

    if ARGS.RUN_LIS is not None:
        lis_job_file(ARGS.CONFIGFILE, ARGS.JOBFILE, ARGS.JOBNAME, ARGS.CWD, str(ARGS.HOURS))

    if ARGS.SCHEDULE_FILE is not None:
        update_job_schedule(ARGS.SCHEDULE_FILE, ARGS.MYID, ARGS.JOBFILE, ARGS.AFTERID)

    python_job_file(ARGS.CONFIGFILE, ARGS.JOBFILE, ARGS.JOBNAME, ARGS.NTASKS, str(ARGS.HOURS), ARGS.CWD, ARGS.group_jobs)
