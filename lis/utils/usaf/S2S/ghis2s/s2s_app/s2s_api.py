#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

'''
This script:
 (1) writes Python/LIS job files
 (2) writes JOB_SCHEDULE table
 (3) runs CFSv2 file checker
'''

import sys
import argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
from ghis2s.shared import utils
from ghis2s.bcsd.bcsd_library import convert_forecast_data_to_netcdf as cfdn

def cfsv2_file_checker(cfsv2_file, yyyymmdir, py_call=False):
    ''' calls CFSv2 file checker '''
    ds = cfdn.wgrib2_to_netcdf(cfsv2_file)
    date_str = np.datetime_as_string(ds['valid_time'][-1].values, unit='s')
    grib_lastdate = datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S')
    fcst_lastdate = datetime.strptime(yyyymmdir, "%Y%m") + relativedelta(months=9)
    del ds
    if grib_lastdate >= fcst_lastdate:
        if py_call:
            return 0
        else:
            sys.exit(0)
    else:
        if py_call:
            return 1
        else:
            sys.exit(1)
    return

def lis_job_file(configfile, jobfile, jobname, cwd, hours):
    ''' calls LIS job script writer '''
    utils.job_script_lis(configfile, jobfile, jobname, cwd, hours=hours)

def python_job_file(configfile, jobfile, jobname, ntasks, hours, cwd,
                    group_jobs, parallel_run=None):
    ''' calls python job script writer '''
    if group_jobs is None:
        utils.job_script(configfile, jobfile, jobname, ntasks, hours, cwd, parallel_run)
    else:
        with open(group_jobs, 'r', encoding="utf-8") as file:
            commands = [line.strip() for line in file if line.strip()]
        if parallel_run is not None and 'SKIP_ARG' in parallel_run:
            ntasks = 1
        else:
            ntasks = len(commands)
        utils.job_script(configfile, jobfile, jobname, ntasks, hours, cwd,
                         parallel_run, group_jobs=commands)

def update_job_schedule(schedule_file, myid, jobfile, afterid):
    ''' updates job schedule '''
    utils.update_job_schedule(schedule_file, myid, jobfile, afterid)

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-f', '--JOBFILE', required=False, help='job file name')
    PARSER.add_argument('-t', '--NTASKS', required=False, help='NTASKS')
    PARSER.add_argument('-c', '--CONFIGFILE', required=False, help='config file name')
    PARSER.add_argument('-C', '--group_jobs', required=False,
                        help='list of commands for group jobs')
    PARSER.add_argument('-H', '--HOURS', required=False, help='time HOURS')
    PARSER.add_argument('-j', '--JOBNAME', required=False, help='job-name')
    PARSER.add_argument('-w', '--CWD', required=False, help='current working directory')
    PARSER.add_argument('-m', '--MYID', required=False, help='my job id')
    PARSER.add_argument('-a', '--AFTERID', required=False, help='after id')
    PARSER.add_argument('-s', '--SCHEDULE_FILE', required=False, help='schedule file')
    PARSER.add_argument('-d', '--YYYYMMDIR', required=False, help='yyyymm directory')
    PARSER.add_argument('-L', '--RUN_LIS', required=False, help='running LISF executables')
    PARSER.add_argument('-i', '--CFSV2_FILE', required=False, help='run CFSv2 file-checker')

    ARGS = PARSER.parse_args()

    if ARGS.CFSV2_FILE is not None:
        cfsv2_file_checker(ARGS.CFSV2_FILE, ARGS.YYYYMMDIR)
        sys.exit()

    if ARGS.RUN_LIS is not None:
        lis_job_file(ARGS.CONFIGFILE, ARGS.JOBFILE, ARGS.JOBNAME, ARGS.CWD, str(ARGS.HOURS))
        sys.exit()

    if ARGS.SCHEDULE_FILE is not None:
        update_job_schedule(ARGS.SCHEDULE_FILE, ARGS.MYID, ARGS.JOBFILE, ARGS.AFTERID)
        sys.exit()

    python_job_file(ARGS.CONFIGFILE, ARGS.JOBFILE, ARGS.JOBNAME, ARGS.NTASKS, str(ARGS.HOURS),
                    ARGS.CWD, ARGS.group_jobs)
