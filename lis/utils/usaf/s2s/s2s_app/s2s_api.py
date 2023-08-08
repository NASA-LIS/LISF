#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

import sys
import argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
import yaml
#pylint: disable=wrong-import-position
#pylint: disable=import-error

'''
This script writes SLURM job files.
'''

PARSER = argparse.ArgumentParser()
PARSER.add_argument('-f', '--JOBFILE', required=False, help='job file name')
PARSER.add_argument('-t', '--NTASKS', required=False, help='NTASKS')
PARSER.add_argument('-c', '--CONFIGFILE', required=False, help='config file name')
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
if ARGS.CONFIGFILE is not None:
    CONFIGFILE = ARGS.CONFIGFILE
    with open(CONFIGFILE, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)
    sys.path.append(config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')

REPORT = ARGS.REPORT
if REPORT is not None:
    from s2s_modules.shared import utils
    CWD = ARGS.CWD
    YYYYMMDIR = ARGS.YYYYMMDIR
    utils.print_status_report (CWD, YYYYMMDIR)
    sys.exit()

if ARGS.CFSV2_FILE is not None:
    from s2s_modules.bcsd_fcst.bcsd_library import convert_forecast_data_to_netcdf as cfdn
    ds = cfdn.wgrib2_to_netcdf(ARGS.CFSV2_FILE)
    date_str = np.datetime_as_string(ds['valid_time'][-1].values, unit='s')
    grib_lastdate = datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S')
    fcst_lastdate = datetime.strptime(ARGS.YYYYMMDIR, "%Y%m") + relativedelta(months=9)
    if grib_lastdate >= fcst_lastdate:
        sys.exit(0)
    else:
        sys.exit(1)

NTASKS = ARGS.NTASKS
JOBFILE = ARGS.JOBFILE
RUN_LIS = ARGS.RUN_LIS

if RUN_LIS is not None:
    from s2s_modules.shared import utils
    HOURS = ARGS.HOURS
    JOBNAME = ARGS.JOBNAME
    CWD = ARGS.CWD
    utils.job_script_lis(CONFIGFILE, JOBFILE, JOBNAME, CWD, hours = str(HOURS))
    sys.exit()

if NTASKS is None:
    from s2s_modules.shared import utils
    SCHEDULE_FILE = ARGS.SCHEDULE_FILE
    MYID = ARGS.MYID
    AFTERID = ARGS.AFTERID
    utils.update_job_schedule(SCHEDULE_FILE, MYID, JOBFILE, AFTERID)
else:
    from s2s_modules.shared import utils
    HOURS = ARGS.HOURS
    JOBNAME = ARGS.JOBNAME
    CWD = ARGS.CWD
    utils.job_script(CONFIGFILE, JOBFILE, JOBNAME, NTASKS, str(HOURS), CWD)
