'''
This script writes SLURM job files.
'''
import os
import sys
import argparse
sys.path.append(os.getenv('S2STOOL'))
from s2s_modules.shared import utils

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

ARGS = PARSER.parse_args()
REPORT = ARGS.REPORT
if REPORT is not None:
    CWD = ARGS.CWD
    YYYYMMDIR = ARGS.YYYYMMDIR
    utils.print_status_report (CWD, YYYYMMDIR)
    sys.exit()

NTASKS = ARGS.NTASKS
JOBFILE = ARGS.JOBFILE

if NTASKS is None:
    SCHEDULE_FILE = ARGS.SCHEDULE_FILE
    MYID = ARGS.MYID
    AFTERID = ARGS.AFTERID
    utils.update_job_schedule(SCHEDULE_FILE, MYID, JOBFILE, AFTERID)
else:
    CONFIGFILE = ARGS.CONFIGFILE
    HOURS = ARGS.HOURS
    JOBNAME = ARGS.JOBNAME
    CWD = ARGS.CWD
    utils.job_script(CONFIGFILE, JOBFILE, JOBNAME, NTASKS, str(HOURS), CWD)
