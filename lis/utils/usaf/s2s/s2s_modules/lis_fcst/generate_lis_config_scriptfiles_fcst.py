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
# SCRIPT: generate_s2sglobal_da_CONFIG_scriptfiles_fcst.py
#
# PURPOSE: Generates both lis.config files and SLURM-based run job script
# files. Based on earlier shell scripts developed/modified by Amy McNally,
# Jossy Jacob, Kristi Arsenault, and Ryan Zamora.
#
# REVISION HISTORY:
# 3 Nov 2021: Eric Kemp/SSAI, first version.
#
#------------------------------------------------------------------------------
"""


# Standard modules
import os
import sys
import platform
import argparse
import datetime
import shutil
import numpy as np
import yaml
# pylint: disable=consider-using-f-string, too-many-arguments, too-many-locals, import-outside-toplevel
# Local constants
_NUM_ENSMEMBERS = {}
_NMME_MODELS = []
FORECAST_YEAR = 0
FORECAST_MONTH = 0
WORKDIR = ''
CONFIGFILE = ''
JOB_NAME = ''

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} -c CONFIGFILE \
            -w cwd -j JOB_NAME -t NTASKS -H hours -y year -m month"
    print(txt)
    print("[INFO] where:")
    print("[INFO]   CONFIGFILE: Path to config file.")
    print("[INFO] JOB_NAME: SLURM JOB_NAME")
    print("[INFO] hours: SLURM time hours")
    print("[INFO] cwd: current working directory")
    print("[INFO] year: forecast start year")
    print("[INFO] month: forecast start month")

def _handle_dates(year, month, input_numfcstmons):
    """Collect and return date information."""
    currentdate = datetime.date(year, month, 1) #today()
    print("[INFO] Current year / month: " + \
          f"{currentdate.year:04d} / {currentdate.month:02d}")
    startdate = datetime.date(currentdate.year, currentdate.month, 1)
    print(f"[INFO] Start Month Name: {startdate.strftime('%b')}")
    prevdate = startdate - datetime.timedelta(days=1)
    print("[INFO] Rst year / month / day: " + \
          f"{prevdate.year:04d} / {prevdate.month:02d} / {prevdate.day:02d}")
    count = 0
    month = startdate.month
    year = startdate.year

    while count < input_numfcstmons:
        month += 1
        if month > 12:
            month = 1
            year += 1
        count += 1
    enddate = datetime.date(year, month, 1)
    dates = {
        "start" : startdate,
        "end" : enddate,
        "prev" : prevdate,
    }
    return dates

def _customize_lisconfig(lisconfig_target, config, dates, \
                         nmme_model, fcstdir, ic_date = None, jobid = None):
    """Create customized lis.config file"""
    with open(lisconfig_target, "rt", encoding='ascii') as file_obj:
        data = file_obj.read()
    data = data.replace("NUMPROCX", f"{config['FCST']['numprocx']}")
    data = data.replace("NUMPROCY", f"{config['FCST']['numprocy']}")
    data = data.replace("START_YEAR", f"{dates['start'].year:04d}")
    data = data.replace("START_MONTH", f"{dates['start'].month:02d}")
    data = data.replace("END_YEAR", f"{dates['end'].year:04d}")
    data = data.replace("END_MONTH", f"{dates['end'].month:02d}")

    if ic_date is None:
        # The entire forecast is run by a single job
        rst_date = f"{dates['prev'].year:04d}{dates['prev'].month:02d}" + \
            f"{dates['prev'].day:02d}"
        data = data.replace("OUTPUTDIR", \
                            f"{dates['start'].year:04d}{dates['start'].month:02d}/{nmme_model}")
        data = data.replace("LISLOGFILES", \
                            f"{dates['start'].year:04d}{dates['start'].month:02d}/{nmme_model}" + \
                            "/logs/lislog")
        data = data.replace("HYMAPLOGFILES", \
                            f"{dates['start'].year:04d}{dates['start'].month:02d}/{nmme_model}" + \
                            "/logs/hymaplog")
        data = data.replace("IC_YYYYMMDD", \
                            f"{dates['start'].year:04d}{dates['start'].month:02d}01")

        rst_monname = dates['start'].strftime("%b")
        lis_rstdir = f"./input/LDT_ICs/{nmme_model}"

        lis_rstfile = \
            f"{lis_rstdir}/LIS_RST_{config['EXP']['lsm'].upper()}_{rst_date}" + \
            f"2345.ICS_{rst_monname}{dates['start'].year:04d}." + \
            f"ens{_NUM_ENSMEMBERS[nmme_model]}.nc"

        hymap_rstdir = f"./input/LDT_ICs/{nmme_model}"
        hymap_rstfile = \
            f"{hymap_rstdir}/LIS_RST_{config['EXP']['routing_name'].upper()}" + \
            f"_router_{rst_date}2345.ICS_{rst_monname}" + \
            f"{dates['start'].year:04d}." + \
            "ens1.nc"

    else:
        rst_date = f"{dates['start'].year:04d}{dates['start'].month:02d}" + \
            f"{dates['start'].day:02d}"
        data = data.replace("OUTPUTDIR", \
                            f"{ic_date['start'].year:04d}{ic_date['start'].month:02d}/{nmme_model}")
        if jobid is None:
            data = data.replace("LISLOGFILES", \
                            f"{ic_date['start'].year:04d}{ic_date['start'].month:02d}/{nmme_model}"
                            f"/logs/lislog")
            data = data.replace("HYMAPLOGFILES", \
                            f"{ic_date['start'].year:04d}{ic_date['start'].month:02d}/{nmme_model}"
                            f"/logs/hymaplog")
        else:
            data = data.replace("LISLOGFILES", \
                            f"{ic_date['start'].year:04d}{ic_date['start'].month:02d}/{nmme_model}"
                            f"/logs/lislog{jobid}")
            data = data.replace("HYMAPLOGFILES", \
                            f"{ic_date['start'].year:04d}{ic_date['start'].month:02d}/{nmme_model}"
                            f"/logs/hymaplog{jobid}")

        lis_rstfile = \
            f"{ic_date['start'].year:04d}{ic_date['start'].month:02d}/{nmme_model}" \
            f"/SURFACEMODEL/{dates['start'].year:04d}{dates['start'].month:02d}/" + \
            f"/LIS_RST_{config['EXP']['lsm'].upper()}_{rst_date}0000.d01.nc"

        hymap_rstfile = \
            f"{ic_date['start'].year:04d}{ic_date['start'].month:02d}/{nmme_model}" \
            f"/ROUTING/{dates['start'].year:04d}{dates['start'].month:02d}/" + \
            f"/LIS_RST_{config['EXP']['routing_name'].upper()}_router_{rst_date}0000.d01.nc"
        data = data.replace("IC_YYYYMMDD", \
                            f"{ic_date['start'].year:04d}{ic_date['start'].month:02d}01")

    data = data.replace("NUMENSMEMBERS", \
                        f"{_NUM_ENSMEMBERS[nmme_model]}")
    data = data.replace("LISRSTFILE", f"{lis_rstfile}")
    data = data.replace("HYMAPRSTFILE", f"{hymap_rstfile}")
    data = data.replace("FCSTDIR", f"{fcstdir}")
    data = data.replace("LDTINPUTFILE", f"./input/{config['SETUP']['ldtinputfile']}")

    with open(lisconfig_target, "wt", encoding='ascii') as file_obj:
        data = file_obj.write(data)

def _driver(config):
    """Main driver."""

    # import local module
    sys.path.append(config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.shared import utils

    domlabel=''
    if config['EXP']['DOMAIN'] == 'AFRICOM':
        domlabel = 's2safricom'
    elif config['EXP']['DOMAIN'] == 'GLOBAL':
        domlabel = 's2sglobal'

    for nmme_model in _NMME_MODELS:
        # Seasonal forecasts
        fcstdir = \
            f"bcsd_fcst/NMME/final/6-Hourly/{nmme_model}"
        lisconfig_nameconv = \
            f"lis.config.{domlabel}.{config['EXP']['lsm']}." + \
            f"{config['EXP']['routing_name']}." + \
            f"da_ics_{config['SETUP']['DATATYPE']}_" + \
            f"{nmme_model}_"

        final_outputdir = f"input/lis.config_files/{nmme_model}"
        if not os.path.exists(final_outputdir):
            os.makedirs(final_outputdir, exist_ok=True)

        # Now customize the lis.config file
        template_dir = "./input/template_files"
        lisconfig_template = \
            f"{template_dir}/template_lis.config.{domlabel}." + \
            f"{config['EXP']['lsm']}." + \
            f"{config['EXP']['routing_name']}.da_ics_" + \
            f"{config['SETUP']['DATATYPE']}"

        print(lisconfig_template)

        nseg = config['FCST']['JOB_SEGMENTS'][0].get(nmme_model)
        input_numfcstmons = int(config['EXP']['lead_months'])
        lseg = -(-input_numfcstmons//nseg)

        if nseg == 1:
            # The entire forecast is run by a single job
            dates = _handle_dates(FORECAST_YEAR, FORECAST_MONTH, input_numfcstmons)
            lisconfig_target = \
                f"{final_outputdir}/{lisconfig_nameconv}" + \
                f"{dates['start'].year:04d}{dates['start'].month:02d}"
            shutil.copy(lisconfig_template, lisconfig_target)
            _customize_lisconfig(lisconfig_target, config, dates, \
                                 nmme_model, fcstdir)
            jobfile = JOB_NAME + '_' + nmme_model + '_run.j'
            jobname = JOB_NAME + '_' + nmme_model + '_'

            if 'discover' in platform.node() or 'borg' in platform.node():
                mpi_cmd = 'mpirun -np $SLURM_NTASKS ./LIS' + ' -f ' + lisconfig_target
            else:
                mpi_cmd = 'srun ./LIS' + ' -f ' + lisconfig_target

            utils.job_script_lis(CONFIGFILE, jobfile, jobname, WORKDIR,
                                      in_command=mpi_cmd)
        else:
            # The forecast is divided to nseg number of jobs
            slen = np.ones (nseg, dtype = np.int32)*lseg
            dates = []
            slen[nseg-1] = slen[nseg-1] - (np.sum(slen) - input_numfcstmons)
            for i in range (nseg):
                jno = "_{:02d}".format(i)
                if i == 0:
                    dates.append(_handle_dates(FORECAST_YEAR, FORECAST_MONTH, slen[i]))
                    lisconfig_target = \
                        f"{final_outputdir}/{lisconfig_nameconv}" + \
                        f"{dates[0]['start'].year:04d}{dates[0]['start'].month:02d}_" + \
                        "{:02d}".format(i)
                    shutil.copy(lisconfig_template, lisconfig_target)
                    _customize_lisconfig(lisconfig_target, config, dates[0], \
                                         nmme_model, fcstdir)
                else:
                    dates.append(_handle_dates(dates[i-1].get('end').year,
                                               dates[i-1].get('end').month, slen[i]))
                    lisconfig_target = \
                        f"{final_outputdir}/{lisconfig_nameconv}" + \
                        f"{dates[0]['start'].year:04d}{dates[0]['start'].month:02d}_" + \
                        "{:02d}".format(i)
                    shutil.copy(lisconfig_template, lisconfig_target)
                    _customize_lisconfig(lisconfig_target, config, dates[i], \
                                         nmme_model, fcstdir, \
                                         ic_date = dates[0], jobid = '_{:02d}'.format(i))

                jobfile = JOB_NAME + '_' + nmme_model + jno + '_run.j'
                jobname = JOB_NAME + '_' + nmme_model + jno + '_'
                print(lisconfig_target)

                if 'discover' in platform.node() or 'borg' in platform.node():
                    mpi_cmd = 'mpirun -np $SLURM_NTASKS ./LIS' + ' -f ' + lisconfig_target
                else:
                    mpi_cmd = 'srun ./LIS' + ' -f ' + lisconfig_target

                utils.job_script_lis(CONFIGFILE, jobfile, jobname, WORKDIR,
                                 in_command=mpi_cmd)

    print("[INFO] Done generating LIS config files and SLURM script files.")

if __name__ == "__main__":

    # input arguments
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-y', '--year', required=True, help='forecast start year')
    PARSER.add_argument('-m', '--month', required=True, help='forecast start month')
    PARSER.add_argument('-w', '--WORKDIR', required=True, help='working directory')
    PARSER.add_argument('-c', '--CONFIGFILE', required=True, help='s2s.config file')
    PARSER.add_argument('-j', '--JOB_NAME', required=True, help='JOB_NAME')
    ARGS = PARSER.parse_args()

    FORECAST_YEAR = int(ARGS.year)
    FORECAST_MONTH = int(ARGS.month)
    WORKDIR = ARGS.WORKDIR
    CONFIGFILE = ARGS.CONFIGFILE
    JOB_NAME = ARGS.JOB_NAME

    # Read s2s.config
    with open(CONFIGFILE, 'r', encoding="utf-8") as file:
        _CONFIG = yaml.safe_load(file)
    _NMME_MODELS = _CONFIG['EXP']['NMME_models']
    _NUM_ENSMEMBERS = _CONFIG['EXP']['ensemble_sizes'][0]
    print(_NMME_MODELS)
    print(_NUM_ENSMEMBERS)
    _driver(_CONFIG)
