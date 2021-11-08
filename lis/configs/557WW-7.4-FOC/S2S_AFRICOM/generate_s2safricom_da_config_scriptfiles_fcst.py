#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: generate_s2safricom_da_config_scriptfiles_fcst.py
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
import configparser
import datetime
import os
import shutil
#import subprocess
import sys

# Local constants
_NUM_ENSMEMBERS = {
    "CCM4" : 10,
    "CCSM4" : 10,
    "CFSv2" : 24,
    "GEOSv2" : 10,
    "GFDL" : 15,
    "GNEMO" : 10,
}

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} configfile nmme_model"
    print(txt)
    print("[INFO] where:")
    print("[INFO]   configfile: Path to config file.")
    print("[INFO]   nmme_model: Name of NMME model for LIS forcing.")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 3:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    configfile = sys.argv[1]
    if not os.path.exists(configfile):
        print(f"[ERR] Config file {configfile} does not exist!")
        sys.exit(1)

    nmme_model = sys.argv[2]
    if nmme_model not in list(_NUM_ENSMEMBERS):
        print(f"[ERR] Invalid NMME model {nmme_model}!")
        sys.exit(1)

    return configfile, nmme_model

def _read_config(configfile):
    """Read from config file."""
    config = configparser.ConfigParser()
    config.read(configfile)
    return config

def _chdir(config):
    """Change to working directory."""
    print("[INFO] Submitting NRT DA -- AFRICOM Run...")
    print(f"[INFO] Working directory: {config['lisda']['workdir']}")
    try:
        os.chdir(config['lisda']['workdir'])
    except OSError:
        txt = \
            "[ERR] Problem changing to directory " + \
            f"{config['lisda']['workdir']}"
        print(txt)
        sys.exit(1)

def _handle_dates(config):
    """Collect and return date information."""
    currentdate = datetime.date.today()
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
    input_numfcstmons = int(config['lisda']['input_numfcstmons'])
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
                         nmme_model, fcstdir):
    """Create customized lis.config file"""
    with open(lisconfig_target, "rt", encoding='ascii') as file_obj:
        data = file_obj.read()
    data = data.replace("NUMPROCX", f"{config['lisda']['numprocx']}")
    data = data.replace("NUMPROCY", f"{config['lisda']['numprocy']}")
    data = data.replace("START_YEAR", f"{dates['start'].year:04d}")
    data = data.replace("START_MONTH", f"{dates['start'].month:02d}")
    data = data.replace("END_YEAR", f"{dates['end'].year:04d}")
    data = data.replace("END_MONTH", f"{dates['end'].month:02d}")
    data = data.replace("OUTPUTDIR", \
                        f"output/{config['lisda']['datatype']}/" + \
                        f"{dates['start'].strftime('%b')}1/{nmme_model}")
    data = data.replace("LISLOGFILES", \
                        f"output/{config['lisda']['datatype']}/" + \
                        f"{dates['start'].strftime('%b')}1/{nmme_model}" + \
                        f"/logs/{dates['start'].year}/lislog")
    data = data.replace("HYMAPLOGFILES", \
                        f"output/{config['lisda']['datatype']}/" + \
                        f"{dates['start'].strftime('%b')}1/{nmme_model}" + \
                        f"/logs/{dates['start'].year}/hymaplog")
    data = data.replace("NUMENSMEMBERS", \
                        f"{_NUM_ENSMEMBERS[nmme_model]}")

    rst_date = f"{dates['prev'].year:04d}{dates['prev'].month:02d}" + \
               f"{dates['prev'].day:02d}"
    rst_monname = dates['start'].strftime("%b")
    lis_rstdir = f"./input/LDT_ICs/{nmme_model}"
    lis_rstfile = \
        f"{lis_rstdir}/LIS_RST_{config['lisda']['lsmname'].upper()}_{rst_date}" + \
        f"2345.ICS_{rst_monname}{dates['start'].year:04d}." + \
        f"ens{_NUM_ENSMEMBERS[nmme_model]}.nc"
    data = data.replace("LISRSTFILE", f"{lis_rstfile}")

    hymap_rstdir = f"./input/LDT_ICs/{nmme_model}"
    hymap_rstfile = \
            f"{hymap_rstdir}/LIS_RST_{config['lisda']['routingname'].upper()}" + \
            f"_router_{rst_date}2345.ICS_{rst_monname}" + \
            f"{dates['start'].year:04d}." + \
            f"ens1.nc"
    data = data.replace("HYMAPRSTFILE", f"{hymap_rstfile}")

    data = data.replace("FCSTDIR", f"{fcstdir}")

    data = data.replace("LDTINPUTFILE", f"{config['lisda']['ldtinputfile']}")
    with open(lisconfig_target, "wt", encoding='ascii') as file_obj:
        data = file_obj.write(data)

def _customize_slurmjob(slurmjob_target, config, dates, lisconfig_target, \
                        nmme_model):
    """Create customized SLURM job script."""

    with open(slurmjob_target, "rt", encoding='ascii') as file_obj:
        data = file_obj.read()
    data = data.replace("FCSTDATEDIR", \
                        f"{dates['start'].year:04d}{dates['start'].month:02d}")
    data = data.replace("TARGETLISCONFIG", \
                        f"{lisconfig_target}")

    totnumproc = int(config['lisda']['numprocx']) * \
                 int(config['lisda']['numprocy'])
    data = data.replace("TOTNUMPROC", \
                        f"{totnumproc}")

    data = data.replace("LISEXEC", f"{config['lisda']['lisexec']}")
    data = data.replace("MODEL", f"{nmme_model}")
    with open(slurmjob_target, "wt", encoding='ascii') as file_obj:
        file_obj.write(data)

def _driver():
    """Main driver."""
    configfile, nmme_model = _read_cmd_args()
    config = _read_config(configfile)
    _chdir(config)
    dates = _handle_dates(config)

    # Seasonal forecasts
    fcstdir = \
        f"../../{config['lisda']['datatype']}/NMME/final/6-Hourly/{nmme_model}"

    lisconfig_nameconv = \
        f"lis.config.s2safricom.{config['lisda']['lsmname']}." + \
        f"{config['lisda']['routingname']}." + \
        f"da_ics_{config['lisda']['datatype']}_" + \
        f"{nmme_model}_"
    slurmjob_nameconv = \
        f"slurm_s2s_DA_{config['lisda']['datatype']}_{nmme_model}"

    final_outputdir = f"input/lis.config_files/{nmme_model}"
    if not os.path.exists(final_outputdir):
        os.makedirs(final_outputdir)

    # Now customize the lis.config file
    template_dir = "./input/template_files"
    lisconfig_template = \
        f"{template_dir}/template_lis.config.s2safricom." + \
        f"{config['lisda']['lsmname']}." + \
        f"{config['lisda']['routingname']}.da_ics_" + \
        f"{config['lisda']['datatype']}"
    print(lisconfig_template)
    lisconfig_target = \
        f"{final_outputdir}/{lisconfig_nameconv}" + \
        f"{dates['start'].year:04d}{dates['start'].month:02d}"
    shutil.copy(lisconfig_template, lisconfig_target)
    _customize_lisconfig(lisconfig_target, config, dates, \
                    nmme_model, fcstdir)

    # Now customize the slurmjob file
    slurmjob_template = \
        f"{template_dir}/template_slurm_s2s_DA_" + \
        f"{config['lisda']['datatype']}_runtime"
    if _NUM_ENSMEMBERS[nmme_model] <= 12:
        slurmjob_template += "-12hr"
    else:
        slurmjob_template += "-24hr"
    slurmjob_target = \
        f"{slurmjob_nameconv}{dates['start'].year:04d}" + \
        f"{dates['start'].month:02d}.job"
    shutil.copy(slurmjob_template, slurmjob_target)
    _customize_slurmjob(slurmjob_target, config, dates, lisconfig_target, \
                        nmme_model)

    # Submit the job
    os.chdir("..")
    cmd = f"sbatch {slurmjob_target}"
    print(cmd)
    #subprocess.run(cmd, shell=True, check=True)

    print("[INFO] Done generating LIS config files and SLURM script files.")

if __name__ == "__main__":
    _driver()
