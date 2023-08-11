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
# SCRIPT: generate_lisda_config_nrt.py
#
# PURPOSE: Generates lis.config file for running NoahMP401 12-ensemble
# data assimilation (DA) based run, along with single member HyMAP2 run.
# Inputs required are: starting YYYYMMDD and ending YYYYMMDD for the run.
#
# REQUIREMENTS as of 14 Oct 2021:
# * Python 3.8 or higher
#
# REVISION HISTORY:
# * 14 Oct 2021: Eric Kemp/SSAI, first version.
# * 25 Oct 2021: Kristi Arsenault/SAIC, adapted for LIS DA run version.
#
#------------------------------------------------------------------------------
"""


# Standard modules
import datetime
import os
import shutil
import sys
import configparser

# Local template directory structure (required):
_TEMPLATE_DIR = "./template_files"
_LISCONFIG_TEMPLATE = f"{_TEMPLATE_DIR}/lis.config_template"

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {(sys.argv[0])} configfile STYYYYMMDD EDYYYYMMDD"
    print(txt)
    print("[INFO] where: ")
    print("[INFO]  configfile: Path to python-based input config file.")
    print("[INFO]  STYYYYMMDD: year/month/day of start of next LIS run.")
    print("[INFO]  EDYYYYMMDD: year/month/day of end of next LIS run.")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        print(len(sys.argv))
        print(sys.argv[:])
        _usage()
        sys.exit(1)

    # Check if python script config file exists.
    configfile = sys.argv[1]
    if not os.path.exists(configfile):
        print(f"[ERR] {configfile} does not exist!")
        sys.exit(1)

    # Get start date of new LIS run.
    yyyymmdd = sys.argv[2]
    if len(yyyymmdd) != 8:
        print("[ERR] Invalid length for start YYYYMMDD, must be 8 characters!")
        sys.exit(1)
    year = int(yyyymmdd[0:4])
    month = int(yyyymmdd[4:6])
    day = int(yyyymmdd[6:8])
    try:
        startdate = datetime.date(year, month, day)
    except ValueError:
        print("[ERR] Invalid start YYYYMMDD passed to script!")
        sys.exit(1)

    # Get end date of new LIS run.
    yyyymmdd = sys.argv[3]
    if len(yyyymmdd) != 8:
        print("[ERR] Invalid length for end YYYYMMDD, must be 8 characters!")
        sys.exit(1)
    year = int(yyyymmdd[0:4])
    month = int(yyyymmdd[4:6])
    day = int(yyyymmdd[6:8])
    try:
        finaldate = datetime.date(year, month, day)
    except ValueError:
        print("[ERR] Invalid end YYYYMMDD passed to script!")
        sys.exit(1)

    return configfile, startdate, finaldate

def _read_config(configfile):
    """Read from LIS DA run config file."""
    config = configparser.ConfigParser()
    config.read(configfile)
    return config

def _create_lisconfig_lsm_target(lisconfig_outdir, startdate):
    """Create name of new lis.config file."""

    lisconfig_lsm_target = f"{lisconfig_outdir}/lis.config_darun_" + \
              f"{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"

    print(f"[INFO] The target LIS DA config file is: {lisconfig_lsm_target}")

    return lisconfig_lsm_target


def _customize_lis_config(config,
                          lisconfig_lsm_target,
                          startdate, finaldate):

    """Customize new lis.config file."""
    shutil.copy(_LISCONFIG_TEMPLATE, lisconfig_lsm_target)

    # Build customized settings for target ldt.config file
    # Assemble LIS runtime log file:
    output_dir = config["lisdarun"]["output_dir"]
    lsm_logfile = f"{output_dir}/logs_{startdate.year:04d}{startdate.month:02d}/lislog"

    input_dir = config["lisdarun"]["input_dir"]

    # Assemble NoahMP401 restart file path:
    noahmprst_fname = f"{input_dir}/SURFACEMODEL/{startdate.year:04d}{startdate.month:02d}" + \
          f"/LIS_RST_NOAHMP401_{startdate.year:04d}{startdate.month:02d}" + \
          f"{startdate.day:02d}0000.d01.nc"
    print(f"[INFO] NoahMP restart file: {noahmprst_fname}")

    # Assemble HyMAP2 restart file path:
    hymaprst_fname = f"{input_dir}/ROUTING/{startdate.year:04d}{startdate.month:02d}" + \
          f"/LIS_RST_HYMAP2_router_{startdate.year:04d}{startdate.month:02d}" + \
          f"{startdate.day:02d}0000.d01.nc"
    print(f"[INFO] HyMAP2 restart file: {hymaprst_fname}")

    # Assemble DA PERT restart file path:
    dapertrst_fname = f"{input_dir}/DAPERT/{startdate.year:04d}{startdate.month:02d}" + \
          f"/LIS_DAPERT_{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}0000.d01.bin"
    print(f"[INFO] DA PERT restart file: {dapertrst_fname}")

    # Now edit the target ldt.config with these customized settings
    with open(lisconfig_lsm_target, "rt", encoding='ascii') as file_obj:
        data = file_obj.read()
    data = data.replace("STARTYR", str(startdate.year))
    data = data.replace("STARTMO", str(startdate.month))
    data = data.replace("STARTDA", str(startdate.day))
    data = data.replace("FINALYR", str(finaldate.year))
    data = data.replace("FINALMO", str(finaldate.month))
    data = data.replace("FINALDA", str(finaldate.day))
    data = data.replace("LSMLISLOGFILE", lsm_logfile)
    pertmode = config["lisdarun"]["pertmode"]
    data = data.replace("PERTMODE", str(pertmode))
    if pertmode == "restart":
        data = data.replace("DAPERTRSTFILE", dapertrst_fname)
    data = data.replace("NOAHMP401RSTFILE", noahmprst_fname)
    data = data.replace("HYMAP2RSTFILE", hymaprst_fname)
    with open(lisconfig_lsm_target, "wt", encoding='ascii') as file_obj:
        file_obj.write(data)

def _recursive_chmod(path, mode):
    """Recursively runs chmod"""
    os.chmod(path, mode)
    for dirpath, dirnames, filenames in os.walk(path):
        for dirname in dirnames:
            os.chmod(os.path.join(dirpath, dirname), mode)
        for filename in filenames:
            os.chmod(os.path.join(dirpath, filename), mode)

def _driver():
    """Main driver"""

    configfile, startdate, finaldate = _read_cmd_args()
    config = _read_config(configfile)
    print(f"[INFO] startdate: {startdate}")
    print(f"[INFO] finaldate: {finaldate}")
    # Make directory path (if doesn't exist) for
    # LIS config files that are generated by this script:
    lisconfig_outdir = config["lisdarun"]["lisconfig_outdir"]
    if not os.path.exists(lisconfig_outdir):
        os.makedirs(lisconfig_outdir)

    lisconfig_lsm_target = \
        _create_lisconfig_lsm_target(lisconfig_outdir, startdate)

    _customize_lis_config(config,
                          lisconfig_lsm_target,
                          startdate, finaldate)

    _recursive_chmod(lisconfig_outdir, 0o755)
    print("[INFO] Done generating LIS config file and running LIS DA run.")

if __name__ == "__main__":
    _driver()
