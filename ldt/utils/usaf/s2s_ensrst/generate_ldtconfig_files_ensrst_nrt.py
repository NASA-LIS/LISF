#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: generate_ldtconfig_files_ensrst_nrt.py
#
# PURPOSE: Generates ldt.config files for adjusting NoahMP401 ensemble restart
# files for the six different NMME models.  Also copies/renames restart files
# for the 1-member HyMAP2 restart for use as ICs for the forecast runs.  Based
# on generate_ldtconfig_files_ensrst_nrt.sh by Kristi Arsenault/SAIC.
#
# REQUIREMENTS as of 14 Oct 2021:
# * Python 3.8 or higher
#
# REVISION HISTORY:
# * 14 Oct 2021: Eric Kemp/SSAI, first version.
# * 01 Nov 2021: Eric Kemp/SSAI, addressed pylint format complaints.
# * 15 Nov 2021: K. Arsenault/SAIC, made minor fixes and added config entries.
# * 19 Nov 2021: Eric Kemp/SSAI, fixed typos reported by K. Arsenault and
#                R. Zamora.
# * 18 Jan 2022: K. Arsenault/SAIC, fixed prev to current year in naming convention
#
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
import datetime
import os
import shutil
import subprocess
import sys

# Private constants
_LDT_EXEC = "./LDT"
_LSM_NAME = "noahmp401"
_ROUTING_NAME = "hymap2"
_INPUT_NUMFCSTMONS = 9

# NMME model names
_NMME_MODELS = ["CCM4", "CCSM4", "CFSv2", "GEOSv2", "GFDL", "GNEMO5"]

# Ensemble members per NMME model
_ENSEMBLE_SIZES = {
    "CCM4" : 10,
    "CCSM4" : 10,
    "GNEMO5" : 10,
    "GEOSv2" : 10,
    "CFSv2" : 24,
    "GFDL" : 15,
}

# Options for scaling to 0.25 deg
_NMME_SCALINGS = {
    "CCM4" : "downscale",
    "CCSM4" : "downscale",
    "GNEMO5" : "downscale",
    "GEOSv2" : "downscale",
    "CFSv2" : "upscale",
    "GFDL" : "upscale",
}

# Main path entries are located in the following python config file:
_PYCONFIGFILE = "./ldt_ics.cfg"
#
# Other local input files
_LDT_INPUT_FILE = \
    f"./input/lis_input.s2s_africom.{_LSM_NAME}_{_ROUTING_NAME}.25km.nc"

_CONFIGS_OUTPUT_DIR = "./ldt.config_files"
_TEMPLATE_DIR = "./template_files"
_LDTCONFIG_LSM_TEMPLATE = \
    f"{_TEMPLATE_DIR}/ldt.config_{_LSM_NAME}_nmme_TEMPLATE"

def _read_config(pyconfigfile):
    """Read from LIS DA run config file."""
    config = configparser.ConfigParser()
    config.read(pyconfigfile)
    return config

def _recursive_chmod(path, mode):
    """Recursively runs chmod"""
    os.chmod(path, mode)
    for dirpath, dirnames, filenames in os.walk(path):
        for dirname in dirnames:
            os.chmod(os.path.join(dirpath, dirname), mode)
        for filename in filenames:
            os.chmod(os.path.join(dirpath, filename), mode)

def _handle_dates():
    """Collect and return date information."""
    currentdate = datetime.date.today()
    print("[INFO] Current year / month: " + \
          f"{currentdate.year:04d} / {currentdate.month:02d}")
    rstdate = datetime.date(currentdate.year, currentdate.month, 1)
    print(f"[INFO] Rst Start Month Name: {rstdate.strftime('%b')}")
    prevdate = rstdate - datetime.timedelta(days=1)
    print("[INFO] Rst year / month / day: " + \
          f"{prevdate.year:04d} / {prevdate.month:02d} / {prevdate.day:02d}")
    return currentdate, prevdate

def _create_ldtconfig_lsm_target(nmme_model, currentdate):
    """Create name of new ldt.config file."""
    lc_nmmemodel = nmme_model.lower()
    ldtconfig_nameconv_lsm = \
        f"ldt.config_{_LSM_NAME}_nmme_{lc_nmmemodel}"
    print("[INFO] Processing LDT ensemble restart file for " + \
          f"{ldtconfig_nameconv_lsm}")
    ldtconfig_lsm_target = \
        f"{_CONFIGS_OUTPUT_DIR}/{ldtconfig_nameconv_lsm}_" + \
        f"{currentdate.year:04d}{currentdate.month:02d}"
    return ldtconfig_lsm_target

def _customize_ldt_config(ldtconfig_lsm_target, lsm_rstdir, currentdate,
                          prevdate, nmme_model):
    """Customize new ldt.config file."""
    shutil.copy(_LDTCONFIG_LSM_TEMPLATE, ldtconfig_lsm_target)

    num_ensmems = _ENSEMBLE_SIZES[nmme_model]
    ldt_rstgen = _NMME_SCALINGS[nmme_model]

    # Build customized settings for target ldt.config file
    input_fname = f"{lsm_rstdir}/LIS_RST_NOAHMP401_" + \
        f"{prevdate.year:04d}{prevdate.month:02d}{prevdate.day:02d}2345.d01.nc"
    rst_date = f"{prevdate.year:04d}{prevdate.month:02d}{prevdate.day:02d}"
    rst_monname = currentdate.strftime("%b")

    output_fname = f"{nmme_model}/LIS_RST_NOAHMP401_{rst_date}2345.ICS_" + \
        f"{rst_monname}{currentdate.year:04d}.ens{num_ensmems:d}.nc"
    lsm_logfile = \
        f"{nmme_model}/ldtlog_{_LSM_NAME}_{rst_monname}{currentdate.year:04d}"

    mask_parmlogfile = f"{nmme_model}/MaskParamFill.log"

    # Now edit the target ldt.config with these customized settings
    with open(ldtconfig_lsm_target, "rt", encoding='ascii') as file_obj:
        data = file_obj.read()
    data = data.replace("LDTINPUTFILE", _LDT_INPUT_FILE)
    data = data.replace("LDTRSTGENOPT", ldt_rstgen)
    data = data.replace("INPUTRSTFILE", input_fname)
    data = data.replace("OUTPUTRSTFILE", f"./{output_fname}")
    data = data.replace("OUTENSMEMS", f"{num_ensmems:d}")
    data = data.replace("LSMLDTLOGFILE", f"./{lsm_logfile}")
    data = data.replace("PARAMLOGFILE", f"./{mask_parmlogfile}")
    data = data.replace("MODELDIR", f"./{nmme_model}/")
    with open(ldtconfig_lsm_target, "wt", encoding="ascii") as file_obj:
        file_obj.write(data)

def _driver():
    """Main driver"""

    # Check if python script config file exists.
    pyconfigfile = _PYCONFIGFILE
    if not os.path.exists(pyconfigfile):
        print(f"[ERR] {pyconfigfile} does not exist!")
        sys.exit(1)

    config = _read_config(pyconfigfile)

    # print(f"[INFO] Working directory is: {_WORK_DIR}")
    # os.chdir(_WORK_DIR)
    workdir = config["ldtics"]["workdir"]
    print(f"[INFO] Working directory is: {workdir}")
    os.chdir(workdir)

    if not os.path.exists(_LDT_EXEC):
        print(f"[ERR] {_LDT_EXEC} does not exist!")
        sys.exit(1)

    if not os.path.exists(_CONFIGS_OUTPUT_DIR):
        os.makedirs(_CONFIGS_OUTPUT_DIR)

    currentdate, prevdate = _handle_dates()

    # Form input restart directory and filenames together
    inputdir = config["ldtics"]["inputdir"]
    lsm_rstdir = f"{inputdir}/SURFACEMODEL/" + \
        f"{prevdate.year:04d}{prevdate.month:02d}"
    hymap_rstdir = f"{inputdir}/ROUTING/" + \
        f"{prevdate.year:04d}{prevdate.month:02d}"

    for nmme_model in _NMME_MODELS:

        if not os.path.exists(nmme_model):
            os.makedirs(nmme_model)

        num_ensmems = _ENSEMBLE_SIZES[nmme_model]
        print(f"[INFO] Number of members: {num_ensmems:d}" + \
              f" for model {nmme_model}")
        ldt_rstgen = _NMME_SCALINGS[nmme_model]
        print("[INFO] LDT ensemble restart generation: " + \
              f"{ldt_rstgen} for model {nmme_model}")

        ldtconfig_lsm_target = \
            _create_ldtconfig_lsm_target(nmme_model, currentdate)

        _customize_ldt_config(ldtconfig_lsm_target, lsm_rstdir,
                              currentdate, prevdate,
                              nmme_model)

        # Run LDT
        cmd = f"mpirun -np 1 {_LDT_EXEC} {ldtconfig_lsm_target}"
        print(f"[INFO] {cmd}")
        subprocess.run(cmd, shell=True, check=True)

        # Copy HyMAP2 restart file and rename, placing 1-member restart
        # file per NMME model directory
        rst_date = f"{prevdate.year:04d}{prevdate.month:02d}{prevdate.day:02d}"
        rst_monname = currentdate.strftime("%b")

        hymap_inrstfile = \
            f"{hymap_rstdir}/LIS_RST_HYMAP2_router_{rst_date}2345.d01.nc"
        hymap_outrstfile = \
            f"./{nmme_model}/LIS_RST_HYMAP2_router_{rst_date}2345.ICS_" + \
            f"{rst_monname}{currentdate.year:04d}.ens1.nc"
        shutil.copy(hymap_inrstfile, hymap_outrstfile)

        # Recursively update permissions in work subdirectory
        _recursive_chmod(nmme_model, 0o755)

    _recursive_chmod(_CONFIGS_OUTPUT_DIR, 0o755)
    print("[INFO] Done generating LDT config files and restart files.")

if __name__ == "__main__":
    _driver()
