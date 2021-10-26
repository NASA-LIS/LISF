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
#
#------------------------------------------------------------------------------
"""

# Standard modules
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
_NMME_MODELS = ["CCM4", "CCSM4", "CFSv2", "GEOSv2", "GFDL", "GNEMO"]

# Ensemble members per NMME model
_ENSEMBLE_SIZES = {
    "CCM4" : 10,
    "CCSM4" : 10,
    "GNEMO" : 10,
    "GEOSv2" : 10,
    "CFSv2" : 24,
    "GFDL" : 30,
}

# Options for scaling to 0.25 deg
_NMME_SCALINGS = {
    "CCM4" : "downscale",
    "CCSM4" : "downscale",
    "GNEMO" : "downscale",
    "GEOSv2" : "downscale",
    "CFSv2" : "upscale",
    "GFDL" : "upscale",
}

# NOTE: Make a master file with all of these paths! FIXME
_WORK_DIR = \
    "/discover/nobackup/projects/lis_aist17/karsenau/GHI_S2S/AFRICOM/LDT_ICs"
_INPUT_DIR = "../SPINUP/Run-E/retro_2densgrid_hymap1mem/"

_WORK_DIR = \
    "/discover/nobackup/projects/lis_aist17/emkemp/AFWA/lis74_s2s_patches/work"

_INPUT_DIR = "./SPINUP/Run-E/retro_2densgrid_hymap1mem/"
_LDT_INPUT_FILE = "./input/lis_input.s2s_africom.%s_%s.25km.nc" \
    %(_LSM_NAME, _ROUTING_NAME)
_CONFIGS_OUTPUT_DIR = "./ldt.config_files"
_TEMPLATE_DIR = "./template_files"
_LDTCONFIG_LSM_TEMPLATE = \
    "%s/ldt.config_%s_nmme_TEMPLATE" %(_TEMPLATE_DIR, _LSM_NAME)

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
    print("[INFO] Current year / month: %4.4d / %2.2d" \
          %(currentdate.year, currentdate.month))
    rstdate = datetime.date(currentdate.year, currentdate.month, 1)
    print("[INFO] Rst Start Month Name: %s" %(rstdate.strftime("%b")))
    prevdate = rstdate - datetime.timedelta(days=1)
    print("[INFO] Rst year / month / day: %4.4d / %2.2d / %2.2d" \
          %(prevdate.year, prevdate.month, prevdate.day))
    return currentdate, prevdate

def _create_ldtconfig_lsm_target(nmme_model, currentdate):
    """Create name of new ldt.config file."""
    lc_nmmemodel = nmme_model.lower()
    ldtconfig_nameconv_lsm = "ldt.config_%s_nmme_%s" %(_LSM_NAME,
                                                       lc_nmmemodel)
    print("[INFO] Processing LDT ensemble restart file for %s" \
          %(ldtconfig_nameconv_lsm))
    ldtconfig_lsm_target = "%s/%s_%4.4d%2.2d" %(_CONFIGS_OUTPUT_DIR,
                                                ldtconfig_nameconv_lsm,
                                                currentdate.year,
                                                currentdate.month)
    return ldtconfig_lsm_target

def _customize_ldt_config(ldtconfig_lsm_target, lsm_rstdir, currentdate,
                          prevdate, nmme_model):
    """Customize new ldt.config file."""
    shutil.copy(_LDTCONFIG_LSM_TEMPLATE, ldtconfig_lsm_target)

    num_ensmems = _ENSEMBLE_SIZES[nmme_model]
    ldt_rstgen = _NMME_SCALINGS[nmme_model]

    # Build customized settings for target ldt.config file
    input_fname = "%s/LIS_RST_NOAHMP401_%4.4d%2.2d%2.2d2345.d01.nc" \
        %(lsm_rstdir, prevdate.year, prevdate.month, prevdate.day)

    rst_date = "%4.4d%2.2d%2.2d" \
        %(prevdate.year, prevdate.month, prevdate.day)
    rst_monname = currentdate.strftime("%b")

    output_fname = "%s/LIS_RST_NOAHMP401_%s2345.ICS_%s%4.4d.ens%d.nc" \
        %(nmme_model, rst_date, rst_monname, prevdate.year, num_ensmems)

    lsm_logfile = "%s/ldtlog_%s_%s%4.4d" \
        %(nmme_model, _LSM_NAME, rst_monname, prevdate.year)

    mask_parmlogfile = "%s/MaskParamFill.log" %(nmme_model)

    # Now edit the target ldt.config with these customized settings
    with open(ldtconfig_lsm_target, "rt") as file_obj:
        data = file_obj.read()
    data = data.replace("LDTINPUTFILE", _LDT_INPUT_FILE)
    data = data.replace("LDTRSTGENOPT", ldt_rstgen)
    data = data.replace("INPUTRSTFILE", input_fname)
    data = data.replace("OUTPUTRSTFILE", "./%s" %(output_fname))
    data = data.replace("OUTENSMEMS", "%s" %(num_ensmems))
    data = data.replace("LSMLDTLOGFILE", "./%s" %(lsm_logfile))
    data = data.replace("PARAMLOGFILE", "./%s" %(mask_parmlogfile))
    data = data.replace("MODELDIR", "./%s/" %(nmme_model))
    with open(ldtconfig_lsm_target, "wt") as file_obj:
        file_obj.write(data)

def _driver():
    """Main driver"""

    print("[INFO] Working directory is: %s" %(_WORK_DIR))
    os.chdir(_WORK_DIR)

    if not os.path.exists(_LDT_EXEC):
        print("[ERR] %s does not exist!" %(_LDT_EXEC))
        sys.exit(1)

    if not os.path.exists(_CONFIGS_OUTPUT_DIR):
        os.makedirs(_CONFIGS_OUTPUT_DIR)

    currentdate, prevdate = _handle_dates()

    lsm_rstdir = "%s/SURFACEMODEL/%4.4d%2.2d" %(_INPUT_DIR,
                                                prevdate.year,
                                                prevdate.month)
    hymap_rstdir = "%s/ROUTING/%4.4d%2.2d" %(_INPUT_DIR,
                                             prevdate.year,
                                             prevdate.month)

    for nmme_model in _NMME_MODELS:

        if not os.path.exists(nmme_model):
            os.makedirs(nmme_model)

        num_ensmems = _ENSEMBLE_SIZES[nmme_model]
        print("[INFO] Number of members: %s for model %s" %(num_ensmems,
                                                            nmme_model))
        ldt_rstgen = _NMME_SCALINGS[nmme_model]
        print("[INFO] LDT ensemble restart generation: %s for model %s" \
              %(ldt_rstgen, nmme_model))

        ldtconfig_lsm_target = \
            _create_ldtconfig_lsm_target(nmme_model, currentdate)

        _customize_ldt_config(ldtconfig_lsm_target, lsm_rstdir,
                              currentdate, prevdate,
                              nmme_model)

        # Run LDT
        cmd = "%s %s" %(_LDT_EXEC, ldtconfig_lsm_target)
        print("[INFO] %s" %(cmd))
        subprocess.run(cmd, shell=True, check=True)

        # Copy HyMAP2 restart file and rename, placing 1-member restart
        # file per NMME model directory
        rst_date = "%4.4d%2.2d%2.2d" \
            %(prevdate.year, prevdate.month, prevdate.day)
        rst_monname = currentdate.strftime("%b")

        hymap_inrstfile = "%s/LIS_RST_HYMAP2_router_%s2345.d01.nc" \
            %(hymap_rstdir, rst_date)
        hymap_outrstfile = \
            "./%s/LIS_RST_HYMAP2_router_%s2345.ICS_%s%4.4d.ens1.nc" \
            %(nmme_model, rst_date, rst_monname, prevdate.year)
        shutil.copy(hymap_inrstfile, hymap_outrstfile)

        # Recursively update permissions in work subdirectory
        _recursive_chmod(nmme_model, 0o755)

    _recursive_chmod(_CONFIGS_OUTPUT_DIR, 0o755)
    print("[INFO] Done generating LDT config files and restart files.")

if __name__ == "__main__":
    _driver()
