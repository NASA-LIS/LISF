#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_05.py
#
# PURPOSE: Computes the bias correction for the NMME dataset. Based on
# FORECAST_TASK_03.sh.
#
# REVISION HISTORY:
# 24 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""

#
# Standard modules
#

import os
import subprocess
import sys

#
# Local constants.  FIXME:  Put in single location for whole system
#

# Path of the main project directory
PROJDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'

# Path of the directory where all the BC codes are kept:
SRCDIR="{}/scripts/code_library".format(PROJDIR)

# Path for where observational files are located:
FORCEDIR="{}/data/forecast".format(PROJDIR)
OBS_CLIM_INDIR="{}/USAF-LIS7.3rc8_25km/raw/Climatology".format(FORCEDIR)

# Mask file
MASK_FILE="{}/supplementary_files/Mask_nafpa.nc".format(SRCDIR)

#  Log file output directory
LOGDIR="{}/scripts/log_files".format(PROJDIR)

#  Calculate bias correction for different variables separately:
OBS_VAR="Rainf_f_tavg"
FCST_VAR="PRECTOT"
UNIT="kg/m^2/s"
VAR_TYPE="PRCP"

#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} FCST_SYR FCST_EYR CLIM_SYR CLIM_EYR month_abbr"\
          "iMonNo lat1 lat2 lon1 lon2 NMME_MODEL lead_months".format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] FCST_SYR: Start year of forecast")
    print("[INFO] FCST_EYR: End year of forecast")
    print("[INFO] CLIM_SYR: Start year of the climatological period")
    print("[INFO] CLIM_EYR: End year of the climatological period")
    print("[INFO] month_abbr: Abbreviation of the initialization month")
    print("[INFO] month_num: Integer number of the initialization month")
    print("[INFO] lat1: Minimum latitudinal extent")
    print("[INFO] lat2: Maximum latitudinal extent")
    print("[INFO] lon1: Minimum longitudinal extent")
    print("[INFO] lon2: Maximum longitudinal extent")
    print("[INFO] NMME_MODEL: NMME model name")
    print("[INFO] lead_months: Number of lead months")

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 13:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # FCST_SYR
    try:
        FCST_SYR = int(sys.argv[1])
    except ValueError:
        print("[ERR] Invalid argument for FCST_SYR!  Received {}" \
            .format(sys.argv[1]))
        _usage()
        sys.exit(1)
    if FCST_SYR < 0:
        print("[ERR] Invalid argument for FCST_SYR!  Received {}" \
              .format(sys.argv[1]))
        _usage()
        sys.exit(1)

    # FCST_EYR
    try:
        FCST_EYR = int(sys.argv[2])
    except ValueError:
        print("[ERR] Invalid argument for FCST_EYR!  Received {}" \
              .format(sys.argv[2]))
        _usage()
        sys.exit(1)
    if FCST_EYR < 0:
        print("[ERR] Invalid argument for FCST_EYR!  Received {}" \
              .format(sys.argv[2]))
        _usage()
        sys.exit(1)

    # CLIM_SYR
    try:
        CLIM_SYR = int(sys.argv[3])
    except ValueError:
        print("[ERR] Invalid argument for CLIM_SYR!  Received {}" \
            .format(sys.argv[3]))
        _usage()
        sys.exit(1)
    if CLIM_SYR < 0:
        print("[ERR] Invalid argument for CLIM_SYR!  Received {}" \
              .format(sys.argv[3]))
        _usage()
        sys.exit(1)

    # CLIM_EYR
    try:
        CLIM_EYR = int(sys.argv[4])
    except ValueError:
        print("[ERR] Invalid argument for CLIM_EYR!  Received {}" \
              .format(sys.argv[4]))
        _usage()
        sys.exit(1)
    if CLIM_EYR < 0:
        print("[ERR] Invalid argument for CLIM_EYR!  Received {}" \
              .format(sys.argv[4]))
        _usage()
        sys.exit(1)

    # MONTH_ABBR
    MONTH_ABBR = str(sys.argv[5])

    # MONTH_NUM
    try:
        MONTH_NUM = int(sys.argv[6])
    except ValueError:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[6]))
        _usage()
        sys.exit(1)
    if MONTH_NUM < 1:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[6]))
        _usage()
        sys.exit(1)
    if MONTH_NUM > 12:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[6]))
        _usage()
        sys.exit(1)

    # LAT1
    try:
        LAT1 = int(sys.argv[7])
    except ValueError:
        print("[ERR] Invalid argument for LAT1!  Received {}" \
              .format(sys.argv[7]))
        _usage()
        sys.exit(1)

    # LAT2
    try:
        LAT2 = int(sys.argv[8])
    except ValueError:
        print("[ERR] Invalid argument for LAT2!  Received {}" \
              .format(sys.argv[8]))
        _usage()
        sys.exit(1)

    # LON1
    try:
        LON1 = int(sys.argv[9])
    except ValueError:
        print("[ERR] Invalid argument for LON1!  Received {}" \
              .format(sys.argv[9]))
        _usage()
        sys.exit(1)

    # LON2
    try:
        LON2 = int(sys.argv[10])
    except ValueError:
        print("[ERR] Invalid argument for LON2!  Received {}" \
              .format(sys.argv[10]))
        _usage()
        sys.exit(1)

    # MONTH_ABBR
    NMME_MODEL = str(sys.argv[11])

    # LEAD_MONTHS
    try:
        LEAD_MONTHS = int(sys.argv[12])
    except ValueError:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[12]))
        _usage()
        sys.exit(1)
    if LEAD_MONTHS < 0:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[12]))
        _usage()
        sys.exit(1)

    return FCST_SYR, FCST_EYR, CLIM_SYR, CLIM_EYR, MONTH_ABBR, MONTH_NUM, \
    LAT1, LAT2, LON1, LON2, NMME_MODEL, LEAD_MONTHS

def _gather_ensemble_info(NMME_MODEL):
    """Gathers ensemble information based on NMME model."""

    # Number of ensembles in the forecast (ENS_NUMF)
    # Number of ensembles in the climatology (ENS_NUMC)
    # Ensemble start index (ENS_START)
    # Ensemble end index (ENS_END)
    if NMME_MODEL == "CFSv2":
        ENS_NUMF=24
        ENS_NUMC=12
        ENS_START=1
        ENS_END=24
    elif NMME_MODEL == "GEOSv2":
        ENS_NUMF=10
        ENS_NUMC=4
        ENS_START=25
        ENS_END=34
    elif NMME_MODEL == "CCM4":
        ENS_NUMF=10
        ENS_NUMC=10
        ENS_START=35
        ENS_END=44
    elif NMME_MODEL == "GNEMO":
        ENS_NUMF=10
        ENS_NUMC=10
        ENS_START=45
        ENS_END=54
    elif NMME_MODEL == "CCSM4":
        ENS_NUMF=10
        ENS_NUMC=10
        ENS_START=55
        ENS_END=64
    elif NMME_MODEL == "GFDL":
        ENS_NUMF=30
        ENS_NUMC=15
        ENS_START=65
        ENS_END=94
    else:
        print("[ERR] Invalid argument for NMME_MODEL!  Received {}" \
            .format(NMME_MODEL))
        sys.exit(1)

    return ENS_NUMF, ENS_NUMC, ENS_START, ENS_END

def _driver():
    """Main driver."""
    FCST_SYR, FCST_EYR, CLIM_SYR, CLIM_EYR, MONTH_ABBR, MONTH_NUM, \
        LAT1, LAT2, LON1, LON2, NMME_MODEL, LEAD_MONTHS = _read_cmd_args()

    # Path for where nmme forecast files are located:
    FCST_CLIM_INDIR="{forcedir}/NMME/raw/Climatology/{month_abbr}01".format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)
    FCST_INDIR="{forcedir}/NMME/raw/Monthly/{month_abbr}01".format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)

    # Path for where output BC forecast file are located:
    OUTDIR="{forcedir}/NMME/bcsd/Monthly/{month_abbr}01".format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    print("[INFO] Processing forecast bias correction of NMME-{} precip" \
        .format(NMME_MODEL))

    ENS_NUMF, ENS_NUMC, ENS_START, ENS_END = _gather_ensemble_info(NMME_MODEL)

    for YEAR in range(FCST_SYR, (FCST_EYR + 1)):
        cmd = "sbatch"
        cmd += " {srcdir}/run_NMME_BCSD_calctest.scr".format(srcdir=SRCDIR)
        cmd += " {srcdir}".format(srcdir=SRCDIR)
        cmd += " {obs_var}".format(obs_var=OBS_VAR)
        cmd += " {fcst_var}".format(fcst_var=FCST_VAR)
        cmd += " {month_num}".format(month_num=MONTH_NUM)
        cmd += " {var_type}".format(var_type=VAR_TYPE)
        cmd += " {unit}".format(unit=UNIT)
        cmd += " {lat1}".format(lat1=LAT1)
        cmd += " {lat2}".format(lat2=LAT2)
        cmd += " {lon1}".format(lon1=LON1)
        cmd += " {lon2}".format(lon2=LON2)
        cmd += " {ens_numc}".format(ens_numc=ENS_NUMC)
        cmd += " {ens_numf}".format(ens_numf=ENS_NUMF)
        cmd += " {nmme_model}".format(nmme_model=NMME_MODEL)
        cmd += " {lead_months}".format(lead_months=LEAD_MONTHS)
        cmd += " {fcst_syr}".format(fcst_syr=YEAR)
        cmd += " {fcst_eyr}".format(fcst_eyr=YEAR)
        cmd += " {clim_syr}".format(clim_syr=CLIM_SYR)
        cmd += " {clim_eyr}".format(clim_eyr=CLIM_EYR)
        cmd += " {mask_file}".format(mask_file=MASK_FILE)
        cmd += " {fcst_clim_indir}".format(fcst_clim_indir=FCST_CLIM_INDIR)
        cmd += " {obs_clim_indir}".format(obs_clim_indir=OBS_CLIM_INDIR)
        cmd += " {fcst_indir}".format(fcst_indir=FCST_INDIR)
        cmd += " {outdir}".format(outdir=OUTDIR)
        cmd += " {logdir}".format(logdir=LOGDIR)
        cmd += " {ens_start}".format(ens_start=ENS_START)
        cmd += " {ens_end}".format(ens_end=ENS_END)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)

    print("[INFO] Completed processing NMME bias correction for: {}"\
        .format(MONTH_ABBR))

#
# Main Method
#
if __name__ == "__main__":
    _driver()
