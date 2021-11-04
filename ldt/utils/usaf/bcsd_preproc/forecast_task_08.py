#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_08.py
#
# PURPOSE: Generate bias-corrected 6-hourly nmme forecasts using raw monthly
# forecasts, bias-corrected monthly forecasts and raw 6-hourly forecasts. Based
# on FORECAST_TASK_08.sh.
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

#  Log file output directory
LOGDIR="{}/scripts/log_files".format(PROJDIR)

# Path for where forecast files are located:
FORCEDIR="{}/data/forecast/NMME".format(PROJDIR)

# Mask file
MASK_FILE_PRECIP="{}/supplementary_files/Mask_nafpa.nc".format(SRCDIR)
MASK_FILE_NONPRECIP="{}/supplementary_files/Mask_nafpa.nc".format(SRCDIR)

#  Calculate bias correction for different variables separately:
OBS_VAR="PRECTOT"
FCST_VAR="PRECTOT"
UNIT="kg/m^2/s"
VAR_TYPE='PRCP'

#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} FCST_SYR FCST_EYR MONTH_ABBR MONTH_NUM LAT1 LAT2 "\
        "LON1 LON2 FCST_TYPE NMME_MODEL LEAD_MONTHS".format(sys.argv[0])
    print(txt)
    print("[INFO] where")
    print("[INFO] FCST_SYR: Start year of forecast")
    print("[INFO] FCST_EYR: End year of forecast")
    print("[INFO] MONTH_ABBR: Abbreviation of the initialization month")
    print("[INFO] MONTH_NUM: Integer number of the initialization month")
    print("[INFO] LAT1: Minimum latitudinal extent")
    print("[INFO] LAT2: Maximum latitudinal extent")
    print("[INFO] LON1: Minimum longitudinal extent")
    print("[INFO] LON2: Maximum longitudinal extent")
    print("[INFO] FCST_TYPE: Forecast type (should be nmme)")
    print("[INFO] NMME_MODEL: NMME model name")
    print("[INFO] LEAD_MONTHS: Number of lead months")

def read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 12:
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

    # MONTH_ABBR
    MONTH_ABBR = str(sys.argv[3])

    # MONTH_NUM
    try:
        MONTH_NUM = int(sys.argv[4])
    except ValueError:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[4]))
        _usage()
        sys.exit(1)
    if MONTH_NUM < 1:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[4]))
        _usage()
        sys.exit(1)
    if MONTH_NUM > 12:
        print("[ERR] Invalid argument for MONTH_NUM!  Received {}" \
              .format(sys.argv[4]))
        _usage()
        sys.exit(1)

    # LAT1
    try:
        LAT1 = int(sys.argv[5])
    except ValueError:
        print("[ERR] Invalid argument for LAT1!  Received {}" \
              .format(sys.argv[5]))
        _usage()
        sys.exit(1)

    # LAT2
    try:
        LAT2 = int(sys.argv[6])
    except ValueError:
        print("[ERR] Invalid argument for LAT2!  Received {}" \
              .format(sys.argv[6]))
        _usage()
        sys.exit(1)

    # LON1
    try:
        LON1 = int(sys.argv[7])
    except ValueError:
        print("[ERR] Invalid argument for LON1!  Received {}" \
              .format(sys.argv[7]))
        _usage()
        sys.exit(1)

    # LON2
    try:
        LON2 = int(sys.argv[8])
    except ValueError:
        print("[ERR] Invalid argument for LON2!  Received {}" \
              .format(sys.argv[8]))
        _usage()
        sys.exit(1)

    # FCST_TYPE
    FCST_TYPE = str(sys.argv[9])

    # NMME_MODEL
    NMME_MODEL = str(sys.argv[10])

    # LEAD_MONTHS
    try:
        LEAD_MONTHS = int(sys.argv[11])
    except ValueError:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[11]))
        _usage()
        sys.exit(1)
    if LEAD_MONTHS < 0:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[11]))
        _usage()
        sys.exit(1)

    return FCST_SYR, FCST_EYR, MONTH_ABBR, MONTH_NUM, LAT1, LAT2, LON1, LON2, \
        FCST_TYPE, NMME_MODEL, LEAD_MONTHS

def gather_ensemble_info(NMME_MODEL):
    """Gathers ensemble information based on NMME model."""

    # Number of ensembles in the forecast (ENS_NUM)
    # Ensemble start index (ENS_START)
    # Ensemble end index (ENS_END)
    if NMME_MODEL == "CFSv2":
        ENS_NUM=24
        ENS_START=1
        ENS_END=24
    elif NMME_MODEL == "GEOSv2":
        ENS_NUM=10
        ENS_START=25
        ENS_END=34
    elif NMME_MODEL == "CCM4":
        ENS_NUM=10
        ENS_START=35
        ENS_END=44
    elif NMME_MODEL == "GNEMO":
        ENS_NUM=10
        ENS_START=45
        ENS_END=54
    elif NMME_MODEL == "CCSM4":
        ENS_NUM=10
        ENS_START=55
        ENS_END=64
    elif NMME_MODEL == "GFDL":
        ENS_NUM=30
        ENS_START=65
        ENS_END=94
    else:
        print("[ERR] Invalid argument for NMME_MODEL!  Received {}" \
            .format(NMME_MODEL))
        sys.exit(1)

    return ENS_NUM, ENS_START, ENS_END

def _driver():
    """Main driver."""
    FCST_SYR, FCST_EYR, MONTH_ABBR, MONTH_NUM, LAT1, LAT2, LON1, LON2, \
    FCST_TYPE, NMME_MODEL, LEAD_MONTHS = read_cmd_args()

    # Path for where forecast and bias corrected files are located:
    SUBDAILY_RAW_FCST_DIR="{forcedir}/linked_cfsv2_precip_files/{month_abbr}01"\
        .format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)
    MONTHLY_RAW_FCST_DIR="{forcedir}/raw/Monthly/{month_abbr}01" \
        .format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)
    MONTHLY_BC_FCST_DIR="{forcedir}/bcsd/Monthly/{month_abbr}01" \
        .format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)

    OUTDIR="{forcedir}/bcsd/6-Hourly/{month_abbr}01/{nmme_model}" \
        .format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR, nmme_model=NMME_MODEL)

    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    ENS_NUM, ENS_START, ENS_END = gather_ensemble_info(NMME_MODEL)
    print("[INFO] Processing temporal disaggregation of CFSv2 variables")
    for YEAR in range(FCST_SYR, (FCST_EYR + 1)):
        cmd = "sbatch"
        cmd += " {srcdir}/run_NMME_Temporal_disagg.scr".format(srcdir=SRCDIR)
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
        cmd += " {nmme_model}".format(nmme_model=NMME_MODEL)
        cmd += " {ens_num}".format(ens_num=ENS_NUM)
        cmd += " {lead_months}".format(lead_months=LEAD_MONTHS)
        cmd += " {year}".format(year=YEAR)
        cmd += " {year}".format(year=YEAR)
        cmd += " {mask_file_p}".format(mask_file_p=MASK_FILE_PRECIP)
        cmd += " {mask_file_np}".format(mask_file_np=MASK_FILE_NONPRECIP)
        cmd += " {monthly_bc_fcst_dir}"\
            .format(monthly_bc_fcst_dir=MONTHLY_BC_FCST_DIR)
        cmd += " {monthly_raw_fcst_dir}"\
            .format(monthly_raw_fcst_dir=MONTHLY_RAW_FCST_DIR)
        cmd += " {subdaily_raw_fcst_dir}"\
            .format(subdaily_raw_fcst_dir=SUBDAILY_RAW_FCST_DIR)
        cmd += " {outdir}".format(outdir=OUTDIR)
        cmd += " {logdir}".format(logdir=LOGDIR)
        cmd += " {ens_start}".format(ens_start=ENS_START)
        cmd += " {ens_end}".format(ens_end=ENS_END)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)

    print("[INFO] Completed NMME temporal disaggregation for: {}"\
        .format(MONTH_ABBR))

#
# Main Method
#
if __name__ == "__main__":
    _driver()

