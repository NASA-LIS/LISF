#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_06.py
#
# PURPOSE: Generate bias-corrected 6-hourly forecasts using raw monthly
# forecasts, bias-corrected monthly forecasts and raw 6-hourly forecasts. Based
# on FORECAST_TASK_06.sh.
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
FORCEDIR="{}/data/forecast/CFSv2_25km".format(PROJDIR)

# Mask file
MASK_FILE_PRECIP="{}/supplementary_files/Mask_nafpa.nc".format(SRCDIR)
MASK_FILE_NONPRECIP="{}/supplementary_files/Mask_nafpa.nc".format(SRCDIR)


#  Calculate bias correction for different variables separately:
OBS_VAR_LIST=["LWGAB", "SWGDN", "PS", "QV2M", "T2M", "U10M"]
FCST_VAR_LIST=["LWS", "SLRSF", "PS", "Q2M", "T2M", "WIND10M"]
UNIT_LIST=["W/m^2", "W/m^2", "Pa", "kg/kg", "K", "m/s"]

#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} FCST_SYR FCST_EYR MONTH_ABBR MONTH_NUM LAT1 LAT2 LON1"\
        "LON2 MODEL_NAME LEAD_MONTHS ENS_NUM".format(sys.argv[0])
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
    print("[INFO] MODEL_NAME: Model name (should be CFSv2)")
    print("[INFO] LEAD_MONTHS: Number of lead months")
    print("[INFO] ENS_NUM: Integer number of ensembles")

def _read_cmd_args():
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

    # MODEL_NAME
    MODEL_NAME = str(sys.argv[9])

    # LEAD_MONTHS
    try:
        LEAD_MONTHS = int(sys.argv[10])
    except ValueError:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[10]))
        _usage()
        sys.exit(1)
    if LEAD_MONTHS < 0:
        print("[ERR] Invalid argument for LEAD_MONTHS!  Received {}" \
              .format(sys.argv[10]))
        _usage()
        sys.exit(1)

    # ENS_NUM
    try:
        ENS_NUM = int(sys.argv[11])
    except ValueError:
        print("[ERR] Invalid argument for ENS_NUM!  Received {}" \
              .format(sys.argv[11]))
        _usage()
        sys.exit(1)
    if ENS_NUM < 0:
        print("[ERR] Invalid argument for ENS_NUM!  Received {}" \
              .format(sys.argv[11]))
        _usage()
        sys.exit(1)

    return FCST_SYR, FCST_EYR, MONTH_ABBR, MONTH_NUM, LAT1, LAT2, LON1, LON2, \
    MODEL_NAME, LEAD_MONTHS, ENS_NUM

def _driver():
    """Main driver."""
    FCST_SYR, FCST_EYR, MONTH_ABBR, MONTH_NUM, LAT1, LAT2, LON1, LON2, \
    MODEL_NAME, LEAD_MONTHS, ENS_NUM = _read_cmd_args()

    # Path for where forecast and bias corrected files are located:
    SUBDAILY_RAW_FCST_DIR="{forcedir}/raw/6-Hourly/{month_abbr}01" \
        .format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)
    MONTHLY_RAW_FCST_DIR="{forcedir}/raw/Monthly/{month_abbr}01" \
        .format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)
    MONTHLY_BC_FCST_DIR="{forcedir}/bcsd/Monthly/{month_abbr}01" \
        .format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)

    OUTDIR="{forcedir}/bcsd/6-Hourly/{month_abbr}01" \
        .format(forcedir=FORCEDIR, month_abbr=MONTH_ABBR)

    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    print("[INFO] Processing temporal disaggregation of CFSv2 variables")
    for YEAR in range(FCST_SYR, (FCST_EYR + 1)):
        for VAR_NUM in range(len(OBS_VAR_LIST)):
            if VAR_NUM == 1:
                VAR_TYPE="PRCP"
            else:
                VAR_TYPE="TEMP"

            OBS_VAR=OBS_VAR_LIST[VAR_NUM]
            FCST_VAR=FCST_VAR_LIST[VAR_NUM]
            UNIT=UNIT_LIST[VAR_NUM]

            cmd = "sbatch"
            cmd += " {srcdir}/run_Temporal_disagg.scr".format(srcdir=SRCDIR)
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
            cmd += " {model_name}".format(model_name=MODEL_NAME)
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
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem calling sbatch!")
                sys.exit(1)

    print("[INFO] Completed CFSv2 temporal disaggregation for: {}"\
        .format(MONTH_ABBR))

#
# Main Method
#
if __name__ == "__main__":
    _driver()
