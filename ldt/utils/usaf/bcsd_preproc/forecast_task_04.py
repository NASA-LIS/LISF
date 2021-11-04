#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: forecast_task_04.py
#
# PURPOSE: Computes the bias correction for the forecast (CFSv2) dataset. Based
# on FORECAST_TASK_04.sh.
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

# Path for where observational & forecast files are located:
FORCEDIR="{}/data/forecast".format(PROJDIR)
OBS_INDIR="{}/USAF-LIS7.3rc8_25km".format(FORCEDIR)
FCST_INDIR="{}/CFSv2_25km".format(FORCEDIR)

# Mask file
MASK_FILE="{}/supplementary_files/Mask_nafpa.nc".format(SRCDIR)

#  Log file output directory
LOGDIR="{}/scripts/log_files".format(PROJDIR)

#  Calculate bias correction for different variables separately:
OBS_VAR_LIST=["Rainf_f_tavg", "LWdown_f_tavg", "SWdown_f_tavg", "Psurf_f_tavg", "Qair_f_tavg", "Tair_f_tavg", "Wind_f_tavg"]
FCST_VAR_LIST=["PRECTOT", "LWS", "SLRSF", "PS", "Q2M", "T2M", "WIND10M"]
UNIT_LIST=["kg/m^2/s", "W/m^2", "W/m^2", "Pa", "kg/kg", "K", "m/s"]

#
# Local methods
#

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: {} FCST_SYR FCST_EYR CLIM_SYR CLIM_EYR month_abbr"\
          "iMonNo lat1 lat2 lon1 lon2 lead_months ens_numc".format(sys.argv[0])
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
    print("[INFO] lead_months: Number of lead months")
    print("[INFO] ens_num: Number of ensembles")

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

    # ENS_NUM
    try:
        ENS_NUM = int(sys.argv[12])
    except ValueError:
        print("[ERR] Invalid argument for ENS_NUM!  Received {}" \
              .format(sys.argv[12]))
        _usage()
        sys.exit(1)
    if ENS_NUM < 0:
        print("[ERR] Invalid argument for ENS_NUM!  Received {}" \
              .format(sys.argv[12]))
        _usage()
        sys.exit(1)

    return FCST_SYR, FCST_EYR, CLIM_SYR, CLIM_EYR, MONTH_ABBR, MONTH_NUM, \
    LAT1, LAT2, LON1, LON2, LEAD_MONTHS, ENS_NUM
def _driver():
    """Main driver."""
    FCST_SYR, FCST_EYR, CLIM_SYR, CLIM_EYR, MONTH_ABBR, MONTH_NUM, \
        LAT1, LAT2, LON1, LON2, LEAD_MONTHS, ENS_NUM = _read_cmd_args()

    # BC output directory for FCST:
    OUTDIR="{fcst_indir}/bcsd/Monthly/{month_abbr}01" \
        .format(fcst_indir=FCST_INDIR, month_abbr=MONTH_ABBR)

    print("[INFO] Processing forecast bias correction of CFSv2 variables")

    for VAR_NUM in range(len(OBS_VAR_LIST)):
        if VAR_NUM in [0, 2]:
            VAR_TYPE="PRCP"
        else:
            VAR_TYPE="TEMP"

        OBS_VAR=OBS_VAR_LIST[VAR_NUM]
        FCST_VAR=FCST_VAR_LIST[VAR_NUM]
        UNIT=UNIT_LIST[VAR_NUM]
        print("{var_num} {fcst_var}".format(var_num=VAR_NUM, fcst_var=FCST_VAR))

        cmd = "sbatch"
        cmd += " {srcdir}/run_BCSD_calctest.scr".format(srcdir=SRCDIR)
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
        cmd += " {ens_num}".format(ens_num=ENS_NUM)
        cmd += " {lead_months}".format(lead_months=LEAD_MONTHS)
        cmd += " {fcst_syr}".format(fcst_syr=FCST_SYR)
        cmd += " {fcst_eyr}".format(fcst_eyr=FCST_EYR)
        cmd += " {clim_syr}".format(clim_syr=CLIM_SYR)
        cmd += " {clim_eyr}".format(clim_eyr=CLIM_EYR)
        cmd += " {mask_file}".format(mask_file=MASK_FILE)
        cmd += " {obs_indir}".format(obs_indir=OBS_INDIR)
        cmd += " {fcst_indir}".format(fcst_indir=FCST_INDIR)
        cmd += " {outdir}".format(outdir=OUTDIR)
        cmd += " {logdir}".format(logdir=LOGDIR)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            print("[ERR] Problem calling sbatch!")
            sys.exit(1)

    print("[INFO] Completed processing forecast bias correction for: {}" \
        .format(MONTH_ABBR))

#
# Main Method
#
if __name__ == "__main__":
    _driver()
