#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: main_forecast.py
#
# PURPOSE: Runs all of the NMME based LIS-HYDRO S2S scripts for NRT
# Hydrological Forecasts. Based on Main_forecast.sh.
#
# REVISION HISTORY:
# 31 Oct 2021: Ryan Zamora, first version
#
#------------------------------------------------------------------------------
"""

#
# Standard modules
#

import subprocess
import sys
from datetime import date

#
# Local constants
#

# In an operational setting these variables will be set based
# on the current month and year. For testing, you can set these
# variables manually.
#
today = date.today()
# MONTH_FORMAL_ABBR = today.strftime("%b") # e.g. Oct
# MONTH_NUM = today.strftime("%m") # e.g. 10
# CURRENT_YEAR = today.strftime("%Y") #e.g. 2021
MONTH_FORMAL_ABBR = "Oct"
MONTH_NUM = 10
CURRENT_YEAR = 2021

MONTH_ABBR = MONTH_FORMAL_ABBR.lower() # e.g. oct

# Current Years
SYR = CURRENT_YEAR
EYR = CURRENT_YEAR

# Climatology Years
CLIM_SYR = 2008
CLIM_EYR = 2020

# Domain extents (from LL to UR):
LAT1 = -40
LAT2 = 40
LON1 = -20
LON2 = 60

# Number of lead months
LEAD_MONTHS = 9

# Number of ensembles in the raw forecast
ENS_NUM = 12

# Forecast Dataset Names
FCST_DATA_TYPE="CFSv2"
NMME_DATA_TYPE="nmme"
NMME_MODELS=["CFSv2", "GEOSv2", "CCSM4", "CCM4", "GNEMO", "GFDL"]

# Config file that sets up directories
CONFIG_FILE="/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/"\
    "scripts/code_library/supplementary_files/bcsd_preproc.cfg"

def run_cmd(cmd, error_msg):
    """Handle running shell command and checking error."""
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print(error_msg)
        sys.exit(1)

def setup_environment():
    """Setup running environment."""
    error_environment = "[ERR] Problem setting up running environment"
    run_cmd("source /usr/share/modules/init/sh", error_environment)
    run_cmd("ulimit -s unlimited", error_environment)
    run_cmd("module load python/GEOSpyD/Ana2019.10_py3.7", error_environment)
    run_cmd("module load cdo/1.9.9rc2", error_environment)
    run_cmd("module load nco", error_environment)

def driver():
    """Main driver."""
    setup_environment()
    print(f"MONTH_ABBR: {MONTH_ABBR}")
    print(f"MONTH_NUM: {MONTH_NUM}")
    print(f"CURRENT_YEAR: {CURRENT_YEAR}")
    #
    ##-----PART A: DOWNLOAD, RESCALE, AND ORGANIZE-----##
    #
    # This section involves downloading the correct forecasts, rescaling them to
    # a unified 25 KM spatial grid, and organizing them in a format for the
    # proceeding tasks
    #

    # Task 1: Generate and rescale 6-hourly files to 25 KM
    #cmd_task_01 = f"python forecast_task_01.py {CURRENT_YEAR} {CURRENT_YEAR} {MONTH_ABBR} {CONFIG_FILE}"
    #run_cmd(cmd_task_01, "[ERR] Problem running forecast_task_01.py")

    # Task 2: Download NMME Data
    # This task has been temporarily removed as the downloading of all raw data
    # will be handled prior to this workflow
    #csh FORECAST_TASK_02.csh $CURRENT_YEAR $MONTH_FORMAL_ABBR $MONTH_NUM

    # Task 3: Rescale and reorganize NMME Data
    #cmd_task_03 = f"python forecast_task_03.py {MONTH_NUM} {CURRENT_YEAR} {CONFIG_FILE}"
    #run_cmd(cmd_task_03, "[ERR] Problem running forecast_task_03.py")

    #
    ##-----PART B: BIAS CORRECTION-----##
    #
    # This sections involves creating monthly bias corrected data for CFSv2
    # (non-precip) and NMME (precip) datasets
    #

    # Task 4: Monthly "BC" step applied to CFSv2
    #cmd_task_04 = f"python forecast_task_04.py {SYR} {EYR} {CLIM_SYR} {CLIM_EYR} {MONTH_ABBR} {MONTH_NUM} {LAT1} {LAT2} {LON1} {LON2} {LEAD_MONTHS} {ENS_NUM} {CONFIG_FILE}"
    #run_cmd(cmd_task_04, "[ERR] Problem running forecast_task_04.py")

    # Task 5: Monthly "BC" step applied to NMME
    #for mod in NMME_MODELS:
    #    cmd_task_05 = f"python forecast_task_05.py {SYR} {EYR} {CLIM_SYR} {CLIM_EYR} {MONTH_ABBR} {MONTH_NUM} {LAT1} {LAT2} {LON1} {LON2} {mod} {LEAD_MONTHS} {CONFIG_FILE}"
    #    run_cmd(cmd_task_05, "[ERR] Problem running forecast_task_05.py")

    ##-----PART C: TEMPORAL DISAGGREGATION-----##
    #
    # This section involves temporally disaggregating the monthly bias corrected
    # forecasts generated in PART B to sub-daily (6-hourly) resolution
    #

    # Task 6: CFSv2 Temporal Disaggregation
    #cmd_task_06 = f"python forecast_task_06.py {SYR} {EYR} {MONTH_ABBR} {MONTH_NUM} {LAT1} {LAT2} {LON1} {LON2} {FCST_DATA_TYPE} {LEAD_MONTHS} {ENS_NUM} {CONFIG_FILE}"
    #run_cmd(cmd_task_06, "[ERR] Problem running forecast_task_06.py")

    #
    # Task 7: Generate symbolic links to sub-daily CFSv2 BC forecasts for NMME
    # temporal disaggregation due to an uneven number of ensembles
    #cmd_task_07 = f"python forecast_task_07.py {CURRENT_YEAR} {MONTH_ABBR} {CONFIG_FILE}"
    #run_cmd(cmd_task_07, "[ERR] Problem running forecast_task_07.py")

    #
    # Task 8: NMME Temporal Disaggregation
    #for mod in NMME_MODELS:
    #    cmd_task_08 = f"python forecast_task_08.py {SYR} {EYR} {MONTH_ABBR} {MONTH_NUM} {LAT1} {LAT2} {LON1} {LON2} {mod} {LEAD_MONTHS} {CONFIG_FILE}"
    #    run_cmd(cmd_task_08, "[ERR] Problem running forecast_task_08.py")

    ##-----PART D: LIS PREPARATION-----##
    #
    # This section involves preparing the bias corrected forecasts for LIS
    #

    # Task 9: Combine the CFSv2 forcing fields into final format for LIS to read
    #cmd_task_09 = f"python forecast_task_09.py {SYR} {EYR} {MONTH_ABBR} {MONTH_NUM} {LAT1} {LAT2} {LON1} {LON2} {FCST_DATA_TYPE} {LEAD_MONTHS} {ENS_NUM} {CONFIG_FILE}"
    #run_cmd(cmd_task_09, "[ERR] Problem running forecast_task_09.py")

    # Task 10: Combine the NMME forcing fields into final format for LIS to read
    # and symbolically link to the reusable CFSv2 met forcings
    #for mod in NMME_MODELS:
    #    cmd_task_10 = f"python forecast_task_10.py {MONTH_ABBR} {MONTH_NUM} {CURRENT_YEAR} {mod} {CONFIG_FILE}"
    #    run_cmd(cmd_task_10, "[ERR] Problem running forecast_task_10.py")


    # Task 11: Copy 9th forecast lead file as 10th forecast lead for LIS runs
    #cmd_task_11 = f"python forecast_task_11.py {MONTH_ABBR} {MONTH_NUM} {CURRENT_YEAR} {LEAD_MONTHS} {CONFIG_FILE}"
    #run_cmd(cmd_task_11, "[ERR] Problem running forecast_task_11.py")


    # Task 12: Temporary task to introduce an all-zero variable V10M due to the
    # way wind is handled in the USAF forcing
    # cmd_task_12 = f"python forecast_task_12.py {MONTH_ABBR} {MONTH_NUM} {CURRENT_YEAR} {ENS_NUM} {LEAD_MONTHS} {CONFIG_FILE}"
    # run_cmd(cmd_task_12, "[ERR] Problem running forecast_task_12.py")


    ##--------------Check Final Preprocessed Files and Folders----------------##
    # A final processing check is available to the user to check that the number
    # of files and file size generated within this workflow meet expectations
    #for mod in NMME_MODELS:
    #    cmd_task_check = f"python check_preprocess_forecast_files.py {MONTH_ABBR} {CURRENT_YEAR} {mod} {LEAD_MONTHS} {CONFIG_FILE}"
    #    run_cmd(cmd_task_check, "[ERR] Problem running check_preprocess_forecast_files.py")

    ##-------------- Completed Preprocessing -----------------##

#
# Main Method
#
if __name__ == "__main__":
    driver()
