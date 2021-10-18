#!/bin/sh
#
# The following main script runs all the NMME based LIS-HYDRO S2S
#  scripts for NRT Hydrological Forecasts:

#
################~~~~~SETUP ENVIRONMENT~~~~~################
#
##-----INPUT: TIME PERIOD-----##
#
# In an operational setting these variables will be set based
# on the current month and year. For testing, you can set these
# variables manually.
#
#currentmon="Jan"
#currentmon1=01
#currentyear=2021
currentmon=`date +%b` # e.g. Oct
currentmon1=`date +%m` # e.g. 10
currentyear=`date +%Y`  #e.g. 2021

currentmon2=${currentmon,} # e.g. oct

echo ${currentmon}
echo ${currentmon1}
echo ${currentmon2}
echo ${currentyear}

##-----INPUT: SETUP ENVIRONMENT-----##
# The remaining variables set here should not be changed
# under normal circumstances

# Climatology Years
CLIM_SYR=2008
CLIM_EYR=2020

SYR=${currentyear}
EYR=${currentyear}

# Forecast Dataset Names
fcstdatatype="CFSv2"
fcstdatatype1="nmme"
declare -a nmmemodel=('CFSv2' 'GEOSv2' 'CCSM4' 'CCM4' 'GNEMO' 'GFDL')

# Domain extents (from LL to UR):
lat1=-40
lat2=40
lon1=-20
lon2=60

# Number of lead months
leadmon=9

# Number of ensembles in the raw forecast
bec=12

# Setup Environment
source /usr/share/modules/init/sh
ulimit -s unlimited

#
################~~~~~PREPROCESSING~~~~~################
#
##-----PART A: DOWNLOAD, RESCALE, AND ORGANIZE-----##
#
# This section involves downloading the correct forecasts, rescaling them to a unified
# 25 KM spatial grid, and organizing them in a format for the proceeding tasks
#

# Task 1: Generate 6-hourly files - "SD" part; CFSv2 1Deg--> 25KM, and USAF 10KM --> 25KM; netcdf
#sh FORECAST_TASK_01.sh ${currentyear} ${currentyear} ${currentmon2}

# Task 2: Download NMME Data
# NOTE: This task has been temporarily removed as the downloading of all raw data will be handled prior to this workflow

# Task 3: Rescale and reorganize NMME Data
#sh FORECAST_TASK_03.sh $currentmon1 $currentyear

#
##-----PART B: BIAS CORRECTION-----##
#
# This sections involves creating monthly bias corrected data for CFSv2 (non-precip)
# and NMME (precip) datasets
#

# Task 4: Monthly "BC" step applied to CFSv2
#sh FORECAST_TASK_04.sh ${SYR} ${EYR} ${CLIM_SYR} ${CLIM_EYR} ${currentmon2} ${currentmon1} ${lat1} ${lat2} ${lon1} ${lon2} ${leadmon} ${bec}

# Task 5: Monthly "BC" step applied to NMME
#for mod in "${nmmemodel[@]}"; do
#     sh FORECAST_TASK_05.sh $SYR $EYR $CLIM_SYR $CLIM_EYR $currentmon2 $currentmon1 ${lat1} ${lat2} ${lon1} ${lon2} ${mod} ${leadmon}
#done

##-----PART C: TEMPORAL DISAGGREGATION-----##
#
# This section involves temporally disaggregating the monthly bias corrected forecasts
# generated in PART B to sub-daily (6-hourly) resolution
#

# Task 6: CFSv2 Temporal Disaggregation
#sh FORECAST_TASK_06.sh ${SYR} ${EYR} ${currentmon} ${currentmon1} ${lat1} ${lat2} ${lon1} ${lon2} ${fcstdatatype} ${leadmon} ${bec}
#
# Task 7: Generate symbolic links to sub-daily CFSv2 BC forecasts for NMME temporal disaggregation due to an uneven number of ensembles between the datasets
#sh FORECAST_TASK_07.sh ${currentyear} ${currentmon2}
#
# Task 8: NMME Temporal Disaggregation
#for mod in "${nmmemodel[@]}"; do
#    sh FORECAST_TASK_08.sh $SYR $EYR $currentmon $currentmon1 ${lat1} ${lat2} ${lon1} ${lon2} ${fcstdatatype1} ${mod} ${leadmon}
#done

##-----PART D: LIS PREPARATION-----##
#
# This section involves preparing the sub-daily bias corrected forecasts for LIS
#

# Task 9: Combine the CFSv2 forcing fields into final format for LIS to read
#sh FORECAST_TASK_09.sh ${SYR} ${EYR} ${currentmon} ${currentmon1} ${lat1} ${lat2} ${lon1} ${lon2} ${fcstdatatype} ${leadmon} ${bec

# Task 10: Combine the NMME forcing fields into final format for LIS to read and symbolically link to the reusable CFSv2 met forcings
#for ((YEAR=$SYR; YEAR<=$EYR; YEAR++)); do
# for mod in "${nmmemodel[@]}"; do
#   sh FORECAST_TASK_10.sh ${currentmon1} ${YEAR} ${currentmon} ${fcstdatatype} ${mod}
# done
#done

# Task 11: Copy 9th forecast lead file as 10th forecast lead for LIS runs
#
#for ((YEAR=$SYR; YEAR<=$EYR; YEAR++)); do
#  sh FORECAST_TASK_11.sh ${currentmon2} ${YEAR}
#done

# Task 12: Temporary task to introduce an all-zero variable V10M due to the way wind is handled in the USAF forcing
#sh FORECAST_TASK_12.sh ${currentmon2} ${currentmon1} ${bec}


##--------------Check Final Preprocessed Files and Folders-----------------##
# A final processing check is available to the user to check that the number of files and file size generated within this workflow meet expectations
#for ((YEAR=$SYR; YEAR<=$EYR; YEAR++)); do
#  for mod in "${nmmemodel[@]}"; do
#     sh Check_Preprocess_Finalfiles.F.sh ${currentmon2} ${YEAR} ${mod}
#  done
#done

##-------------- Completed Preprocessing -----------------##
