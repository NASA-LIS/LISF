#!/bin/sh
# The following function identifiies the associated initial condition dates
#  for the different ensembles as used by NMME given an initialization month:

#month_number=$1
month_number=11

#currentmon2=$1
#currentmon=Jan
#currentmon1=01
declare -a month_abbreviations=('Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec')
currentmon=${month_abbreviations[$month_number - 1]}  # e.g. Jan
currentmon1=$(printf "%02d" $month_number)            # e.g. 01
currentmon2=${currentmon,}                            # e.g. jan

CLIM_SYR=2008
CLIM_EYR=2020

SYR=2008
EYR=2020
#SYR=$2
#EYR=$2

fcstdatatype="CFSv2"
#basemodelname="CFS"
fcstdatatype1="nmme"
declare -a nmmemodel=('CFSv2' 'GEOSv2' 'CCSM4' 'CCM4' 'GNEMO' 'GFDL')

#  Domain extents (from LL to UR): FAME Domain
lat1=-40
lat2=40
lon1=-20
lon2=60

leadmon=9
bec=12
#bef=4

source /usr/share/modules/init/sh
ulimit -s unlimited
module load python/GEOSpyD/Ana2019.10_py3.7

################~~~~~PREPROCESSING~~~~~################


##-----PART 1: DOWNSCALE, RESCALE, AND REORGANIZE-----##

# Task 1: Generate 6-hourly files - "SD" part; CFSv2 1Deg--> 25KM, and USAF 10KM --> 25KM; netcdf
#sh HINDCAST_TASK_01.sh ${SYR} ${EYR} ${currentmon2}

# Task 2: Generate monthly files for USAF (3-hrly) and CFSv2 (6-hrly); netcdf output
# *** This was done outside of these scripts as a preprocess for the hindcast data

# Task 3: Reorganize inputs for clim / BC steps
# *** Based on the changes to Task 1 and 2, there doesn't appear to be anything to do in this step


##-----PART 2: CLIMATOLOGY-----##

# Task 4: Generate monthly climatology files for USAF forcing (For Jan-2008 to Dec-2020)
#sh HINDCAST_TASK_02.sh $currentmon $currentmon1 $CLIM_SYR $CLIM_EYR ${bec} ${leadmon} ${lat1} ${lat2} ${lon1} ${lon2}

# Task 5: Generate monthly climatology files for CFSv2 (For Jan-2008 to Dec-2020)
#sh HINDCAST_TASK_03.sh ${SYR} ${EYR} ${CLIM_SYR} ${CLIM_EYR} ${currentmon} ${currentmon1} ${bec} ${leadmon} ${lat1} ${lat2} ${lon1} ${lon2}

# Task 6: Generate monthly NMME precip climatology files (For Jan-2008 to Dec-2020)
#for mod in "${nmmemodel[@]}"; do
#    sh HINDCAST_TASK_04.sh ${SYR} ${EYR} ${CLIM_SYR} ${CLIM_EYR} ${currentmon} ${currentmon1} ${lat1} ${lat2} ${lon1} ${lon2} ${mod} ${leadmon}
#done;

##-----PART 3: BIAS CORRECTION-----##

# Task 7: Monthly "BC" step applied for CFSv2
#sh HINDCAST_TASK_05.sh ${SYR} ${EYR} ${CLIM_SYR} ${CLIM_EYR} ${currentmon2} ${currentmon1} ${lat1} ${lat2} ${lon1} ${lon2} ${leadmon} ${bec}

# Task 8: Monthly "BC" step applied to NMME
#for mod in "${nmmemodel[@]}"; do
#    sh HINDCAST_TASK_06.sh $SYR $EYR $CLIM_SYR $CLIM_EYR $currentmon2 $currentmon1 ${lat1} ${lat2} ${lon1} ${lon2} ${mod} ${leadmon}
#done

##-----PART 4: TEMPORAL DISAGGREGATION-----##

# Task 9: Temporal dissagg. of CFS/NMME fields and members to 6-hourly netcdf output
#
# A: CFSv2 Temporal Disaggregation
#sh HINDCAST_TASK_07.sh ${SYR} ${EYR} ${currentmon} ${currentmon1} ${lat1} ${lat2} ${lon1} ${lon2} ${fcstdatatype} ${leadmon} ${bec}
#
# B: Symbolic link to daily data for NMME Temporal Disaggregation
#sh HINDCAST_TASK_08.sh
#
# C: NMME Temporal Disaggregation
#
#for mod in "${nmmemodel[@]}"; do
#    sh HINDCAST_TASK_09.sh $SYR $EYR $currentmon $currentmon1 ${lat1} ${lat2} ${lon1} ${lon2} ${fcstdatatype1} ${mod} ${leadmon}
#done


##-----PART 5: COMBINE FORCING FIELDS-----##

# Task 10: Combine the forcing fields into final format for LIS to read
#
#5.a CFSv2 Combine Variables
#
#sh HINDCAST_TASK_10.sh ${SYR} ${EYR} ${currentmon} ${currentmon1} ${lat1} ${lat2} ${lon1} ${lon2} ${fcstdatatype} ${leadmon} ${bec}

#5.b NMME Combine Variables (symbolically link to CFSv2 met forcings)
#
#for ((YEAR=$SYR; YEAR<=$EYR; YEAR++)); do
# for mod in "${nmmemodel[@]}"; do
#   sh HINDCAST_TASK_11.sh ${currentmon1} ${YEAR} ${currentmon} ${fcstdatatype} ${mod}
# done
#done

#5.c Copy 9th forecast lead file as 10th forecast lead for LIS runs
#
#for ((YEAR=$SYR; YEAR<=$EYR; YEAR++)); do
#  sh HINDCAST_TASK_12.sh ${currentmon2} ${YEAR}
#done


##--------------Check Final Preprocessed Files and Folders-----------------##

#for ((YEAR=$SYR; YEAR<=$EYR; YEAR++)); do
#  for mod in "${nmmemodel[@]}"; do
#     sh Check_Preprocess_Finalfiles.H.sh ${currentmon2} ${YEAR} ${mod}
#  done
#done
