#!/bin/sh
#
# This task processes the CFSv2 forecast data and outputs in 6-hourly and
# monthly time resolutions.
#
#------------------------------------------------------------------------------
#
# HELPER FUNCTIONS
#
# Generates forecast initialization dates based on the initialization month
function cal_ic_dates() {

    icmon=$1
    if [ $icmon == jan01 ]; then
        ic[0]='1217'
        ic[1]='1222'
        ic[2]='1227'
    elif [ $icmon == feb01 ]; then
        ic[0]='0121'
        ic[1]='0126'
        ic[2]='0131'
    elif [ $icmon == mar01 ]; then
        ic[0]='0215'
        ic[1]='0220'
        ic[2]='0225'
    elif [ $icmon == apr01 ]; then
        ic[0]='0317'
        ic[1]='0322'
        ic[2]='0327'
    elif [ $icmon == may01 ]; then
        ic[0]='0416'
        ic[1]='0421'
        ic[2]='0426'
    elif [ $icmon == jun01 ]; then
        ic[0]='0521'
        ic[1]='0526'
        ic[2]='0531'
    elif [ $icmon == jul01 ]; then
        ic[0]='0620'
        ic[1]='0625'
        ic[2]='0630'
    elif [ $icmon == aug01 ]; then
        ic[0]='0720'
        ic[1]='0725'
        ic[2]='0730'
    elif [ $icmon == sep01 ]; then
        ic[0]='0819'
        ic[1]='0824'
        ic[2]='0829'
    elif [ $icmon == oct01 ]; then
        ic[0]='0918'
        ic[1]='0923'
        ic[2]='0928'
    elif [ $icmon == nov01 ]; then
        ic[0]='1018'
        ic[1]='1023'
        ic[2]='1028'
    elif [ $icmon == dec01 ]; then
        ic[0]='1117'
        ic[1]='1122'
        ic[2]='1127'
    fi
}
#
#------------------------------------------------------------------------------
#
# INPUT ARGUMENTS AND PATH DIRECTORIES
#
# Input Arguments
#
    FCST_SYR=${1}
    FCST_EYR=${2}
    month_abbr=${3}
#
# Other Arguments
#
    iMon=${month_abbr}"01"
#
# Path of the main project directory
#
    PROJDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'
#
# Path of the directory where all the download codes are kept:
#
    SRCDIR=${PROJDIR}'/scripts/code_library'
#
# Paths for the daily forecast data (input and output paths):
#
    FORCEDIR='/discover/nobackup/projects/lis/MET_FORCING/CFSv2'
    OUTDIR=${PROJDIR}'/data/CFSv2_25km/raw'
    GRIDDESC=${SRCDIR}'/supplementary_files/CFSv2_25km_AFRICOM_grid_description.txt'
#
#  Log file output directory
#
    LOGDIR=${PROJDIR}'/scripts/log_files'
    mkdir -p ${LOGDIR}
#
#------------------------------------------------------------------------------
#
# Calls function 'cal_ic_dates' for the forecast initialization dates
#
    cal_ic_dates $iMon
#
# Source the correct path
#
#    source /usr/share/modules/init/sh
#
#------------------------------------------------------------------------------
#
# Process 3-hourly CFSv2 forecasts and output in monthly and 6-hourly formats
#
#   Process and convert raw 3-hourly CFSv2 forecasts:
    echo " -- Processing CFSv2 3-hourly forecast variables -- "
    for ((YEAR=$FCST_SYR; YEAR<=$FCST_EYR; YEAR++)); do
        sbatch $SRCDIR/run_process_forecast_data.scr $YEAR $YEAR $iMon $SRCDIR $OUTDIR $FORCEDIR $GRIDDESC ${ic[@]}
    done
#
#------------------------------------------------------------------------------
#
    echo " -- Jobs submitted to process CFSv2 forecast files for: "${iMon}" -- "
#
#------------------------------------------------------------------------------
#
