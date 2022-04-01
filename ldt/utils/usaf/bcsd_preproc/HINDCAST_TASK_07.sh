#!/bin/sh

ulimit -s unlimited

# Part3-a: Generate bias-corrected 6-hourly forecasts using raw monthly forecasts, bias-corrected
#          monthly forecasts and raw 6-hourly forecasts.
#

# 1. User-specified entries:
   FCST_SYR=${1}
   FCST_EYR=${2}
   fcstdatatype=${9}

# GEOS5 forecast information:
   m=${3}
   mon=${m,}
   iMon=${mon}"01"
   iMonNo=${4}
   lead_months=${10}
   ens_num=${11}

# Run domain:
#  Domain extents (from LL to UR):
#  FAME DOMAIN ...
   lat1=${5}
   lat2=${6}
   lon1=${7}
   lon2=${8}
#
# Path of the main project directory
#
    PROJDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'
#
# Path of the directory where all the BC codes are kept:
#
    SRCDIR=${PROJDIR}'/scripts/code_library'
#
#  Log file output directory
#
   LOGDIR=${SRCDIR}'/Log_Files'
   mkdir -p ${LOGDIR}
#
# Path for where raw and bias corrected forecast files are located:
#
   FORCEDIR=${PROJDIR}'/data/CFSv2_25km'
   SUBDAILY_RAW_FCST_DIR=${FORCEDIR}'/raw/6-Hourly/'${iMon}
   MONTHLY_RAW_FCST_DIR=${FORCEDIR}'/raw/Monthly/'${iMon}
   MONTHLY_BC_FCST_DIR=${FORCEDIR}'/bcsd/Monthly/'${iMon}

# (2) MERRA2 and CHIRPS masks

   MASK_FILE_PRECIP=${SRCDIR}'/supplementary_files/Mask_nafpa.nc'
   MASK_FILE_NONPRECIP=${SRCDIR}'/supplementary_files/Mask_nafpa.nc'

   #please change path for OUTDIR
   OUTDIR=${FORCEDIR}'/bcsd/6-Hourly/'${iMon}
   mkdir -p ${OUTDIR}

#
#------------------------------------------------------------------------------
#
#  Temporally downscale the monthly bias-corrected forecasts to daily
#   and then to sub-daily output files.

   echo " -- Downscale the monthly BCSD forecasts to daily and sub-daily output --"

   #OBS_VAR_LIST=(LWdown_f_tavg SWdown_f_tavg Psurf_f_tavg Qair_f_tavg Tair_f_tavg Wind_f_tavg)
   OBS_VAR_LIST=(LWGAB SWGDN PS QV2M T2M U10M)
   FCST_VAR_LIST=(LWS SLRSF PS Q2M T2M WIND10M)
   UNITS=('kg/m^2/s' 'W/m^2' 'W/m^2' 'Pa' 'kg/kg' 'K' 'm/s')
   for ((YEAR=$FCST_SYR; YEAR<=$FCST_EYR; YEAR++)); do
     for VAR_NUM in 0 1 2 3 4 5; do

       if [ $VAR_NUM == 1 ]; then
 	     VAR_TYPE='PRCP'
       else
	     VAR_TYPE='TEMP'
       fi

       sbatch $SRCDIR/run_Temporal_disagg.scr ${SRCDIR} ${OBS_VAR_LIST[$VAR_NUM]} ${FCST_VAR_LIST[$VAR_NUM]} $iMonNo $VAR_TYPE ${UNITS[$VAR_NUM]} $lat1 $lat2 $lon1 $lon2 $fcstdatatype $ens_num $lead_months $YEAR $YEAR $MASK_FILE_PRECIP $MASK_FILE_NONPRECIP $MONTHLY_BC_FCST_DIR $MONTHLY_RAW_FCST_DIR $SUBDAILY_RAW_FCST_DIR $OUTDIR $LOGDIR
       #sh $SRCDIR/run_Temporal_disagg.scr ${SRCDIR} ${OBS_VAR_LIST[$VAR_NUM]} ${FCST_VAR_LIST[$VAR_NUM]} $iMonNo $VAR_TYPE ${UNITS[$VAR_NUM]} $lat1 $lat2 $lon1 $lon2 $fcstdatatype $ens_num $lead_months $YEAR $YEAR $MASK_FILE_PRECIP $MASK_FILE_NONPRECIP $MONTHLY_BC_FCST_DIR $MONTHLY_RAW_FCST_DIR $SUBDAILY_RAW_FCST_DIR $OUTDIR $LOGDIR

     done;
   done;

   echo " -- Completed Submitting Sbatch Scripts for Temporal Disaggregation -- "

#------------------------------------------------------------------------------
