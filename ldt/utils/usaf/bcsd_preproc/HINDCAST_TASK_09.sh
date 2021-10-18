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
   lead_months=${11}

# Run domain:
#  Domain extents (from LL to UR):
#  FAME DOMAIN ...
   lat1=${5}
   lat2=${6}
   lon1=${7}
   lon2=${8}
   nmmemod=${10}
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
# Path for where raw and bias corrected nmme forecast files are located:
#
    FORCEDIR=${PROJDIR}'/data/NMME/'
    SUBDAILY_RAW_FCST_DIR=${FORCEDIR}'/linked_cfsv2_precip_files/'${iMon}
    MONTHLY_RAW_FCST_DIR=${FORCEDIR}'/raw/Monthly/'${iMon}
    MONTHLY_BC_FCST_DIR=${FORCEDIR}'/bcsd/Monthly/'${iMon}

# (2) MERRA2 and CHIRPS masks

   MASK_FILE_PRECIP=${SRCDIR}'/supplementary_files/Mask_nafpa.nc'
   MASK_FILE_NONPRECIP=${SRCDIR}'/supplementary_files/Mask_nafpa.nc'

#
#------------------------------------------------------------------------------
#
#  Temporally downscale the monthly bias-corrected forecasts to daily
#   and then to sub-daily output files.

   echo " -- Downscale the monthly BCSD forecasts to daily and sub-daily output --"

   #OBS_VAR='Rainf_f_tavg'
   OBS_VAR='PRECTOT'
   FCST_VAR='PRECTOT'
   UNITS='kg/m^2/s'

   #OBS_VAR_LIST=(Rainf_f_tavg LWdown_f_tavg SWdown_f_tavg Psurf_f_tavg Qair_f_tavg Tair_f_tavg Wind_f_tavg)
   #FCST_VAR_LIST=(PRECTOT LWS SLRSF PS Q2M T2M WIND10M)
   #UNITS=('kg/m^2/s' 'W/m^2' 'W/m^2' 'Pa' 'kg/kg' 'K' 'm/s')

   VAR_TYPE='PRCP'
   VAR_NUM=0

     if [ $VAR_NUM == 0 ] || [ $VAR_NUM == 1 ] || [ $VAR_NUM == 3 ]; then
 	VAR_TYPE='PRCP'
     else
	VAR_TYPE='TEMP'
     fi

     echo ${VAR_NUM}" "${FCST_VAR_LIST[$VAR_NUM]}
     for model in ${nmmemod}; do
       echo ${model}

       if [ $model == CFSv2 ]; then
         ens_num=12
	     #enss=1
	     #ensf=12
	   elif [ $model == GEOSv2 ]; then
         ens_num=4
	     #enss=1
	     #ensf=4
	   elif [ $model == CCM4 ]; then
         ens_num=10
	     #enss=1
	     #ensf=10
	   elif [ $model == GNEMO ]; then
         ens_num=10
	     #enss=1
	     #ensf=10
	   elif [ $model == CCSM4 ]; then
         ens_num=10
	     #enss=1
	     #ensf=10
	   elif [ $model == GFDL ]; then
	     ens_num=15
	     #enss=1
	     #ensf=15
	   fi


	   echo $ens_num
	   #
	   #please change path for OUTDIR
	   OUTDIR=${FORCEDIR}'/bcsd/6-Hourly/'${iMon}'/'${model}
	   mkdir -p ${OUTDIR}

	   for ((YEAR=$FCST_SYR; YEAR<=$FCST_EYR; YEAR++)); do

        sbatch $SRCDIR/run_Temporal_disagg.scr ${SRCDIR} ${OBS_VAR} ${FCST_VAR} $iMonNo $VAR_TYPE ${UNITS} $lat1 $lat2 $lon1 $lon2 $model $ens_num $lead_months $YEAR $YEAR $MASK_FILE_PRECIP $MASK_FILE_NONPRECIP $MONTHLY_BC_FCST_DIR $MONTHLY_RAW_FCST_DIR $SUBDAILY_RAW_FCST_DIR $OUTDIR $LOGDIR

	   done;
 done;

   echo " -- Completed Submitting Sbatch Scripts for Temporal Disaggregation -- "

#------------------------------------------------------------------------------
