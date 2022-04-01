#!/bin/sh
#
#  PLEASE READ FILE 'README_PART2'
#
#  COMPUTE THE BIAS CORRECTION FOR CFSv2.
#
# -----------------------------------------------------------------------

#  Forecast years to be processed:
   FCST_SYR=${1}
   FCST_EYR=${2}

#  Years specified for generating the climatologies:
   CLIM_SYR=${3}
   CLIM_EYR=${4}

#  Forecast data specifications:
   iMon=${5}"01"
   iMonNo=${6}
   lead_months=${11}
   ens_numc=${12}

   echo ${ens_numc}

#  Domain extents (from LL to UR)
   lat1=${7}
   lat2=${8}
   lon1=${9}
   lon2=${10}

#
# Path of the main project directory
#
    PROJDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'
#
# Path of the directory where all the BC codes are kept:
#
    SRCDIR=${PROJDIR}'/scripts/code_library'
#
# Path for where observational & forecast files are located:
#
    FORCEDIR=${PROJDIR}'/data'
    OBS_INDIR=${FORCEDIR}'/USAF-LIS7.3rc8_25km'
    FCST_INDIR=${FORCEDIR}'/CFSv2_25km'
#
# BC output directory for FCST:
#
    OUTDIR=${FCST_INDIR}'/bcsd/Monthly/'${iMon}
#
# Mask file
#
   MASK_FILE=${SRCDIR}'/supplementary_files/Mask_nafpa.nc'
#
#  Log file output directory
#
   LOGDIR=${SRCDIR}'/Log_Files'
   mkdir -p ${LOGDIR}
#
#
#------------------------------------------------------------------------------
#
#   Perform bias corrections, using observed and forecast sorted
#    climatologies, and target forecasts
#
#   Note that below script uses CHIRPS_0.25_MASK.nc and CHIRPS_MASK.nc as
#   masks for PRECTOT and the other variables, respectively.
#
   echo " -- Processing forecast bias correction of GEOS5.0 variables -- "

#  Calculate bias correction for different variables separately:
   OBS_VAR_LIST=(Rainf_f_tavg  LWdown_f_tavg  SWdown_f_tavg  Psurf_f_tavg  Qair_f_tavg  Tair_f_tavg  Wind_f_tavg)
   FCST_VAR_LIST=(PRECTOT LWS SLRSF PS Q2M T2M WIND10M)
   UNIT=('kg/m^2/s' 'W/m^2' 'W/m^2' 'Pa' 'kg/kg' 'K' 'm/s')

   for VAR_NUM in 0 1 2 3 4 5 6; do
      if [ $VAR_NUM == 0 ] || [ $VAR_NUM == 2 ]; then
        VAR_TYPE='PRCP'
      else
        VAR_TYPE='TEMP'
      fi
       echo ${VAR_NUM}" "${FCST_VAR_LIST[$VAR_NUM]}

      for ((YEAR=$FCST_SYR; YEAR<=$FCST_EYR; YEAR++)); do
         echo 'year='$YEAR

         sbatch $SRCDIR/run_BCSD_calctest.scr ${SRCDIR} ${OBS_VAR_LIST[$VAR_NUM]} ${FCST_VAR_LIST[$VAR_NUM]} ${iMonNo} ${VAR_TYPE} ${UNIT[$VAR_NUM]} ${lat1} ${lat2} ${lon1} ${lon2} ${ens_numc} ${lead_months} ${YEAR} ${YEAR} ${CLIM_SYR} ${CLIM_EYR} ${MASK_FILE} ${OBS_INDIR} ${FCST_INDIR} ${OUTDIR} ${LOGDIR}
         #sh $SRCDIR/run_BCSD_calctest.scr ${SRCDIR} ${OBS_VAR_LIST[$VAR_NUM]} ${FCST_VAR_LIST[$VAR_NUM]} ${iMonNo} ${VAR_TYPE} ${UNIT[$VAR_NUM]} ${lat1} ${lat2} ${lon1} ${lon2} ${ens_numc} ${lead_months} ${YEAR} ${YEAR} ${CLIM_SYR} ${CLIM_EYR} ${MASK_FILE} ${OBS_INDIR} ${FCST_INDIR} ${OUTDIR} ${LOGDIR}
      done

   done

# -------

   echo " -- Completed processing BCSD forcing files for: "${iMon}" -- "

#------------------------------------------------------------------------------

