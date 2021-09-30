#!/bin/sh
#
# This task generates climatologies for the CFSv2 forecast data
#
# -----------------------------------------------------------------------
#
# INPUT ARGUMENTS AND PATH DIRECTORIES
#
#  Forecast years to be processed:
#
   FCST_SYR=${1}
   FCST_EYR=${2}
   m=${5}
   mon=${m,}
   iMon=${mon}"01"
   iMonNo=${6}

#  Years specified for generating the climatologies:
   CLIM_SYR=${3}
   CLIM_EYR=${4}

#  Forecast data specifications:
   ens_num=${7}
   lead_months=${8}

#  Domain extents
   lat1=${9}
   lat2=${10}
   lon1=${11}
   lon2=${12}

#
# Path of the main project directory
#
    PROJDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM'
#
# Path of the directory where all the climatology and BCSD codes are kept:
#
    SRCDIR=${PROJDIR}'/scripts/code_library'
#
# Path for where observational files are located:
#
    FORCEDIR=${PROJDIR}'/data/CFSv2_25km/raw'
#
#  Climatology output directory for CFSv2:
#
    OUTDIR=${FORCEDIR}'/Climatology/'${iMon}
#
# Paths to the CFSv2 data to calculate the climatologies
#
    FCST_INDIR=${FORCEDIR}'/Monthly'
    FCST_MASK=${SRCDIR}'/supplementary_files/Mask_nafpa.nc'
#
#  Log file output directory
#
   LOGDIR=${PROJDIR}'/scripts/log_files'
   mkdir -p ${LOGDIR}
#
#------------------------------------------------------------------------------
#
#   Calculate the climatology for CFSv2 forecast fields:
#   Note that below script uses chirps-v2.0.198201_p25.nc and
#   Mask_merra2.nc as masks to get climatology over land grid cells only.
#
   echo " -- Processing forecast climatology of CFSv2 variables -- "

   #for VAR in PRECTOT  LWS  PS  Q2M  SLRSF  T2M  U10M  V10M; do
   for VAR in PRECTOT  LWS  PS  Q2M  SLRSF  T2M  WIND10M; do
     echo $VAR
     python $SRCDIR/Calc_and_Write_forecast_climatology.py $VAR $lat1 $lat2 $lon1 $lon2 $iMonNo $FCST_SYR $FCST_EYR $lead_months $ens_num $CLIM_SYR $CLIM_EYR $FCST_INDIR $FCST_MASK $OUTDIR

#> $LOGDIR/Calc_fcst_clim_$VAR3.log
   done;
#
#------------------------------------------------------------------------------

  echo " -- Completed Processing Climatology Forcing Files -- "

#------------------------------------------------------------------------------

