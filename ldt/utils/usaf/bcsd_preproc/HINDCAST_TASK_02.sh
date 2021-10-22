#!/bin/sh
#
# This task generates climatologies for the USAF-LIS7.3rc8 observation data
#
#------------------------------------------------------------------------------
#
# INPUT ARGUMENTS AND PATH DIRECTORIES
#
# Forecast years to be processed:
#
    m=${1}
    mon=${m,}
    iMon=${mon}"01"
    iMonNo=${2}

#  Years specified for generating the climatologies:
   CLIM_SYR=${3}
   CLIM_EYR=${4}

#  Forecast data specifications:
   ens_num=${5}
   lead_months=${6}

#  Domain extents
   lat1=${7}
   lat2=${8}
   lon1=${9}
   lon2=${10}

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
    FORCEDIR=${PROJDIR}'/data/USAF-LIS7.3rc8_25km/raw'
#
#  Climatology output directory for USAF-LIS7.3rc8:
#
    OUTDIR=${FORCEDIR}'/Climatology'
#
# Paths to the USAF-LIS7.3rc8 data to calculate the climatologies:
#
    OBS_INDIR=${FORCEDIR}'/Monthly/SURFACEMODEL'
    OBS_MASK=${SRCDIR}'/supplementary_files/Mask_nafpa.nc'
#
#  Log file output directory
#
    LOGDIR=${PROJDIR}'/scripts/log_files'
    mkdir -p ${LOGDIR}
#
#------------------------------------------------------------------------------
#
#   Calculate observed climatology for USAF-LIS7.3rc8 fields:
#   Note that below script uses Mask_merra2.nc as a mask to get climatology
#   over land grid cells only.
#
   echo " -- Processing observed climatology of USAF-LIS7.3rc8 variables -- "

   for VAR in LWdown_f_tavg  Rainf_f_tavg  Psurf_f_tavg  Qair_f_tavg  SWdown_f_tavg  Tair_f_tavg  Wind_f_tavg; do
     echo $VAR
     python $SRCDIR/Calc_and_Write_observational_climatology.py $VAR $lat1 $lat2 $lon1 $lon2 $OBS_INDIR $CLIM_SYR $CLIM_EYR $OBS_MASK $OUTDIR

#> $LOGDIR/Calc_obs_${VAR}_clim.log
   done;
#
#------------------------------------------------------------------------------
#
    echo " -- Completed Processing GEOS Climatology Forcing Files -- "
#
#------------------------------------------------------------------------------
#

#for mon in Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec; do
# echo $mon
# cp /discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/FAME_Mar_data/GEOS5v2/CLIM_DATA/OBS/* /discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/FAME_${mon}_data/GEOS5v2/CLIM_DATA/OBS/.
# cp /discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/FAME_Mar_data/GEOS5v2/CLIM_DATA/OBS/PRECTOT_obs_clim.nc /gpfsm/dnb02/projects/p63/FORECASTS/GEOS5/BCSD_Test/EXPERIMENTS/NMME/data/FAME_${mon}_data/nmme/CLIM_DATA/OBS/.
#done;


