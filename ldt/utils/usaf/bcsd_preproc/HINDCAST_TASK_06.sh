#!/bin/sh
#
#  PLEASE READ FILE 'README_PART2'
#
#  COMPUTE THE BIAS CORRECTION FOR GEOS5.
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
   currentmon=${5}
   lead_months=${12}

#  Domain extents (from LL to UR):
#  FAME DOMAIN ...
   lat1=${7}
   lat2=${8}
   lon1=${9}
   lon2=${10}

   nmmemod=${11}

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
  OBS_CLIM_INDIR=${FORCEDIR}'/USAF-LIS7.3rc8_25km/raw/Climatology'
  FCST_CLIM_INDIR=${FORCEDIR}'/NMME/raw/Climatology/'${iMon}
  FCST_INDIR=${FORCEDIR}'/NMME/raw/Monthly/'${iMon}
#
# BC output directory for FCST:
#
  OUTDIR=${FORCEDIR}'/NMME/bcsd/Monthly/'${iMon}
  mkdir -p ${OUTDIR}
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
# Source and load required modules:
source /usr/share/modules/init/sh
module load python/GEOSpyD/Ana2019.10_py3.7
ulimit -s unlimited

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
   OBS_VAR_LIST=(PRECCON Rainf_f_tavg LWGAB SWGDN PS QV2M T2M U10M V10M)
   FCST_VAR_LIST=(CNPRCP PRECTOT LWS SLRSF PS Q2M T2M U10M V10M)
   UNIT=('kg/m^2/s' 'kg/m^2/s' 'W/m^2' 'W/m^2' 'Pa' 'kg/kg' 'K' 'm/s' 'm/s')

   for VAR_NUM in 1; do
      if [ $VAR_NUM == 0 ] || [ $VAR_NUM == 1 ] || [ $VAR_NUM == 3 ]; then
        VAR_TYPE='PRCP'
      else
        VAR_TYPE='TEMP'
      fi
      echo ${VAR_NUM}" "${FCST_VAR_LIST[$VAR_NUM]}
      for model in ${nmmemod}; do
      #for model in CFSv2 CCM4 GNEMO CCSM4 GFDL GEOSv2; do
          echo ${model}
	  if [ $model == CFSv2 ]; then
            ens_numf=12
	    ens_numc=12
	    enss=1
	    ensf=12
	  elif [ $model == GEOSv2 ]; then
            ens_numf=4
	    ens_numc=4
	    enss=13
	    ensf=16
	  elif [ $model == CCM4 ]; then
            ens_numf=10
	    ens_numc=10
	    enss=17
	    ensf=26
	  elif [ $model == GNEMO ]; then
            ens_numf=10
	    ens_numc=10
	    enss=27
	    ensf=36
	  elif [ $model == CCSM4 ]; then
            ens_numf=10
	    ens_numc=10
	    enss=37
	    ensf=46
	  elif [ $model == GFDL ]; then
	    ens_numf=15
            ens_numc=15
	    enss=47
	    ensf=61
	  fi

      echo $ens_numf

      script='Bias_correction_NMME_modulefast.py'


    echo ${MASK_FILE}
    echo ${FCST_CLIM_INDIR}
    echo ${FCST_INDIR}
    echo ${OUTDIR}
    echo ${script}


    for ((YEAR=$FCST_SYR; YEAR<=$FCST_EYR; YEAR++)); do
         echo 'year='$YEAR
         sbatch $SRCDIR/run_NMME_BCSD_calctest.scr ${SRCDIR} ${OBS_VAR_LIST[$VAR_NUM]} ${FCST_VAR_LIST[$VAR_NUM]} ${iMonNo} ${VAR_TYPE} ${UNIT[$VAR_NUM]} ${lat1} ${lat2} ${lon1} ${lon2} ${ens_numc} ${ens_numf} ${model} ${lead_months} ${YEAR} ${YEAR} ${CLIM_SYR} ${CLIM_EYR} ${MASK_FILE} ${FCST_CLIM_INDIR} ${OBS_CLIM_INDIR} ${FCST_INDIR} ${OUTDIR} ${LOGDIR} ${script} ${enss} ${ensf}

    done

  done

done


# -------

   echo " -- Completed processing BCSD forcing files for: "${iMon}" -- "

#------------------------------------------------------------------------------
