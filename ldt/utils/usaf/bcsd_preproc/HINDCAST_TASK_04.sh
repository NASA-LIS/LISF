#!/bin/sh
#
#  This task generates climatologies for the NMME precipitation data
#
# -----------------------------------------------------------------------
#
# INPUT ARGUMENTS AND PATH DIRECTORIES
#
#  Forecast years to be processed:
   FCST_SYR=${1}
   FCST_EYR=${2}
   m=${5}
   mon=${m,}
   iMon=${mon}"01"
   iMonNo=${6}

#  Years specified for generating the climatologies:
   CLIM_SYR=${3}
   CLIM_EYR=${4}

#  Domain extents (from LL to UR)
   lat1=${7}
   lat2=${8}
   lon1=${9}
   lon2=${10}

#  Forecast data specifications:
   nmmemodel=${11}
   lead_months=${12}

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
    FORCEDIR=${PROJDIR}'/data/NMME'
#
# Paths to the CHIRPS, MERRA and GEOS5 data to calculate their climatologies:
#
   NMME_MASK=${SRCDIR}'/supplementary_files/Mask_nafpa.nc'
#
# Paths to the NMME data to calculate the climatologies:
#
   NMME_INDIR=${FORCEDIR}'/raw/Monthly/'${iMon}
#
#  Log file output directory
#
   LOGDIR=${SRCDIR}'/Log_Files'
   mkdir -p ${LOGDIR}
#
#------------------------------------------------------------------------------
#
#   Calculate the climatology for NMME forecast precipitation:
#   Note that below script uses chirps-v2.0.198201_p25.nc
#   as masks to get climatology over land grid cells only.
#
  echo " -- Processing Precipitation forecast climatology of NMME variables -- "

    for VAR in PRECTOT; do
      echo $VAR
      for model in ${nmmemodel}; do
          echo ${model}
	  if [ $model == CFSv2 ]; then
            ens_numf=12
	    enss=1
	    ensf=12
	  elif [ $model == GEOSv2 ]; then
            ens_numf=4
	    enss=13
	    ensf=16
	  elif [ $model == CCM4 ]; then
            ens_numf=10
	    enss=17
	    ensf=26
	  elif [ $model == GNEMO ]; then
            ens_numf=10
	    enss=27
	    ensf=36
	  elif [ $model == CCSM4 ]; then
            ens_numf=10
	    enss=37
	    ensf=46
	  elif [ $model == GFDL ]; then
            ens_numf=15
	    enss=47
	    ensf=61
	  fi

      echo $ens_numf
#
#  Climatology output directory for NMME:
#
      OUTDIR=${FORCEDIR}'/raw/Climatology/'${iMon}'/'${model}
#

      python $SRCDIR/Calc_and_Write_NMME_forecast_climatology.py $VAR $lat1 $lat2 $lon1 $lon2 $iMonNo $FCST_SYR $FCST_EYR $lead_months $ens_numf $CLIM_SYR $CLIM_EYR $NMME_INDIR $NMME_MASK $model $OUTDIR $enss $ensf

#> $LOGDIR/Calc_fcst_clim_$VAR.log
    done
  done

  echo " -- Completed Processing Climatology Forcing Files -- "

#------------------------------------------------------------------------------
