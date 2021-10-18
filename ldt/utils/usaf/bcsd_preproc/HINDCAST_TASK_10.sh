#!/bin/sh
#
# Part3-b: Combine all non-precip 6-hourly files into one file.
#  and copy BCSD precip files in to the same directory
#
# User-specifications:
   FCST_SYR=${1}
   FCST_EYR=${2}
   m=${3}
   mon=${m,}
   iMon=${mon}"01"
   iMonNo=${4}
   lead_months=${10}
   ens_num=${11}

#  Domain:
#  FAME DOMAIN ...
   lat1=${5}
   lat2=${6}
   lon1=${7}
   lon2=${8}
#
# MODEL

   fcstdatatype=${9}
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

#  Final 6-hour forcing directory:
   FORCEDIR=${PROJDIR}'/data/CFSv2_25km/'

# Source and load required modules:
#source /usr/share/modules/init/sh
#module load other/comp/gcc-5.3-sp3
#module load other/cdo-1.7.1
#
#------------------------------------------------------------------------------
#  Combine all non-precip 6-hourly files into one file.
#  and copy BCSD precip files in to the same directory
#

cd $SRCDIR
ulimit -s unlimited

  for ((YEAR=$FCST_SYR; YEAR<=$FCST_EYR; YEAR++)); do
     echo "Combining subdaily BCSD forecast files"
     sbatch $SRCDIR/run_Combining.scr $SRCDIR $YEAR $iMonNo $ens_num $lead_months $FORCEDIR ${fcstdatatype} ${LOGDIR};
  done

#------------------------------------------------------------------------------

