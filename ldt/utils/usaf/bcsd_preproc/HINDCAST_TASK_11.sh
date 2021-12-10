#!/bin/sh
#
# User-specifications:
  currentmon=${1}
  currentyear=${2}
  currentmon2=${3}
  currentmon1=${currentmon2,}

  echo " == Current Date: "$currentmon", "$currentyear", "$currentmon1", "$currentmon2

  fcst0=`date -d "${currentmon}/01/${currentyear} + 0 month" +%Y%m`
  fcst1=`date -d "${currentmon}/01/${currentyear} + 1 month" +%Y%m`
  fcst2=`date -d "${currentmon}/01/${currentyear} + 2 month" +%Y%m`
  fcst3=`date -d "${currentmon}/01/${currentyear} + 3 month" +%Y%m`
  fcst4=`date -d "${currentmon}/01/${currentyear} + 4 month" +%Y%m`
  fcst5=`date -d "${currentmon}/01/${currentyear} + 5 month" +%Y%m`
  fcst6=`date -d "${currentmon}/01/${currentyear} + 6 month" +%Y%m`
  fcst7=`date -d "${currentmon}/01/${currentyear} + 7 month" +%Y%m`
  fcst8=`date -d "${currentmon}/01/${currentyear} + 8 month" +%Y%m`
  fcst9=`date -d "${currentmon}/01/${currentyear} + 9 month" +%Y%m`
  echo " == FCSTDATE: "$fcst0", "$fcst1", "$fcst2", "$fcst3", "$fcst4", "$fcst5", "$fcst6, "$fcst7", "$fcst8", "$fcst9"

#  Forecast years to be processed:
   FCST_SYR=${currentyear}
   FCST_EYR=${currentyear}
   echo ${FCST_SYR}
   echo ${FCST_EYR}

   iMon=${currentmon1}"01"
   echo "iMon ="$iMon
   iMonNo=$currentmon
   echo "iMonNo ="$iMonNo

# MODEL

   fcstdatatype=${4}
   nmmemod=${5}
   basemodname="CFSv2"
   basemodname2="GEOS5"
   echo "fcstdatatype:"${fcstdatatype}
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
   FORCEDIR_NMME=${PROJDIR}'/data/NMME'
   FORCEDIR_CFSv2=${PROJDIR}'/data/CFSv2_25km'

#
# Source and load required modules:
source /usr/share/modules/init/sh
#
#------------------------------------------------------------------------------
#  Combine all non-precip 6-hourly files into one file.
#  and copy BCSD precip files in to the same directory
#

cd $SRCDIR
ulimit -s unlimited

for model in ${nmmemod}; do
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

   INDIR_NMME=${FORCEDIR_NMME}'/bcsd/6-Hourly/'${iMon}'/'${model}
   INDIR_CFSv2=${FORCEDIR_CFSv2}'/final/6-Hourly'

   #please change path for OUTDIR
   OUTDIR=${FORCEDIR_NMME}'/final/6-Hourly/'${model}

  for ((YEAR=$FCST_SYR; YEAR<=$FCST_EYR; YEAR++)); do
     echo "Copying subdaily BCSD Precip forecast files"
     sh $SRCDIR/Copy_Precipitation_sub-daily.scr $YEAR $iMon $ens_numf $INDIR_NMME $OUTDIR;

     echo "Copying combined subdaily BCSD forecast files"

  done

 if [ $model == CCM4 ] || [ $model == GNEMO ] || [ $model == CCSM4 ]; then
  k=1
  k1=1
  ens_numf=10
  #for i in 1 2 3 4 5 6 7 8 9 10; do
  for i in {1..10}; do
  echo $i;
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst0}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst0}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst1}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst1}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst2}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst2}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst3}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst3}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst4}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst4}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst5}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst5}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst6}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst6}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst7}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst7}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst8}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst8}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst8}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst9}'.nc4'
  k=`expr $k + $k1`
  done
 fi

 if [ $model == CFSv2 ]; then
  k=1
  k1=1
  ens_numf=12
  #for i in 1 2 3 4 5 6 7 8 9 10 11 12; do
  for i in {1..12}; do
  echo $i;
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst0}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst0}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst1}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst1}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst2}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst2}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst3}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst3}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst4}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst4}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst5}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst5}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst6}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst6}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst7}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst7}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst8}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst8}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst8}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst9}'.nc4'
  k=`expr $k + $k1`
  done
 fi


 if [ $model == GEOSv2 ]; then
  k=1
  k1=1
  ens_numf=4
  #for i in 1 2 3 4; do
  for i in {1..4}; do
  echo $i;
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst0}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst0}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst1}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst1}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst2}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst2}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst3}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst3}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst4}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst4}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst5}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst5}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst6}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst6}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst7}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst7}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst8}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst8}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst8}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst9}'.nc4'
  k=`expr $k + $k1`
  done
 fi


 if [ $model == GFDL ]; then
  k=1
  k1=1
  ens_numf=15
  for i in 1 2 3 4 5 6 7 8 9 10 11 12 2 3 4; do
  echo $i;
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst0}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst0}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst1}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst1}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst2}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst2}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst3}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst3}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst4}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst4}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst5}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst5}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst6}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst6}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst7}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst7}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst8}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst8}'.nc4'
  ln -sfn ${INDIR_CFSv2}'/'${currentyear}'/'${iMon}'/ens'${i}'/'${basemodname}'.'${fcst8}'.nc4' ${OUTDIR}'/'${currentyear}'/'${iMon}'/ens'${k}'/'${basemodname2}'.'${fcst9}'.nc4'
  k=`expr $k + $k1`
  done
 fi

done

#------------------------------------------------------------------------------

