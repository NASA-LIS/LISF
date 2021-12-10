#!/bin/sh
#
# Combine all non-precip 6-hourly files into one file.
# and copy BCSD precip files in to the same directory
#

currentyear=${1}
mon=${2}

for iMon in ${mon}01; do
  echo $iMon

  cd $SRCDIR
  ulimit -s unlimited

  for year in ${currentyear}; do

    INDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/forecasts/CFSv2_25km/raw/6-Hourly/'
    #please change path for OUTDIR
    OUTDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/forecasts/NMME/linked_cfsv2_precip_files/'
    echo $year

    cd $OUTDIR
    mkdir -p $iMon/$year

    cd ${OUTDIR}'/'${iMon}'/'${year}'/'

    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens1' 'ens1'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens2' 'ens2'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens3' 'ens3'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens4' 'ens4'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens5' 'ens5'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens6' 'ens6'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens7' 'ens7'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens8' 'ens8'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens9' 'ens9'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens10' 'ens10'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens11' 'ens11'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens12' 'ens12'

    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens1' 'ens13'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens2' 'ens14'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens3' 'ens15'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens4' 'ens16'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens5' 'ens17'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens6' 'ens18'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens7' 'ens19'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens8' 'ens20'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens9' 'ens21'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens10' 'ens22'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens11' 'ens23'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens12' 'ens24'

    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens1' 'ens25'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens2' 'ens26'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens3' 'ens27'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens4' 'ens28'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens5' 'ens29'
    ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens6' 'ens30'

    echo "done sym link"

  done

done

#------------------------------------------------------------------------------

