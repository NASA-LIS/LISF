 #!/bin/sh
#
# This script creates symbolic link of Daily GEOS forecasts for NMME
#

for iMon in jan01 feb01 mar01 apr01 may01 jun01 jul01 aug01 sep01 oct01 nov01 dec01; do

 echo $iMon

 cd $SRCDIR
 ulimit -s unlimited

 for year in 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020; do

   INDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/CFSv2_25km/raw/6-Hourly/'
   #please change path for OUTDIR
   OUTDIR='/discover/nobackup/projects/usaf_lis/razamora/GHI_S2S/AFRICOM/data/NMME/linked_cfsv2_precip_files/'
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

   ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens2' 'ens13'
   ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens3' 'ens14'
   ln -sfn ${INDIR}'/'${iMon}'/'${year}'/ens4' 'ens15'

   echo "done sym link"

 done

done

#------------------------------------------------------------------------------

