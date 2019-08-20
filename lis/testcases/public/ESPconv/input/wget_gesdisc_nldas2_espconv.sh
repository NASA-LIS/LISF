#!/bin/bash
#
#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center Land Information System (LIS) v7.2
#
# Copyright (c) 2015 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
# Settings for wgetting GES DISC NLDAS-2 gridded data files 
#
# -- Script can run simply within your target directory, or
#     can include a target directory (uncomment {targetdir} lines below)
# -- This script is specially designed to download just days of the
#     the years of interest for ESP forecast run mode.
#
#  02.17.2012 -- Kristi Arsenault: v1.0 
#  07.04.2014 -- David Mocko:      v1.1: Public testcase for
#                                        LIS-7.0 Noah-3.3
#  03.17.2017 -- Kristi Arsenault: v1.2: Updated for ESP convention case
#
### Necessary directory paths ###
# Where actual full archive resides - Need a URS Earthdata Portal account to access:
#  serverhost=https://hydro1.gesdisc.eosdis.nasa.gov
#  sourcedira=data/NLDAS/NLDAS_FORA0125_H.002
#  sourcedirb=data/NLDAS/NLDAS_FORB0125_H.002

# Test case data:
  serverhost=https://portal.nccs.nasa.gov
  sourcedira=lisdata_pub/data/MET_FORCING/NLDAS2.FORCING
  sourcedirb=lisdata_pub/data/MET_FORCING/NLDAS2.FORCING

  rundir=$PWD
  mkdir NLDAS2.FORCING
  cd NLDAS2.FORCING
#  targetdir=./
#  cd $targetdir

# Date inputs:
  startyr=2010
  startmo=jun    #jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec
  startdy=01
  finalyr=2013
  finalmo=jun
  finaldy=05

################################################################

  fn=`date --date=$startdy' '$startmo' '$startyr +%j`
  startdoy=`echo $fn | awk '{print substr($1,1,3)}'`
  echo $startdoy", "$startyr

  fn=`date --date=$finaldy' '$finalmo' '$finalyr +%j`
  finaldoy=`echo $fn | awk '{print substr($1,1,3)}'`
  echo $finaldoy", "$finalyr

### WGET FILES:: ###

# Year loop:
  for (( year = $startyr ; year <= $finalyr; year++ )); do
     if [ -d $year ]; then
       echo " ... Year directory exists: "$year
     else
       echo " ... Making Year directory: "$year
       mkdir $year
     fi

   # DOY loop:
     for (( doy = $startdoy ; doy <= $finaldoy; doy++ )); do
        if [ $doy -le 9 ]; then
           doydir="00"$doy
        elif [ $doy -ge 10 -a $doy -le 99 ]; then
           doydir="0"$doy
        else
           doydir=$doy
        fi
        echo "... Doy: "$doydir", "$doy

     ## NLDAS-2 A files:
        file=${serverhost}/${sourcedira}/$year/$doydir/
        wget --no-check-certificate -nv -r -l1 -nH --reject "txt" --accept "grb" ${file} --cut-dirs=4  -a ${rundir}/download_nldas-2a.log

     ## NLDAS-2 B files:
        file=${serverhost}/${sourcedirb}/$year/$doydir/
        wget --no-check-certificate -nv -r -l1 -nH --reject "txt" --accept "grb" ${file} --cut-dirs=4  -a ${rundir}/download_nldas-2b.log
 
     done   # end doy loop
  done      # end year loop

### COMPLETED!! ###

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

###
