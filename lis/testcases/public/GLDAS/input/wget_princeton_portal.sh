#!/bin/bash
#
#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center Land Information System (LIS) v7.2
#
# Copyright (c) 2015 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
#
# Settings for wgetting UMD data files from the 
#  LIS Data Portal (https://portal.nccs.nasa.gov/lisdata_pub).
#
#  07.16.2014 -- James Geiger: v1.0: Initial version
#
### Necessary directory paths ###
 rundir=$PWD
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=lisdata_pub/data/MET_FORCING/PRINCETON

  echo "== Downloading Princeton forcing data == "

### WGET FILES:: ###

## Create directory structure:
   targetdir=$rundir/MET_FORCING
   mkdir -p $targetdir
   cd $targetdir

   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "*README*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_princeton.log 
   wget --no-check-certificate -r --accept "*CITATION*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_princeton.log 

   file="${serverhost}/${sourcedir}/1948"
   wget --no-check-certificate -r --accept "*.nc" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_princeton.log 

### COMPLETED!! ###

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
