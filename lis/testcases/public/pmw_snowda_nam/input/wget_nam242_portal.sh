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
# Settings for wgetting GDAS data files from the 
#  LIS Data Portal (https://portal.nccs.nasa.gov/lisdata_pub).
#
#  02.07.2012 -- Kristi Arsenault: v1.0: Initial version
#  06.20.2012 -- Kristi Arsenault: v1.1: Fixed download bug and excess file removal
#  07.17.2014 -- James Geiger: v1.2: Modified for NAM242 files
#
### Necessary directory paths ###
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/MET_FORCING/NAM242

 startyr=2011
 finalyr=2011

### WGET FILES:: ###

  for (( year = $startyr ; year <= $finalyr; year++ )); do
    #for mon in 01 02 03 04 05 06 07 08 09 10 11 12; do
    for mon in 03 04 05; do
       for day in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31; do
          for hr in 00 06 12 18; do

          ## Create directory structure:
          yyyymmddhr=$year$mon$day/$hr
          file=${serverhost}/${sourcedir}/${yyyymmddhr}

          wget --no-check-certificate -r --reject "index.html*" --reject "*gif" --accept "fh*" -l1 -nv ${file} -nH --cut-dirs=2 -a download_nam242.log

          done  # end hr loop
       done  # end day loop
    done  # end month loop
  done    # end year loop


### COMPLETED!! ###

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "
  exit 0

###

