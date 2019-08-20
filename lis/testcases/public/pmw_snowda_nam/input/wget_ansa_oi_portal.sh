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
#  07.17.2014 -- James Geiger: v1.0: Initial version
#
### Necessary directory paths ###
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/RS_DATA/ANSA_OI

  echo "== Downloading the ANSA snow depth observations == "

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "*README*" -l1 -nv ${file} -nH --cut-dirs=3 -a download_ansa_oi.log 
   wget --no-check-certificate -r --accept "*CITATION*" -l1 -nv ${file} -nH --cut-dirs=3 -a download_ansa_oi.log 

   file=${serverhost}/${sourcedir}/2011

   wget --no-check-certificate -r --reject "index.html*" --reject "*gif" -l1 -nv ${file} -nH --cut-dirs=3 -a download_ansa_oi.log 

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
