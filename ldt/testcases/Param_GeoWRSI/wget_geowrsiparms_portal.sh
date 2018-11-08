#!/bin/bash
#
#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center Land Information System (LIS) v7.1
#
# Copyright (c) 2015 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
#
# Settings for wgetting GeoWRSI NLDAS parameter files from the 
#  LIS Data Portal (http://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  07.05.2014 -- David Mocko:      v1.1: Public testcase for LIS
#  03.05.2017 -- David Mocko:      v1.2: Public testcase for LDT
#
### Necessary directory paths ###

# Enter your destination directory (Default below matches LDT inputs files):
  targetdir="./input/"

  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/PARAMETERS/GeoWRSI_PARAMS/data/Africa

  mkdir -p ${targetdir}
  cd ${targetdir}

  echo "== Downloading the parameter files for GeoWRSI == "

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --reject "index.html*","*gif","*css" -l2 -nv ${file} -nH --cut-dirs=3 -a download_geowrsiparms.log 

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
