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
# Settings for wgetting UMD data files from the 
#  LIS Data Portal (http://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  03.02.2017 -- Kristi Arsenault: v2.0: Updated for lispub
#
### Necessary directory paths ###

# Enter your destination directory (Default below matches LDT inputs files):
  targetdir="./input/"

  resolution=10KM
  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/PARAMETERS/UMD

  echo "== Downloading the LDT parameter files for UMD Landcover  == "

  mkdir -p ${targetdir}
  cd ${targetdir}

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}/${resolution}

   wget --no-check-certificate -r --accept "*landcover*UMD*" --reject "*.tgz" -l1 -nv ${file} -nH --cut-dirs=3 -a download_umdlcparm.log

   cd UMD/${resolution}
   gunzip *.gz

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: "$sourcedir/${resolution}" -----  "

###
