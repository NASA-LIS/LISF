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
 sourcedir=/lisdata_pub/data/TESTCASES/Public7.1/PMWSnowDA

  echo "== Downloading the LIS restart file == "

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "LIS_RST_NOAH33_201102282300.ens20.d01.nc" -l1 -nv ${file} -nH --cut-dirs=5 -a download_lisrestart.log 

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
