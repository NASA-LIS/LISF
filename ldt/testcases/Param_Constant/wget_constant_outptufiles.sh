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
# Settings for wgetting MODIS44w output data files from the 
#  LIS Data Portal (http://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  03.02.2017 -- Kristi Arsenault: v2.0: Updated for lispub
#
### Necessary directory paths ###
  rundir=$PWD

# Enter your destination directory (Default below matches LDT inputs files):
  targetdir="./input/"

  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/TESTCASES/LDT/Param_Constant

  echo "== Downloading the LDT Constant Mask Output Files == "

  mkdir -p ${targetdir}
  cd ${targetdir}

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

 # Download the original output files to compare with:
   cd ${rundir}
   wget --no-check-certificate -r --accept "*.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a ${targetdir}/download_constant_output.log

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: "$sourcedir/${resolution}" -----  "

###
