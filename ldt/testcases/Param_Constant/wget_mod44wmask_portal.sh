#!/bin/bash
#
# Settings for wgetting MODIS 44w mask data file from the 
#  LIS Data Portal (http://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  03.05.2017 -- Kristi Arsenault: v1.1: Initial version
#
### Necessary directory paths ###

# Enter your destination directory (Default below matches LDT inputs files):
  targetdir="./input/"

  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/PARAMETERS/mask_parms/MOD44W_V5

  mkdir -p $targetdir
  cd $targetdir

## Create directory structure:
  file=${serverhost}/${sourcedir}

  wget --no-check-certificate -r --reject "index.html*","*gif","*css" -l2 -nv ${file} -nH --cut-dirs=3 -a download_maskfile.log

# Unpack files to target directory:
  cd mask_parms/MOD44W_V5
  gunzip global-1km.1gd4r.gz

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

##
