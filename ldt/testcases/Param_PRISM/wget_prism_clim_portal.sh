#!/bin/bash
#
# Settings for wgetting PRISM data files from the 
#  LIS Data Portal (http://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  03.06.2017 -- Kristi Arsenault: v1.1: Updated
#
### Necessary directory paths ###

# Enter your destination directory:
  targetdir="./input/"

  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/PARAMETERS/climate_maps/PRISM/1KM

  mkdir -p ${targetdir}
  cd ${targetdir}

## Create directory structure:
   file=${serverhost}/${sourcedir}/

   wget --no-check-certificate -r --accept "*txt" -l1 -nv ${file} -nH --cut-dirs=5 -a download_prism1kmfiles.log

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

###
