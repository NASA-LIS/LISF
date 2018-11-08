#!/bin/bash
#
# Settings for wgetting GDAS data files from the 
#  LIS Data Portal (http://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  03.06.2014 -- Kristi Arsenault: v1.1: Updated
#  03.06.2017 -- Kristi Arsenault: v1.2: Updated
#
### Necessary directory paths ###

# Enter your destination directory:
  targetdir="./input/"

  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/PARAMETERS/metforcing_parms/GDAS

  mkdir -p ${targetdir}
  cd ${targetdir}

## Create directory structure:
   file=${serverhost}/${sourcedir}/

   wget --no-check-certificate -r --accept "*global*","*README*","*CITATION*" -l1 -nv ${file} -nH --cut-dirs=5 -a download_gdaselevfiles.log

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

###
