#!/bin/bash
#
# Settings for wgetting NLDAS data files from the 
#  LIS Data Portal (http://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  03.06.2014 -- Kristi Arsenault: v1.1: Updated
#
### Necessary directory paths ###

# Enter your destination directory:
  targetdir="./input/"

  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/PARAMETERS/NLDAS_0.125/
  rundir=$PWD

  mkdir -p ${targetdir}
  cd ${targetdir}

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}/NARR_elevation.1gd4r
   wget --no-check-certificate -r --reject "index.html*" --reject "*gif" -l1 -nv ${file} -nH --cut-dirs=5 -a download_log.file 
   file=${serverhost}/${sourcedir}/NARR_elev-diff.1gd4r
   wget --no-check-certificate -r --reject "index.html*" --reject "*gif" -l1 -nv ${file} -nH --cut-dirs=5 -a download_log.file 
   file=${serverhost}/${sourcedir}/NLDAS_elevation.1gd4r
   wget --no-check-certificate -r --reject "index.html*" --reject "*gif" -l1 -nv ${file} -nH --cut-dirs=5 -a download_log.file 

### COMPLETED!! ###

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

###
