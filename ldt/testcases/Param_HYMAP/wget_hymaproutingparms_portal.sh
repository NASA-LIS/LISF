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
# Settings for wgetting HYMAP routing data files from the 
#  LIS Data Portal (http://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  07.04.2014 -- David Mocko:      v1.1: Public testcase for LIS-7.0
#  03.04.2017 -- Kristi Arsenault: v1.2: Public testcase for LDT 7.1
#
### Necessary directory paths ###

# Enter your destination directory (Default below matches LDT inputs files):
  targetdir="./input/"

  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/PARAMETERS/HYMAP_parms/NLDAS_12.5KM

  mkdir -p ${targetdir}
  cd ${targetdir}

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept bin,"*hydrographfiles*" -l1 -nv ${file} -nH --cut-dirs=3 -a download_hymaproutingparms.log 

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "
#  exit 0

###
