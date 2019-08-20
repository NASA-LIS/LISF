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
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  07.04.2014 -- David Mocko:      v1.1: Public testcase for
#                                        LIS-7.0 Noah-3.3
#  07.17.2014 -- James Geiger:     v1.2: Public testcase for PMW Snow DA
#
### Necessary directory paths ###
 rundir=$PWD
 resolution=5KM
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/PARAMETERS/UMD

  echo "== Downloading the LIS parameter files for Noah LSM == "

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}/${resolution}

   wget --no-check-certificate -r --accept "*NCEP*","*FAO*","elev_GTOPO30.1gd4r","*landcover*UMD*","mxsnoalb_MODIS.1gd4r" -l1 -nv ${file} -nH --cut-dirs=2 -a $rundir/download_lisnoahparms.log 

   cd PARAMETERS
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "*README*","*CITATION*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisnoahparms.log 


## Clean up:
   cd UMD
   rmdir --ignore-fail-on-non-empty */

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir/${resolution} " -----  "

#  exit 0

###
