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
# Settings for wgetting Real-time (RT) SMOPS data files from the 
#  LIS Data Portal (https://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  07.04.2014 -- David Mocko:      v1.1: Public testcase for
#                                        LIS-7.0 Noah-3.3
#
### Necessary directory paths ###

 targetdir=./RT_SMOPS

 rundir=$PWD

 mkdir -p $targetdir
 cd $targetdir

 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/RS_DATA/RT_SMOPS

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --reject "index.html*" --reject "*gif" -l1 -nv ${file} -nH --cut-dirs=4 -a ${rundir}/download_rtsmop.log 

   tar -xzvf *tgz >> ${rundir}/download_rtsmop.log

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

###
