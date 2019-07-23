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
# Settings for wgetting CLSM-F2.5 NLDAS parameter files from the 
#  LIS Data Portal (https://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  07.05.2014 -- David Mocko:      v1.1: Public testcase for
#                                        LIS-7.0 CLSM-F2.5
#
### Necessary directory paths ###
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/PARAMETERS/clsm25_parms/NLDAS_0.125-deg

  echo "== Downloading the parameter files for CLSM-F2.5 == "

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --reject "index.html*","*gif" -l2 -nv ${file} -nH --cut-dirs=2 -a download_clsmf25parms.log 

  cd PARAMETERS/clsm25_parms
  rmdir --ignore-fail-on-non-empty */

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
