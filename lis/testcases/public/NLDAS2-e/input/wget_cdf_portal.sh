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
# Settings for wgetting CDF data files from the 
#  LIS Data Portal (https://portal.nccs.nasa.gov/lisdata_pub).
#
#  07.16.2014 -- James Geiger: v1.0: Initial version
#  10.31.2014 -- David Mocko:  v1.1: Public testcase for
#                                        LIS-7.0 Noah-3.6
#
### Necessary directory paths ###
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/TESTCASES/Public7.2/NLDAS2-a

  echo "== Downloading the CDF files for NLDAS2-a testcase == "

### WGET FILES:: ###

## Create directory structure:
   file="${serverhost}/${sourcedir}"

   wget --no-check-certificate -r --accept "Public-Noah-3.6_input_v7.tar.gz" -l1 -nv ${file} -nH --cut-dirs=4 -a download_cdf.log 

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "
  echo " -----  UNZIPPING FILES FOR:: " $sourcedir " -----  "

  cd NLDAS2-a
  rmdir --ignore-fail-on-non-empty */
  gzip -dc Public-Noah-3.6_input_v7.tar.gz | tar xf -

  echo " -----  COMPLETED UNZIPPING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
