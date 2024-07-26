#!/bin/bash
#
#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
#
# Settings for wgetting CDF data files from the 
#  LIS Data Portal (https://portal.nccs.nasa.gov/lisdata_pub).
#
#  07.16.2014 -- James Geiger: v1.0: Initial version
#
### Necessary directory paths ###
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/TESTCASES/Public7.2/FEWSNET

  echo "== Downloading the sample output files for FEWSNET-GeoWRSI2.0 testcase == "

### WGET FILES:: ###

## Create directory structure:
   file="${serverhost}/${sourcedir}"

   wget --no-check-certificate -r --accept "geowrsi2lsmtest_ldt_v7.tar.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a download_output.log 
   wget --no-check-certificate -r --accept "geowrsi2lsmtest_input_v7.tar.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a download_output.log 
   wget --no-check-certificate -r --accept "geowrsi2lsmtest_output_v72-dev.tar.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a download_output.log 

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "
  echo " -----  UNZIPPING FILES FOR:: " $sourcedir " -----  "

  gzip -dc geowrsi2lsmtest_ldt_v7.tar.gz | tar xf -
  gzip -dc geowrsi2lsmtest_input_v7.tar.gz | tar xf -
  gzip -dc geowrsi2lsmtest_output_v72-dev.tar.gz | tar xf -

  echo " -----  COMPLETED UNZIPPING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
