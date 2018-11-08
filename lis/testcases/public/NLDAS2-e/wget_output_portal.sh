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
#
### Necessary directory paths ###
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/TESTCASES/Public7.2/NLDAS2-a

  echo "== Downloading the sample output files for NLDAS2-a testcase == "

### WGET FILES:: ###

## Create directory structure:
   file="${serverhost}/${sourcedir}"

   wget --no-check-certificate -r --accept "Public-Noah-3.6-NLDAS-OpenLoop_output_v71.tar.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a download_output.log 

   wget --no-check-certificate -r --accept "Public-Noah-3.6-NLDAS-DA_output_v71.tar.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a download_output.log 

   wget --no-check-certificate -r --accept "Public-Noah-3.6-HYMAP-OpenLoop_output_v71.tar.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a download_output.log 

   wget --no-check-certificate -r --accept "Public-Noah-3.6-HYMAP-DA_output_v71.tar.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a download_output.log 

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "
  echo " -----  UNZIPPING FILES FOR:: " $sourcedir " -----  "

  gzip -dc Public-Noah-3.6-NLDAS-OpenLoop_output_v71.tar.gz | tar xf -
  gzip -dc Public-Noah-3.6-NLDAS-DA_output_v71.tar.gz | tar xf -
  gzip -dc Public-Noah-3.6-HYMAP-OpenLoop_output_v71.tar.gz | tar xf -
  gzip -dc Public-Noah-3.6-HYMAP-DA_output_v71.tar.gz | tar xf -

  echo " -----  COMPLETED UNZIPPING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
