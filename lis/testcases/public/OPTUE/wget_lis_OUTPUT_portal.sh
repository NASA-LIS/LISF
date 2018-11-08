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
#
### Necessary directory paths ###
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/TESTCASES/Public7.2/OPTUE
# user=username
# pass=password

  echo "== Downloading the OUTPUT data for OPT/UE Public Testcase  == "

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "Public-OPTUE_output_v72-dev.tar.gz" -l1 -nv ${file} -nH --cut-dirs=5 -a download_amsre_obs.log # --http-user=$user --http-password=$pass

### COMPLETED!! ###

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "
  echo " -----  UNZIPPING FILES FOR:: " $sourcedir " -----  "

  gzip -dc Public-OPTUE_output_v72-dev.tar.gz | tar xf -
  #rm Public-OPTUE_output_v72.tar.gz

  echo " -----  COMPLETED UNZIPPING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
