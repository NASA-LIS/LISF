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
# Settings for wgetting lai data files from the 
#  LIS Data Portal (https://portal.nccs.nasa.gov/lisdata_pub).
#
#  06.27.2014 -- Kristi Arsenault: v1.0: Initial version
#  07.04.2014 -- David Mocko:      v1.1: Public testcase for
#                                        LIS-7.0 Noah-3.3
#
### Necessary directory paths ###
 resolution=25KM
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/PARAMETERS/UMD
# user=username
# pass=password

  echo "== Downloading the LIS lai files needed for RTM  == "

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}/${resolution}

   wget --no-check-certificate -r --accept "lai_AVHRR*" -l1 -nv ${file} -nH --cut-dirs=2 -a download_lisnoahparms.log # --http-user=$user --http-password=$pass

   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "*README*" -l1 -nv ${file} -nH --cut-dirs=3 -a download_lisnoahparms.log # --http-user=$user --http-password=$pass

### COMPLETED!! ###

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir/${resolution} " -----  "
  echo " -----  UNZIPPING FILES FOR:: " $sourcedir/${resolution} " -----  "

  cd PARAMETERS/UMD
  rmdir --ignore-fail-on-non-empty */
  cd ${resolution}
  gunzip *.gz

  echo " -----  COMPLETED UNZIPPING FILES FOR:: " $sourcedir/${resolution} " -----  "

#  exit 0

###
