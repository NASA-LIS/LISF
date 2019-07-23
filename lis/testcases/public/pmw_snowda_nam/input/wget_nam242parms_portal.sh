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
#  07.17.2014 -- James Geiger: v1.0: Initial version
#
### Necessary directory paths ###
 rundir=$PWD
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=/lisdata_pub/data/PARAMETERS/metforcing_parms/NAM

  echo "== Downloading the NAM242 parameter files == "

### WGET FILES:: ###

## Create directory structure:
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "terrain.242.grb" -l1 -nv ${file} -nH --cut-dirs=2 -a $rundir/download_nam242parms.log 

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir/${resolution} " -----  "

#  exit 0

###
