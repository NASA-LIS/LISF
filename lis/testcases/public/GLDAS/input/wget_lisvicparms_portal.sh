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
#  07.16.2014 -- James Geiger: v1.0: Initial version
#
### Necessary directory paths ###
 rundir=$PWD
 serverhost=https://portal.nccs.nasa.gov
 sourcedir=lisdata_pub/data/PARAMETERS/UMD/100KM

  echo "== Downloading LIS/VIC input parameter data == "

### WGET FILES:: ###

## Create directory structure:
   targetdir=$rundir/LS_PARAMETERS
   mkdir -p $targetdir
   cd $targetdir

   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "*README*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 
   wget --no-check-certificate -r --accept "*CITATION*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 
   wget --no-check-certificate -r --accept "elev_GTOPO30.1gd4r" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 

   cd $rundir
   #targetdir=$rundir/metforcing_parms
   #mkdir -p $targetdir
   #cd $targetdir

   sourcedir=lisdata_pub/data/PARAMETERS/metforcing_parms/PRINCETON
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "*README*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 
   wget --no-check-certificate -r --accept "*CITATION*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 
   wget --no-check-certificate -r --accept "hydro1k_elev_mean_1d.asc" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 

   sourcedir=lisdata_pub/data/PARAMETERS/VIC_PARAMETERS/1deg 
   file=${serverhost}/${sourcedir}

   wget --no-check-certificate -r --accept "*README*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 
   wget --no-check-certificate -r --accept "*CITATION*" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 
   wget --no-check-certificate -r --accept "g1d.landcover.1gd4r","snowband_1.0_25","soilparam_1.0","vegparam_1.0","vic_global_mask_1.0.bin","world_veg_lib.txt" -l1 -nv ${file} -nH --cut-dirs=3 -a $rundir/download_lisvicparms.log 
### COMPLETED!! ###

  echo " -----  COMPLETED DOWNLOADING FILES FOR:: " $sourcedir " -----  "

#  exit 0

###
