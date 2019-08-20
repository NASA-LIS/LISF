#!/bin/sh
#
#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center Land Information System (LIS) v7.2
#
# Copyright (c) 2015 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
# 
#  Script:  wget_native_srtm30.sh
#
#  Description:  
#   Downloads all required 2-D tiled-based 30arcsec SRTM v2.1 elevation files
#   When all files are fully downloaded/unzipped, the target directory
#    containing the files should be about ~2.8 GB.
#
#  Sources:
#    http://dds.cr.usgs.gov/srtm/version2_1/SRTM30/
#    File Date: 8/23/2010
#
#  Script Written by: K. Arsenault;  Aug 30, 2013
#  07.04.2014 -- David Mocko:      v1.1: Public testcase for LIS-7.0
#
################################################################

  rundir=$PWD
  targetdir=./PARAMETERS/topo_parms/SRTM

  echo "== Downloading the SRTM-30sec Elevation Tile Files == "

# Change directory to target directory where files are to be downloaded to:
  mkdir -p $targetdir
  cd $targetdir
  echo "- Change to target directory: "$targetdir

# Loop over each gridded tile and download *zip file:
  for nstile in n90 n40 s10; do
    for wetile in w180 w140 w100 w060 w020 e020 e060 e100 e140; do
       wget --no-check-certificate http://dds.cr.usgs.gov/srtm/version2_1/SRTM30/${wetile}${nstile}/${wetile}${nstile}.dem.zip -nv -a ${rundir}/download_srtm30native.log
       unzip ${wetile}${nstile}.dem.zip >> ${rundir}/download_srtm30native.log
    done
  done

# Obtain SRTM30 documentation and version release info:
  wget --no-check-certificate http://dds.cr.usgs.gov/srtm/version2_1/SRTM30/srtm30_documentation.pdf -nv -a ${rundir}/download_srtm30native.log
  wget --no-check-certificate http://dds.cr.usgs.gov/srtm/version2_1/SRTM30/srtm30_version_history.pdf -nv -a ${rundir}/download_srtm30native.log

 echo "== Done downloading SRTM30 tile fields."
