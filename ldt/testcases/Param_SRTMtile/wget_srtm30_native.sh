#!/bin/sh
#
#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center Land Information System (LIS) v7.1
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
#  Updated: 07.04.2014 -- v1.1, David Mocko: Public testcase 
#                          for LIS-7.0
#  Updated: 03.04.2017 -- v1.2, K. Arsenault; latest LDT release
# ____________________________________________________________

  rundir=$PWD
  targetdir=./input/

  echo "== Downloading the SRTM-30sec Elevation Tile Files == "

# Change directory to target directory where files are to be downloaded to:
  mkdir -p $targetdir
  mkdir -p ${targetdir}/SRTM30
  cd $targetdir/SRTM30
  echo "- Change to target directory: "${targetdir}/SRTM30

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
