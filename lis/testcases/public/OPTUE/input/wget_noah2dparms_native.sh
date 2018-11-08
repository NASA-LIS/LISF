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
#  Script:  wget_noah2dparms.sh
#
#  Description:  
#   Downloads all required 2-D base parameters for the Noah LSM.
#   When all files are fully downloaded/unzipped, the target directory
#    containing the files should be about ~2.8 GB.
#
#  Sources:
#    ftp://ftp.emc.ncep.noaa.gov
#    http://www.ral.ucar.edu/
#
#  Script Written by: K. Arsenault;  Aug 30, 2013
#  Script Written by: David Mocko;   July 3, 2014
# _________________________________________________________________________

  rundir=$PWD
# Enter your destination directory (Default below matches LDT inputs files):
  targetdir="./noah_2dparms/"

  echo "== Downloading the required Noah LSM 2D Parameters == "

# Change directory to target directory where files are to be downloaded to:
  mkdir -p $targetdir
  cd $targetdir
  echo "- Change to target directory: "$targetdir

# -------------------------------------------------------------------

# Obtain Modified IGBP MODIS 20-category vegetation (land-use) data:
#
#  The data come from Boston University, with some modifications by NCEP. 
#   NCEP has added new categories based on some other data sets, and has 
#   remapped certain categories to other indices. 

  echo " ... Getting the NCEP-modified MODIS-IGBP Landcover Map  "
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/ldas/noahlsm/igbp.bin.Z
  gunzip igbp.bin.Z

# Download README File:
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/ldas/noahlsm/README

# Obtain Hybrid STATSGO/FAO (30-second for CONUS /5-minute elsewhere) Soil Texture
# Download Main Data File (both top and bottom soil layers):
  echo " ... Downloading the Hybrid (NCAR) STATSGO/FAO Soil Texture Map  "
  wget --no-check-certificate http://www.ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/statsgo/topsoil30snew.gz
  wget --no-check-certificate http://www.ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/statsgo/botsoil30snew.gz
  gunzip topsoil30snew.gz
  gunzip botsoil30snew.gz
# Note: Only top layer currently used in LIS-Noah model runs.

# Download README for STATSGO+FAO soil texture File:
  wget --no-check-certificate http://www.ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/statsgo/READ_ME

# Obtain monthly greenness fraction (gfrac) and non-snow albedo source files:
  echo " ... Downloading the NCEP Greenness Fraction and Albedo Maps  "
  for eachmon in jan feb mar apr may jun jul aug sep oct nov dec
  do
    wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/test/fixed/albedo_${eachmon}.asc.Z
    wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/test/fixed/gfrac_${eachmon}.asc.Z
  done

# Obtain max and min gfrac:
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/test/fixed/gfrac_min.asc.Z
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/test/fixed/gfrac_max.asc.Z
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/test/fixed/gfrac_min_mon.asc.Z
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/test/fixed/gfrac_max_mon.asc.Z

# Obtain GFRAC-Albedo README file:
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/test/fixed/README_albedo_gfrac.txt

# Obtain quarterly (seasonal) non-snow albedo (global) values:
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/oper/fixed/albedo

# Obtain maximum snow albedo file and doc:
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/oper/fixed/maxsnoalb.asc.Z
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/oper/fixed/README_maxsnowalb.txt

# Obtain 1x1 deg ISLSCP-1 (non-elevation adjusted) bottom temperature (K) data:
  echo " ... Downloading the ISLSCP-1 Bottom Temperature Map  "
  wget --no-check-certificate https://ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/tbot/SOILTEMP.60
  wget --no-check-certificate https://ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/tbot/READ_ME
  mv READ_ME README.tbot_ISLSCP1 

# Obtain NCEP 1x1 deg slope type map values:
  echo " ... Downloading the NCEP Slope-type Map "
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/sfcflds/oper/fixed/islope

  gunzip *.Z

  echo "== Done downloading Noah LSM 2D parameter fields."

