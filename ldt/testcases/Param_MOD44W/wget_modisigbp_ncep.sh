#!/bin/sh
#
# Modified IGBP MODIS 20-category vegetation (land-use) data
#
# The data come from Boston University, with some modifications by NCEP. 
#  NCEP has added new categories based on some other data sets, and has 
#  remapped certain categories to other indices.  Information from NCEP 
#  regarding their recategorization.
#
# Additional information from Boston University ...
# ---------------------------------------------------------------------

# Enter your destination directory (Default below matches LDT inputs files):
  targetdir="./input/"

  mkdir -p $targetdir
  cd $targetdir

 echo "== Downloading the NCEP-modified MODIS-IGBP Landcover Map == "

# Download Main Data File:
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/ldas/noahlsm/igbp.bin.Z -a download_igbp.log
  gunzip igbp.bin.Z

# Download README File:
  wget ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/ldas/noahlsm/README -a download_igbp.log
  mv README README_igbp.txt

 echo " Done! "
