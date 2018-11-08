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
#  USGS 30-second global 24-category vegetation (land-use) data
#
#  03.02.2017 -- Kristi Arsenault: v1.0: Updated for lispub
#
### Necessary directory paths ###

# Enter your destination directory (Default below matches LDT inputs files):
  targetdir="./input/"

  mkdir -p ${targetdir}
  cd ${targetdir}
#
# ---------------------------------------------------------------------

 echo "== Downloading the USGS 24-category Landcover Map == "

# Download Main Data File:
  wget http://www.ral.ucar.edu/research/land/technology/lsm/sfc_fields/USGS/veg30susgs.gz
  gunzip veg30susgs.gz

# Download README File:
  wget http://www.ral.ucar.edu/research/land/technology/lsm/sfc_fields/USGS/READ_ME

# Download Figure:
  wget http://www.ral.ucar.edu/research/land/technology/lsm/sfc_fields/USGS/VEG.gif

 echo " Done! "
 
