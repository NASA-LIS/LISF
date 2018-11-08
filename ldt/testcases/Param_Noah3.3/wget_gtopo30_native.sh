#!/bin/sh
#
#  Script:  wget_gtopo30_native.sh
#
#  Description:  
#   Downloads all required 2-D tiled-based GTOPO30 elevation files
#   When all files are fully downloaded/unzipped, the target directory
#    containing the files should be about ~2.8 GB.
#
#  Sources:
#    http://edcftp.cr.usgs.gov/pub/data/gtopo30/global/
#  GTOPO30 file date:  2/28/1997
#
#  Script Written by: K. Arsenault;  Aug 30, 2013
#  Script updated by: K. Arsenault;  Mar 10, 2017
# ____________________________________________________________________

# Enter your destination directory (Default below matches LDT inputs files):
  targetdir=./input/

  echo "== Downloading the GTOPO-30sec Elevation Tile Files == "

# Change directory to target directory where files are to be downloaded to:
  mkdir -p ${targetdir}
  cd ${targetdir}

  echo "- Change to target directory: "${targetdir}

## Create directory structure:
  serverhost=http://portal.nccs.nasa.gov
  sourcedir=/lisdata_pub/data/PARAMETERS/topo_parms/GTOPO30_native
  file=${serverhost}/${sourcedir}

  wget --no-check-certificate -r --accept "*.gz" -l1 -nv ${file} -nH --cut-dirs=4 -a download_gtopo30_native.log

# Loop over each gridded tile and download *zip file:
  cd GTOPO30_native
  for nstile in n90 n40 s10 s60; do
    for wetile in w180 w140 w120 w100 w060 w020 e020 e060 e100 e140; do
       tar -xzvf ${wetile}${nstile}.tar.gz
    done
  done

  echo "== Done downloading GTOPO30 tile fields."

