#!/bin/bash
#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.3
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
#
# @mainpage NASA's Land Information System (LIS) NUOPC Coupling Test
# @author Daniel Rosen (daniel.rosen@noaa.gov)
# @author ESMF Support (esmf_support@ucar.edu)
# @date 08/10/2022 LIS NUOPC Cap Added to GitHub
#
# usage instructions
usage () {
  printf "Usage: $0 [OPTIONS]...\n"
  printf "\n"
  printf "OPTIONS\n"
  printf "  --system=SYSTEM\n"
  printf "      name of machine (e.g. 'discover', 'cheyenne'\n"
  printf "  --compiler=COMPILER\n"
  printf "      compiler to use; valid options are 'intel.X.Y.Z', \n"
  printf "      'gnu.X.Y.Z'; default is system dependent.\n"
  printf "\n"
}

# print settings
settings () {
  printf "Settings:\n"
  printf "\n"
  printf "  LISF_DIR=${LISF_DIR}\n"
  printf "  SYSTEM=${SYSTEM}\n"
  printf "  COMPILER=${MYCOMPILER}\n"
  printf "\n"
}

# find system name
find_system () {
    local sysname=`hostname`
    sysname="${sysname//[[:digit:]]/}"
    echo "$sysname"
}

# default settings
LISF_DIR="$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )")/../../../.." && pwd -P)"
SYSTEM=""
MYCOMPILER=""

# required arguments
if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
  usage
  exit 0
fi

# process arguments
while :; do
  case $1 in
    --help|-h) usage; exit 0 ;;
    --system=?*) SYSTEM=${1#*=} ;;
    --system) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --system=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --compiler=?*) MYCOMPILER=${1#*=} ;;
    --compiler) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    --compiler=) printf "ERROR: $1 requires an argument.\n"; usage; exit 1 ;;
    -?*) printf "ERROR: Unknown option $1\n"; usage; exit 1 ;;
    *) break
  esac
  shift
done

set -eu

# automatically determine system
if [ -z "${SYSTEM}" ] ; then
  SYSTEM=$(find_system)
fi

# automatically determine compiler
if [ -z "${MYCOMPILER}" ] ; then
  if [ "${SYSTEM}" = "discover" ]; then
    MYCOMPILER="lisf_7_intel_2021.4.0"
  else
    printf "ERROR: no default compiler for ${SYSTEM}\n"
    printf "\n"
    exit 1
  fi
fi

# print settings
settings

# load environment for this system/compiler combination
MODFILE="${LISF_DIR}/env/${SYSTEM}/${MYCOMPILER}"
if [ ! -f "${MODFILE}" ]; then
  printf "ERROR: module file does not exist for ${SYSTEM} ${MYCOMPILER}\n"
  printf "\n"
  exit 1
fi
. /etc/profile.d/modules.sh
module purge
module use ${LISF_DIR}/env/discover/
module load ${MYCOMPILER}
export ESMFMKFILE="${LIS_LIBESMF}/esmf.mk"
export WRF_HYDRO="1"

# configure
cd ${LISF_DIR}/lis && ./configure

exit 0
