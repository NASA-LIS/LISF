#!/bin/sh
#SBATCH --job-name=s2sglb
#SBATCH --time=1:00:00
#SBATCH --account s1189
#SBATCH --output s2sglb.slurm.out
#SBATCH --ntasks=560
##SBATCH --mail-user=[user_email@address]
#SBATCH --mail-type=ALL
##SBATCH --qos=debug
#------------------------------------------------------------------------------
#
# SCRIPT: run_lis_global_usaf_forcing.sh
#
# DESCRIPTION: Batch script for running LIS in AGRMET Ops mode to generate
# USAF global atmospheric forcing for LIS S2S runs. LIS is run for a single day
# starting at 00Z.
#
# USAGE: sbatch run_lis_usaf_s2s_forcing.sh $YYYYMMDD
#         where $YYYYMMDD is the start date of the LIS run.
#
# REVISION HISTORY:
# 22 Sep 2021: Eric Kemp (SSAI), first version.
# 26 Sep 2021: Eric Kemp (SSAI), renamed to clarify forcing.
#------------------------------------------------------------------------------

ulimit -s unlimited

# When a batch script is started, it starts in user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Enviroment
module purge
unset LD_LIBRARY_PATH
## USER INPUTS:  ##
module use --append /home/emkemp/privatemodules
module load lisf_7.5_intel_2021.4.0_s2s

# Paths on local system
SCRIPTDIR=/discover/nobackup/projects/usaf_lis/GHI_S2S/GLOBAL/LISF_557WW_7.5/lis/utils/usaf/s2s/global_usaf_forc

CFGTMPL=input/lis.config.template
OUTDIR=./output
RSTDIR=./input/restarts
GRIBDIR=./output/grib

## END OF USER INPUTS ##

# Get the command line arguments.
if [ -z "$1" ] ; then
    echo "[ERR] Missing start date of LIS run!"
    exit 1
fi
YYYYMMDD=$1

# Customize lis.config file
if [ ! -e $SCRIPTDIR/customize_lis_config.py ] ; then
    echo "[ERR], $SCRIPTDIR/customize_lis_config.py does not exist!" && exit 1
fi
echo "[INFO] Customizing lis.config..."
$SCRIPTDIR/customize_lis_config.py $CFGTMPL $RSTDIR $YYYYMMDD || exit 1


# Run LIS
if [ ! -e ./LIS ] ; then
    echo "[ERR] ./LIS does not exist!" && exit 1
fi
echo "[INFO] Running LIS..."
mpirun -np $SLURM_NTASKS ./LIS || exit 1


# Clean up
if [ ! -e $SCRIPTDIR/store_lis_output.py ] ; then
    echo "[ERR] $SCRIPTDIR/store_lis_output.py does not exist!" && exit 1
fi
echo "[INFO] Saving LIS output..."
$SCRIPTDIR/store_lis_output.py $OUTDIR $RSTDIR $GRIBDIR $YYYYMMDD || exit 1


# The end
echo "[INFO] Completed LIS run!"
exit 0
