#!/bin/sh
#SBATCH --job-name=s2sglb
#SBATCH --time=0:30:00
#SBATCH --account s1189
#SBATCH --output s2sglb.slurm.out
#SBATCH --ntasks=560
#SBATCH --mail-user=eric.kemp@nasa.gov
#SBATCH --mail-type=ALL
##SBATCH --qos=debug
#------------------------------------------------------------------------------
#
# SCRIPT: run_lis_global_s2s_forcing.sh
#
# DESCRIPTION: Batch script for running LIS in AGRMET Ops mode to generate
# global atmospheric forcing for LIS S2S runs. LIS is run for a single day
# starting at 00Z.
#
# USAGE: sbatch run_lis_global_s2s_forcing.sh $YYYYMMDD
#         where $YYYYMMDD is the start date of the LIS run.
#
# REVISION HISTORY:
# 22 Sep 2021: Eric Kemp (SSAI), first version.
#
#------------------------------------------------------------------------------

ulimit -s unlimited

# When a batch script is started, it starts in user's home directory.
# Changet to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Enviroment
module purge
unset LD_LIBRARY_PATH
module use --append /home/emkemp/privatemodules
module load lisf_7_intel_19_1_3_304

# Paths on local system
SCRIPTDIR=/discover/nobackup/projects/lis_aist17/emkemp/AFWA/lis74_reanalysis_test/LISF/lis/configs/557WW-7.4-FOC/S2S_GLOBAL
CFGTMPL=$SCRIPTDIR/lis.config.template
OUTDIR=./output
RSTDIR=./restarts
GRIBDIR=./grib

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
