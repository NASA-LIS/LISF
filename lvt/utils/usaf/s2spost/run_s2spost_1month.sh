#!/bin/sh
#SBATCH --job-name=s2scf
#SBATCH --time=0:45:00
#SBATCH --account s1189
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint="sky|hasw|cas"
#------------------------------------------------------------------------------
#
# SCRIPT: run_s2spost_1month.sh
#
# DESCRIPTION: Batch job for processing one month of LIS output, producing
# CF-convention daily and monthly files.  Most work performed by Python
# scripts.
#
# USAGE: sbatch run_s2spost_1month.sh $CONFIGFILE $SCRIPTDIR $TOPDATADIR \
#           $YYYYMM $MODEL_FORCING
#         where $CONFIGFILE is path to s2spost config file
#               $SCRIPTDIR is path to s2spost scripts
#               $TOPDATADIR is top directory with LIS output
#               $YYYYMM is year and month of data to process
#               $MODEL_FORCING is ID for source of LIS atmospheric forcing
#
# REVISION HISTORY:
# 23 Sep 2021: Eric Kemp (SSAI), first version.
# 03 Nov 2021: Eric Kemp/SSAI, revised for new miniconda3 environment on
#   Discover, and to use s2spost config file.
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
module load lisf_7_intel_19_1_3_304_s2s

# Get the command line arguments.
if [ -z "$1" ] ; then
    echo "[ERR] Missing config file!" && exit 1
fi
CONFIGFILE=$1

if [ -z "$2" ] ; then
    echo "[ERR] Missing script directory" && exit 1
fi
SCRIPTDIR=$2

if [ -z "$3" ] ; then
    echo "[ERR] Missing top data directory with LIS output!" && exit 1
fi
TOPDATADIR=$3

if [ -z "$4" ] ; then
    echo "[ERR] Missing year and month (YYYYMM) to process!" && exit 1
fi
YYYYMM=$4

if [ -z "$5" ] ; then
    echo "[ERR] Missing ID of atmospheric forcing for LIS!" && exit 1
fi
MODEL_FORCING=$5

# Run the Python driver script for this month
if [ ! -e "$SCRIPTDIR/run_s2spost_1month.py" ] ; then
    echo "[ERR] $SCRIPTDIR/run_s2spost_1month.py does not exist!" && exit 1
fi
echo "[INFO] Processing $YYYYMM"
$SCRIPTDIR/run_s2spost_1month.py $CONFIGFILE $TOPDATADIR $YYYYMM \
                                 $MODEL_FORCING || exit 1


# The end
echo "[INFO] Finished CF processing of $YYYYMM"
exit 0
