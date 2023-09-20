#!/bin/bash
# REQUIRED SET UP:
#SBATCH --job-name=s2sglb
#SBATCH --ntasks=560
#SBATCH --ntasks-per-core=1
#SBATCH --time=1:00:00
#SBATCH --output s2sglb.slurm.out
#
## USER INPUTS HERE:
## Discover
#SBATCH --account s1189
#SBATCH --constraint="sky|cas"
#SBATCH --mail-user=[USERNAME]
#SBATCH --mail-type=ALL
#
### Cray
##SBATCH --account=
##SBATCH --cluster-constraint=green
##SBATCH --partition=batch
#
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
# 07 Mar 2022: S. Mahanama (SAIC), added uname option for different systems.
# 08 Mar 2023: K. Arsenault (SAIC), additional modifications for Cray env.
#------------------------------------------------------------------------------

export NODE_NAME=`uname -n`
ulimit -s unlimited

# When a batch script is started, it starts in user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Enviroment
if [[ $NODE_NAME =~ discover* ]] || [[ $NODE_NAME =~ borg* ]]; then
  module purge
  unset LD_LIBRARY_PATH
fi

## USER INPUTS:  ##
if [[ $NODE_NAME =~ discover* ]] || [[ $NODE_NAME =~ borg* ]]; then
  module use --append ~/privatemodules
  module load lisf_7.5_intel_2021.4.0_s2s 
# e.g., Cray environment
else
  module load lisf_7.5_prgenv_cray_8.3.3_s2s
fi

# Paths on local system
SCRIPTDIR=[local_script_path]

CFGTMPL=./input/lis.config.template
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
python $SCRIPTDIR/customize_lis_config.py $CFGTMPL $RSTDIR $YYYYMMDD || exit 1

# Run LIS
if [ ! -e ./LIS ] ; then
    echo "[ERR] ./LIS does not exist!" && exit 1
fi
echo "[INFO] Running LIS..."
if [[ $NODE_NAME =~ discover* ]] || [[ $NODE_NAME =~ borg* ]]; then
    mpirun -np $SLURM_NTASKS ./LIS || exit 1
else
   /usr/bin/time srun --cpu-bind=rank_ldom ./LIS || exit 1
fi

# Clean up
if [ ! -e $SCRIPTDIR/store_lis_output.py ] ; then
    echo "[ERR] $SCRIPTDIR/store_lis_output.py does not exist!" && exit 1
fi
echo "[INFO] Saving LIS output..."
python $SCRIPTDIR/store_lis_output.py $OUTDIR $RSTDIR $GRIBDIR $YYYYMMDD || exit 1

# The end
echo "[INFO] Completed LIS run!"
exit 0
