#!/bin/bash
## REQUIRED SET UP:
#SBATCH --job-name=s2sglb
#SBATCH --output s2sglb.slurm.out
#SBATCH --time=2:00:00
#
## USER INPUTS HERE:
### USING SRUN ON DISCOVER
## Discover
#SBATCH --ntasks=480
#SBATCH --ntasks-per-socket=12
#SBATCH --ntasks-per-core=1
#SBATCH --account s1189
#SBATCH --constraint="mil"
#SBATCH --mail-user=[USERNAME]
#SBATCH --mail-type=ALL
#
### Cray -- 480 tasks
##SBATCH  -N 40
##SBATCH --ntasks-per-node=12
##SBATCH --cluster-constraint=blue
##SBATCH --partition=batch
##SBATCH --exclusive
##SBATCH --mem=0
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
# 20 May 2024: K. Arsenault (SAIC), updated timeframe and case.
# 17 Nov 2025: K. Arsenault (SAIC), updated timeframe and case for LISV7.7.
#------------------------------------------------------------------------------

# Paths on local system
SCRIPTDIR=/discover/nobackup/projects/ghilis/S2S/GLOBAL/use-cases/LISV7.7/S2S_Daily_Forc/

CFGTMPL=./input/lis.config.template
OUTDIR=./output
RSTDIR=./input/restarts
GRIBDIR=./output/grib

## END OF USER INPUTS ##

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
  export I_MPI_PMI_LIBRARY=/usr/slurm/lib64/libpmi2.so
  export I_MPI_PMI_VALUE_LENGTH_MAX=480
else
  export I_MPI_PMI_VALUE_LENGTH_MAX=480
fi

## USER INPUTS:  ##
if [[ $NODE_NAME =~ discover* ]] || [[ $NODE_NAME =~ borg* ]]; then
  module use --append ~/privatemodules
  module --ignore-cache load lisf_7.5_intel_2023.2.1_s2s
else
  # e.g., Cray environment -- Latest for LISV7.7 S2S
  module use --append ~/privatemodules
  module load lisf_7.6_prgenv_cray_8.6.0_cpe_25.03_cce_19.0.0_s2s
fi
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
if [[ $NODE_NAME =~ discover* ]] || [[ $NODE_NAME =~ borg* ]]; then
   echo "[INFO] -- Running on Discover -- "
   srun --mpi=pmi2 --ntasks=$SLURM_NTASKS \
         --ntasks-per-socket=$SLURM_NTASKS_PER_SOCKET \
         --ntasks-per-core=$SLURM_NTASKS_PER_CORE \
         --cpu-bind="none" \
         ./LIS || exit 1
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
