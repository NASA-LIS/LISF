#!/bin/sh
#SBATCH --job-name=noah
#SBATCH --time=2:00:00
#SBATCH --account=NWP601
#SBATCH --output noah.slurm.out
##SBATCH --ntasks=504 --ntasks-per-socket=42 --ntasks-per-core=1
##SBATCH --ntasks=504 --ntasks-per-socket=8 --ntasks-per-core=1
#SBATCH --ntasks=256 --ntasks-per-socket=8 --ntasks-per-core=1
#SBATCH --cluster-constraint=blue
#SBATCH --exclusive
#SBATCH --mem=0
#------------------------------------------------------------------------------
#
# SCRIPT: run_lis_mr_noah39_rapid_hpc11.sh
#
# Batch script for running Noah39 and RAPID-based (and MERIT-based)
# MR LIS run.
#
# REVISION HISTORY:
# 05 Nov 2024: Eric Kemp, SSAI. Initial specification.
#------------------------------------------------------------------------------

# Avoid core files
export FI_VERBS_PREFER_XRC=0

ulimit -s unlimited

# When a batch script is started, it starts in the user's home directory.
# Change to the directory where job was submitted.
if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Environment
module use --append /ccs/home/emkemp/hpc11/privatemodules
module load lisf_7.6_prgenv_cray_8.5.0_cpe_23.12
module load afw-python/3.11-202406

lisconfig=lis.config.mr.noah39.rapid.galwem.17km.76

# Sanity checks
if [ ! -e ./$lisconfig ] ; then
    echo "ERROR, ./$lisconfig not found!"
    exit 1
fi

if [ ! -e ./LIS ] ; then
   echo "ERROR, ./LIS does not exist!" && exit 1
fi

srun ./LIS -f $lisconfig || exit 1

# The end
exit 0

