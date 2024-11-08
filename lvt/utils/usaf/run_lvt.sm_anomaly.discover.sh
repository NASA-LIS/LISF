#!/bin/sh
#SBATCH --job-name=smanom
#SBATCH --time=1:00:00
#SBATCH --account s1189
#SBATCH --output smanom.slurm.out
#SBATCH --ntasks=1 --constraint="[mil]"
##SBATCH --qos=debug
#------------------------------------------------------------------------------
#
# SCRIPT: run_lvt.sm_anomaly.discover.sh
#
# PURPOSE: Run LVT in a batch job on Discover, to update 00Z soil moisture
# climatology and to output climatology for the current month w/ soil moisture
# anomaly in netCDF format.  Assumes LVT has been compiled and lvt.config file
# has been customized.
#
# REVISION HISTORY:
# 28 Jun 2021: Eric Kemp (SSAI), first version.
# 15 Oct 2024: Eric Kemp (SSAI), updated for Milan nodes.
#
#------------------------------------------------------------------------------

if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

module purge
module use --append /home/emkemp/privatemodules/sles15
module load lisf_7.6_intel_2023.2.1_emk

if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

if [ ! -e configs/lvt.config.sm_anomaly ] ; then
   echo "ERROR, configs/lvt.config.sm_anomaly does not exist!" && exit 1
fi

# LVT doesn't create RST directory before writing.  This is a work around.
if [ ! -e OUTPUT/STATS.sm_anomaly/RST ] ; then
    mkdir -p OUTPUT/STATS.sm_anomaly/RST || exit 1
fi

time mpirun -np 1 ./LVT configs/lvt.config.sm_anomaly || exit 1

exit 0
