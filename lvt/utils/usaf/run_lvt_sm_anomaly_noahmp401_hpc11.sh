#!/bin/sh
#SBATCH --job-name=smanom
#SBATCH --time=1:00:00
#SBATCH --account nwp601
#SBATCH --output smanom.slurm.out
#SBATCH --ntasks=1
#SBATCH --cluster-constraint=blue
#SBATCH --exclusive
#SBATCH --mem=0

#------------------------------------------------------------------------------
#
# SCRIPT: run_lvt_sm_anomaly_noahmp401_hpc11.sh
#
# PURPOSE: Run LVT in a batch job on HPC11, to update 00Z NoahMP401 soil
# moisture climatology and to output climatology for the current month w/ soil
# moisture anomaly in netCDF format.  Assumes LVT has been compiled and
# lvt.config file has been customized.
#
# REVISION HISTORY:
# 05 Nov 2024: Eric Kemp, SSAI. Initial specification.
#
#------------------------------------------------------------------------------

if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Environment
module use --append /ccs/home/emkemp/hpc11/privatemodules
module load lisf_7.6_prgenv_cray_8.6.0_cpe_25.03_cce_19.0.0
module load afw-python/3.11-202511

if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

#lvtconfig=configs/lvt.config.template_sm_anomaly_noahmp401
lvtconfig=configs/lvt.config.sm_anomaly

if [ ! -e $lvtconfig ] ; then
   echo "ERROR, $lvtconfig does not exist!" && exit 1
fi

# LVT doesn't create RST directory before writing.  This is a work around.
if [ ! -e OUTPUT/STATS.sm_anomaly/RST ] ; then
    mkdir -p OUTPUT/STATS.sm_anomaly/RST || exit 1
fi

srun -n 1 ./LVT $lvtconfig || exit 1

exit 0
