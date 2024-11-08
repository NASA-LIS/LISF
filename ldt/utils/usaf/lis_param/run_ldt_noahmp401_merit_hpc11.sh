#!/bin/sh
#SBATCH --job-name=ldt
#SBATCH --time=1:00:00
#SBATCH --account=NWP601
#SBATCH --output ldt.slurm.out
#SBATCH --ntasks=1
#SBATCH --cluster-constraint=blue
#SBATCH --exclusive
#SBATCH --mem=0
#------------------------------------------------------------------------------
#
# SCRIPT: run_ldt_noahmp401_merit_hpc11.sh
#
# Batch script for running LDT to generate MERIT-based LIS parameter file
# for NoahMP401 LSM on HPC11. Output is used by NoahMP401-based NRT-Streamflow
# and MR LIS runs.
#
# REVISION HISTORY:
# 05 Nov 2024: Eric Kemp, SSAI. Initial specification.
#------------------------------------------------------------------------------


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

ldtconfig=ldt.config.foc.nrt_sf.noahmp401.param.76

# Sanity checks
if [ ! -e ./$ldtconfig ] ; then
    echo "ERROR, ./$ldtconfig not found!" && exit 1
fi
if [ ! -e ./LDT ] ; then
    echo "ERROR, ./LDT does not exist!" && exit 1
fi

echo `date`
srun ./LDT $ldtconfig || exit 1
echo `date`

# The end
exit 0
