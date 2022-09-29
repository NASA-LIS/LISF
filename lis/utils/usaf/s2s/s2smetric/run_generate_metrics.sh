#!/bin/sh
#SBATCH --job-name=s2s_metrics
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --time=03:00:00
#SBATCH --account=s1189
#SBATCH --constraint="sky|hasw|cas"
#SBATCH --mail-type=END
#SBATCH --output=s2s_metrics-%j.out
#SBATCH --error=s2s_metrics-%j.err
#------------------------------------------------------------------------------
#
# SCRIPT: run_generate_metrics.sh
#
# PURPOSE: Batch job script for calculating S2S metrics.  Based on
# run_Convert_Dyn_FCST_postproc.scr
#
# REVISION HISTORY:
# 25 Oct 2021: Eric Kemp/SSAI, first version.
# 30 Oct 2021: Eric Kemp/SSAI, revised to leverage s2smetric config file.
#
#------------------------------------------------------------------------------

# Read command line arguments.
SCRIPT=${1}
FCST_INIT_MON=${2}
TARGET_YEAR=${3}
MODEL=${4}
CONFIGFILE=${5}
PYDIR=${6}

# Change to directory with Python script.  FIXME: Add script directory
# as an argument.
cd $PYDIR || (echo "[ERR] Cannot change to directory $PYDIR" && exit 1)

# Set up environment
ulimit -s unlimited
source /usr/share/modules/init/sh
module use --append /home/emkemp/privatemodules
module load lisf_7_intel_19_1_3_304_s2s

# Execute Python script
python3 $SCRIPT $FCST_INIT_MON $TARGET_YEAR $MODEL $CONFIGFILE


