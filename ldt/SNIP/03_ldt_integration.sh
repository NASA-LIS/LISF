#!/bin/bash
#SBATCH --job-name=03_ldt_integration
#SBATCH --output=./log/03_ldt_integration_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --constraint="mil"
#SBATCH --mail-type=ALL
#SBATCH --account=s1189
#SBATCH --qos=debug

# load modules
module purge
unset LD_LIBRARY_PATH
module use --append /home/emkemp/privatemodules/sles15
module load lisf_7.6_intel_2023.2.1_emk_aiml

# Ensure we are in the correct directory
cd $SLURM_SUBMIT_DIR

TARGET_DATETIME="$1"
TARGET_DATETIME_10="${TARGET_DATETIME:0:10}"
LDT_CONFIG="./03_ldt_integration.config"

if [ ! -f "$LDT_CONFIG" ]; then
  echo "ERROR: LDT config file not found: $LDT_CONFIG"
  exit 1
fi

# 1. Replace the Date (Using the robust .* regex and the correct target file)
sed -i "s/^SNIP valid date (YYYYMMDDHH):.*/SNIP valid date (YYYYMMDDHH):                    ${TARGET_DATETIME_10}/g" "$LDT_CONFIG"

# 2. Run LDT (single process here)
echo "Starting LDT post-processing for ${TARGET_DATETIME_10}..."
mpirun -np 1 ./LDT "$LDT_CONFIG"

# 3. Check for success
if [ $? -ne 0 ]; then
  echo "ERROR: LDT post-processing crashed!"
  exit 1
fi

echo "LDT post-processing complete"