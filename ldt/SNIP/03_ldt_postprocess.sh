#!/usr/bin/env bash
#SBATCH --job-name=03_ldt_postprocess
#SBATCH --output=./log/03_ldt_postprocess_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --constraint="mil"
#SBATCH --mail-type=ALL
#SBATCH --account=s1189
#SBATCH --qos=debug

module purge
unset LD_LIBRARY_PATH
module use --append /home/emkemp/privatemodules/sles15
module load lisf_7.6_intel_2023.2.1_emk_aiml

set -euo pipefail

cd "$(dirname "$0")"

LDT_CONFIG="./03_ldt.config"

if [ ! -f "$LDT_CONFIG" ]; then
  echo "ERROR: LDT config file not found: $LDT_CONFIG"
  exit 1
fi

# Run LDT (single process here)
mpirun -np 1 ./LDT "$LDT_CONFIG"

echo "LDT post-processing complete"
