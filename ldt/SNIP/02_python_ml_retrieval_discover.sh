#!/bin/bash
#SBATCH --job-name=02_python_ml_retrieval
#SBATCH --output=logs/02/02_python_ml_retrieval_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=s1189
#SBATCH --qos=debug

# load modules
module purge
unset LD_LIBRARY_PATH
module use --append /home/emkemp/privatemodules/sles15
module load lisf_7.8_intel_2023.2.1_emk_aiml

cd $SLURM_SUBMIT_DIR

# TARGET_DATETIME must be in the format: YYYYMMDDHHMM (e.g. 202501200600)
if [ $# -lt 1 ]; then
  echo "Usage: $0 YYYYMMDDHHMM"
  echo "Example: $0 202501200600"
  exit 1
fi

TARGET_DATETIME="$1"

# read SNIP python config file
CONFIG="./02_python_ml_retrieval.json"
if [ ! -f "$CONFIG" ]; then
  echo "ERROR: SNIP config file not found: $CONFIG"
  exit 1
fi

# Run the python script
python ./SNIP_ops/main.py "$CONFIG" --input "WSF" --target-datetime "$TARGET_DATETIME"
