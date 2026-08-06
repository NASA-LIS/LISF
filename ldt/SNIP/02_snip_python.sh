#!/usr/bin/env bash
#SBATCH --job-name=02_snip_python
#SBATCH --output=./log/02_snip_python_%j.log
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
module load lisf_7.6_intel_2023.2.1_emk_aiml

cd "$(dirname "$0")"

# TARGET_DATETIME must be in the format: YYYYMMDDHHMM (e.g. 202501200600)
if [ $# -lt 1 ]; then
  echo "Usage: $0 YYYYMMDDHHMM [WSF|AMSR2]"
  echo "Example: $0 202501200600 WSF"
  exit 1
fi

TARGET_DATETIME="$1"
INPUT_TYPE="${2:-WSF}"

# read SNIP python config file
CONFIG="./SNIP_ops/config/SNIP_config.json"
if [ ! -f "$CONFIG" ]; then
  echo "ERROR: SNIP config file not found: $CONFIG"
  exit 1
fi

python ./SNIP_ops/main.py "$CONFIG" --input "$INPUT_TYPE" --target-datetime "$TARGET_DATETIME"
