#!/bin/sh
#SBATCH --job-name=01_ldt_wsf_resample
#SBATCH --output=./logs/01/01_ldt_wsf_resample_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:30:00
#SBATCH --account=NWP601
#SBATCH --partition=batch
#SBATCH --cluster-constraint=blue
#SBATCH --exclusive
#SBATCH --mem=0


cd $SLURM_SUBMIT_DIR

ulimit -s unlimited

CFG="./01_ldt_wsf_resample.json"

if [ ! -f "$CFG" ]; then
  echo "ERROR: WSF config file not found: $CFG"
  exit 1
fi

# User must provide the target datetime as the first argument
if [ $# -lt 1 ]; then
  echo "Usage: $0 YYYYMMDDHHMM"
  echo "Example: $0 202501200600"
  exit 1
fi

TARGET_DATETIME="$1"


./SNIP_ops/data_processing/WSF_reader.py --resample "$TARGET_DATETIME" --config "$CFG"
