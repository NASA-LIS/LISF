#!/bin/sh
#SBATCH --job-name=wsf_sm_retrieval
#SBATCH --output=logs/wsf_sm_retrieval_%j.out
#SBATCH --error=logs/wsf_sm_retrieval_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --account=NWP601
#SBATCH --partition=batch
#SBATCH --cluster-constraint=blue
#SBATCH --exclusive
#SBATCH --mem=0
# =============================================================
# WSF Soil Moisture Retrieval, minimal SLURM wrapper
#
# Runs one job's worth of work — no self-submission, so it drops
# straight into a Cylc task. Cylc (or you) owns the cycling.
#
# Steps:
#   STEP=1  -> resampling only
#   STEP=2  -> prediction only
#   (unset) -> both steps
# Pass FORCE=1 to reprocess hours that already succeeded.
#
# Usage:
#   sbatch --export=START_DT=2024041300,END_DT=2024041500 wsf_retrieval_min.slurm
#   sbatch --export=START_DT=2024041300,END_DT=2024041500,STEP=1 wsf_retrieval_min.slurm
#   sbatch --export=START_DT=2024041300,END_DT=2024041500,STEP=2 wsf_retrieval_min.slurm
#   sbatch --export=START_DT=2024041300,END_DT=2024041500,FORCE=1 wsf_retrieval_min.slurm
# =============================================================
echo "=========================================="
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURM_NODELIST"
echo "CPUs:     $SLURM_CPUS_PER_TASK"
echo "Period:   ${START_DT} -> ${END_DT}"
echo "Step:     ${STEP:-both}"
echo "Force:    ${FORCE:-no}"
echo "Started:  $(date)"
echo "=========================================="

: "${START_DT:?START_DT not set (format YYYYMMDDHH)}"
: "${END_DT:?END_DT not set (format YYYYMMDDHH)}"

mkdir -p logs

# --- Environment ---
module use --append ~emkemp/hpc11/privatemodules
module load lisf_7.6_prgenv_cray_8.6.0_cpe_25.03_cce_19.0.0
module load afw-python/3.13-202605

ulimit -s unlimited
export JOBLIB_TEMP_FOLDER=${JOBLIB_TEMP_FOLDER:-/tmp}

./wsf_opl.py "$START_DT" "$END_DT" \
    --config wsf_opl.json \
    --batch-size "${MAX_PARALLEL:-$SLURM_CPUS_PER_TASK}" \
    ${STEP:+--step $STEP} \
    ${FORCE:+--force}

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
