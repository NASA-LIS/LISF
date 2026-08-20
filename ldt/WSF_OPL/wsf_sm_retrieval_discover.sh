#!/bin/sh
#SBATCH --job-name=wsf_sm_retrieval
#SBATCH --output=logs/wsf_sm_retrieval_%j.out
#SBATCH --error=logs/wsf_sm_retrieval_%j.err
#SBATCH --time=00:29:00
#SBATCH --constraint=mil
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=480G
#SBATCH --account=s1189
#SBATCH --qos=debug
#SBATCH --no-requeue
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
#source ~/.bashrc
#module purge
#module use --append $HOME/.privatemodules
#module load lisf_7.5_intel_2023.2.1
#module load miniforge
#conda activate wsf_sm

# On Discover
module purge
module use --append /home/emkemp/privatemodules/sles15
module load lisf_7.8_intel_2023.2.1_emk_aiml

ulimit -s unlimited
export JOBLIB_TEMP_FOLDER=${JOBLIB_TEMP_FOLDER:-/tmp}

./WSF_OPL/prediction/wsf_opl.py "$START_DT" "$END_DT" \
    --config wsf_opl.json \
    --batch-size "${MAX_PARALLEL:-$SLURM_CPUS_PER_TASK}" \
    ${STEP:+--step $STEP} \
    ${FORCE:+--force}

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
