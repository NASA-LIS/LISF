#!/bin/bash
#SBATCH --job-name=pred30m
#SBATCH --output=pred30m_%j.out
#SBATCH --error=pred30m_%j.err
#SBATCH --time=12:00:00
#SBATCH --time-min=8:00:00
#SBATCH --constraint=mil

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=480G
#SBATCH --account=s1189
#SBATCH --qos=allnccs
#SBATCH --no-requeue

# =============================================================
# WSF Stratified RF Prediction — SLURM Submission
# =============================================================
#
# Usage:
#   # Process a single month:
#   sbatch --export=MONTH=202404 submit_predict.slurm
#
#   # Process all months (one job):
#   sbatch --export=MODE=hourly submit_predict.slurm
#
#   # Submit all months as separate jobs (run from login node):
#   bash submit_predict.slurm --submit-all
#
#   # Zarr modes (unchanged):
#   sbatch --export=MODE=zarr submit_predict.slurm
#   sbatch --export=MODE=zarr,INPUT_NAME=202504_DES submit_predict.slurm
# =============================================================

# ---- Multi-job submission mode ----
# If called with --submit-all, submit one job per month and exit
if [[ "$1" == "--submit-all" ]]; then
    INPUT_BASE="/discover/nobackup/projects/usaf_lis/MET_FORCING/WSF/resampled"
    count=0
    for dir in "$INPUT_BASE"/[0-9][0-9][0-9][0-9][0-9][0-9]; do
        if [ -d "$dir" ]; then
            month=$(basename "$dir")
            echo "Submitting month: $month"
            sbatch --export=MONTH=$month --job-name="wsf_${month}" \
                   --output="wsf_${month}_%j.out" \
                   --error="wsf_${month}_%j.err" \
                   "$0"
            count=$((count + 1))
        fi
    done
    echo "Submitted $count jobs"
    exit 0
fi

# ---- Normal SLURM execution below ----
echo "=========================================="
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURM_NODELIST"
echo "CPUs:      $SLURM_CPUS_PER_TASK"
echo "MONTH:     ${MONTH:-all}"
echo "MODE:      ${MODE:-hourly}"
echo "Started:   $(date)"
echo "=========================================="

source ~/.bashrc
ulimit -s unlimited

# Error handler
error_exit() {
    echo "Error: $1" >&2
    exit 1
}

ml load anaconda || error_exit "Failed to load anaconda"
source activate dask || error_exit "Failed to activate env"

export PYTHONUNBUFFERED=1
export JOBLIB_TEMP_FOLDER=/discover/nobackup/ejalilva/tmp
mkdir -p $JOBLIB_TEMP_FOLDER

echo "Python: $(which python)"
echo "Working dir: $(pwd)"

# =============================================================
# Decide what to run
# =============================================================
MODE=${MODE:-hourly}

# If MONTH is set, override to hourly mode with that month
if [ -n "$MONTH" ]; then
    MODE="hourly_month"
fi

if [ "$MODE" = "zarr" ] && [ -n "$INPUT_NAME" ]; then
    MODE="zarr_single"
fi

case $MODE in

    # ----- Process single month (primary mode for parallel jobs) -----
    hourly_month)
        CMD="python -u run_hourly_stratified.py --month ${MONTH}"
        echo "Running: $CMD"
        $CMD
        ;;

    # ----- Process all hourly NC files -----
    hourly)
        CMD="python -u run_hourly_stratified.py"
        echo "Running: $CMD"
        $CMD
        ;;

    # ----- Process single zarr file -----
    zarr_single)
        CMD="python -u predict_sm_stratified.py zarr ${INPUT_NAME}"

        if [ "$HARD" = "1" ]; then
            CMD="$CMD --hard"
        else
            CMD="$CMD --top-k ${TOP_K:-3}"
        fi

        echo "Running: $CMD"
        $CMD
        ;;

    # ----- Process all zarr files -----
    zarr)
        CMD="python -u run_zarr_stratified.py"
        echo "Running: $CMD"
        $CMD
        ;;

    *)
        echo "ERROR: Unknown MODE=$MODE (use hourly, zarr)"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
