#!/bin/bash
#SBATCH --job-name=SNIP_{input_type}_{dt}
#SBATCH --output=./log/SNIP_{input_type}_{dt}.log
{sbatch_header}

ulimit -s unlimited

# ── Setup ──────────────────────────────────────────────────────────────
INPUT_TYPE="{input_type}"           # Will be "AMSR2" or "WSF"
INPUT_LOWER="${INPUT_TYPE,,}"       # Automatically converts to lowercase: "amsr2" or "wsf"

TARGET_DATETIME="{dt}"
TARGET_DATETIME_10="{dt10}"

PMW_DIR="./data/input/PMW_SD"
SNIP_OUT="./SNIP_ops/data/output/${INPUT_LOWER}_snip_0p1deg_${TARGET_DATETIME}_AFgrid.nc"
TARGET_NC="$PMW_DIR/${INPUT_LOWER}_snip_0p1deg_${TARGET_DATETIME}_AFgrid.nc"

LDT_OUT="./data/output/snip/SNIP_${TARGET_DATETIME_10}.nc"
LDTLOG="./data/output/logs/ldtlog.0000"

mkdir -p ./log "$PMW_DIR"
[ -n "$SLURM_SUBMIT_DIR" ] && { cd "$SLURM_SUBMIT_DIR" || exit 1; }

log() { echo "$*"; }
die() { log "ERROR: $*"; exit 1; }

log "=== SNIP + LDT started $(date) | system={system} | input=${INPUT_TYPE} | dt=${TARGET_DATETIME} ==="

# ── Load modules ───────────────────────────────────────────────────────
{modules}

# ── STEP 0: WSF Resampling (Only runs if input is WSF) ─────────────────
if [ "$INPUT_TYPE" = "WSF" ]; then
    YYYY=${TARGET_DATETIME:0:4}
    MM=${TARGET_DATETIME:4:2}
    DD=${TARGET_DATETIME:6:2}
    HH=${TARGET_DATETIME:8:2}

    # Calculate the start time (6 hours ago)
    START_DATETIME=$(date -u -d "${YYYY}-${MM}-${DD} ${HH}:00:00Z 6 hours ago" +"%Y%m%d%H%M")
    
    log "--- STEP 0: Running WSF resampling from ${START_DATETIME} to ${TARGET_DATETIME} ---"
    python submit_job.py --input WSF --resample "${START_DATETIME}" "${TARGET_DATETIME}" || die "WSF resampling failed"
    log "STEP 0: WSF resampling finished ---"
fi

# ── STEP 1: Python SNIP ────────────────────────────────────────────────
run_snip() {
    local tmp_cfg
    tmp_cfg=$(mktemp ./config/SNIP_config_${TARGET_DATETIME}_XXXXXX.json) || return 1
    cp "./config/SNIP_config.json" "$tmp_cfg"
    
    # Update target_datetime in the JSON (your original code)
    sed -i '/\"target_datetime\"/s/: \"[^\"]*\"/: \"'${TARGET_DATETIME}'\"/g' "$tmp_cfg"
    
    # Run main.py with the config file AND the new --input flag
    python main.py "$tmp_cfg" --input "${INPUT_TYPE}" 2>&1
    
    local rc=$?  # Safely capture the exit code of python
    rm -f "$tmp_cfg"
    return $rc
}

snip_done=false
if [ -e "$TARGET_NC" ] || [ -e "$SNIP_OUT" ]; then
    log "STEP 1: Passive Microwave derived snow depth already exists – skipping"
    snip_done=true
else
    log "STEP 1: Running SNIP Python model..."
    if cd ./SNIP_ops; then
        run_snip && snip_done=true || log "WARNING: SNIP failed – continuing to LDT"
        cd ..
    else
        log "WARNING: Could not cd into ./SNIP_ops – skipping SNIP"
    fi
fi

# Move SNIP output into the PMW input dir for LDT
if [ -e "$SNIP_OUT" ] && [ ! -e "$TARGET_NC" ]; then
    mv "$SNIP_OUT" "$PMW_DIR"
fi

# ── STEP 2: LDT ────────────────────────────────────────────────────────
log "STEP 2: Running SNIP LDT..."
[ ! -f ./ldt.config ] && die "ldt.config not found"
[ ! -x ./LDT ]   && die "SNIP LDT executable not found"

if [ -f "$LDT_OUT" ]; then
    log "LDT output already exists – skipping"
else
    # 1. Replace the Date
    sed -i "s/SNIP valid date (YYYYMMDDHH):[[:space:]]*[0-9]\{10\}/SNIP valid date (YYYYMMDDHH):                    ${TARGET_DATETIME_10}/g" ldt.config

    # 2. Replace the Input Type (amsr2 or wsf)
    sed -i "s/\(SNIP PMW snow depth filename prefix:[[:space:]]*\)[a-zA-Z0-9_]*/\1${INPUT_LOWER}/g" ldt.config
    
    
    mpirun -np 1 ./LDT ldt.config 2>&1
    [ -f "$LDT_OUT" ] || die "SNIP LDT finished but output file not found"

    # Rename ldtlog.0000 to preserve it across runs
    if [ -f "$LDTLOG" ]; then
        mv "$LDTLOG" "${LDTLOG}_${TARGET_DATETIME_10}"
        log "Renamed ldtlog.0000 to ldtlog.0000_${TARGET_DATETIME_10}"
    fi
fi

# ── Summary ────────────────────────────────────────────────────────────
log "=== DONE $(date) ==="
$snip_done && log "SNIP: SUCCESS" || log "SNIP: FAILED"
log "LDT: SUCCESS"