#!/bin/bash
#SBATCH --job-name=SNIP_{dt}
#SBATCH --output=./log/SNIP_{dt}.log
{sbatch_header}

ulimit -s unlimited

# ── Setup ──────────────────────────────────────────────────────────────
TARGET_DATETIME="{dt}"
TARGET_DATETIME_10="{dt10}"
SNIP_OUT="./SNIP_ops/data/output/wsf_snip_0p1deg_{dt}_AFgrid.nc"
WSF_DIR="./data/input/WSF_SD"
LDT_OUT="./data/output/snip/SNIP_{dt10}.nc"
LDTLOG="./data/output/logs/ldtlog.0000"

mkdir -p ./log "$WSF_DIR"
[ -n "$SLURM_SUBMIT_DIR" ] && { cd "$SLURM_SUBMIT_DIR" || exit 1; }

log() { echo "$*"; }
die() { log "ERROR: $*"; exit 1; }

log "=== SNIP + LDT started $(date) | system={system} | dt={dt} ==="

# ── Load modules ───────────────────────────────────────────────────────
{modules}

# ── STEP 0: LDT WSF resampling ────────────────────────────────────────────────
# 1. Extract year, month, day, and hour from TARGET_DATETIME
YYYY=${TARGET_DATETIME:0:4}
MM=${TARGET_DATETIME:4:2}
DD=${TARGET_DATETIME:6:2}
HH=${TARGET_DATETIME:8:2}

# 2. Calculate the start time (6 hours ago) using the Linux date command
START_DATETIME=$(date -u -d "${YYYY}-${MM}-${DD} ${HH}:00:00Z 6 hours ago" +"%Y%m%d%H%M")

log "--- STEP 0: Running WSF resampling from ${START_DATETIME} to ${TARGET_DATETIME} ---"

# 3. Call the Python script in "Worker Mode" using the --resample flag
python submit_job_WSF.py --resample "${START_DATETIME}" "${TARGET_DATETIME}" || die "WSF resampling failed"


# ── STEP 1: Python SNIP ────────────────────────────────────────────────
run_snip() {
    local cfg_src="./config/SNIP_config.json"
    local tmp_cfg
    tmp_cfg=$(mktemp ./config/SNIP_config_{dt}_XXXXXX.json) || {
        log "WARNING: mktemp failed"
        return 1
    }
    cp "$cfg_src" "$tmp_cfg"
    sed -i '/\"target_datetime\"/s/: \"[^\"]*\"/: \"{dt}\"/g' "$tmp_cfg"
    python main.py "$tmp_cfg" 2>&1
    local rc=$PIPESTATUS
    rm -f "$tmp_cfg"
    return $rc
}

snip_done=false
if [ -e "$WSF_DIR/wsf_snip_0p1deg_{dt}_AFgrid.nc" ] || \
   [ -e "$SNIP_OUT" ]; then
    log "STEP 1: SNIP output already exists – skipping"
    snip_done=true
else
    log "STEP 1: Running Python SNIP..."
    if cd ./SNIP_ops; then
        run_snip && snip_done=true || \
            log "WARNING: SNIP failed – continuing to LDT"
        cd ..
    else
        log "WARNING: Could not cd into ./SNIP_ops – skipping SNIP"
    fi
fi

# Move SNIP output into the WSF input dir for LDT
if [ -e "$SNIP_OUT" ] && \
   [ ! -e "$WSF_DIR/$(basename "$SNIP_OUT")" ]; then
    mv "$SNIP_OUT" "$WSF_DIR"
fi

# ── STEP 2: LDT ────────────────────────────────────────────────────────
log "STEP 2: Running LDT..."
[ ! -f ./ldt.config ] && die "ldt.config not found"
[ ! -x ./LDT_SNIP ]        && die "LDT executable not found"

if [ -f "$LDT_OUT" ]; then
    log "LDT output already exists – skipping"
else
    sed -i "s/SNIP valid date (YYYYMMDDHH):[[:space:]]*[0-9]\{10\}/SNIP valid date (YYYYMMDDHH):                    {dt10}/g" \
        ldt.config
    mpirun -np 1 ./LDT_SNIP ldt.config 2>&1
    [ -f "$LDT_OUT" ] || die "LDT finished but output file not found"

     # ── Rename ldtlog.0000 to ldtlog.0000_{dt10} ──────────────────────
    
    if [ -f "$LDTLOG" ]; then
        mv "$LDTLOG" "${LDTLOG}_{dt10}"
        log "Renamed ldtlog.0000 to ldtlog.0000_{dt10}"
    else
        log "WARNING: ldtlog.0000 not found – skipping rename"
    fi
    
fi

# ── Summary ────────────────────────────────────────────────────────────
log "=== DONE $(date) ==="
$snip_done && log "SNIP: SUCCESS" || log "SNIP: FAILED"
log "LDT:  SUCCESS"