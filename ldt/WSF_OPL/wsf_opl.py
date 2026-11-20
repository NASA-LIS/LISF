#!/usr/bin/env python3
"""
WSF Soil Moisture Retrieval — Steps 1 & 2

  Step 1: OPL resampling  — batched parallel LDT
  Step 2: ML prediction   — sequential hours, all cores per hour

Usage:
  python wsf_opl.py START_DT END_DT [--config FILE] [--batch-size N] [--force] [--step {1,2}]

  START_DT / END_DT : period to process, format YYYYMMDDHH (e.g. 2024041300)
  --config          : path to pipeline_config.json (default: ./pipeline_config.json)
  --batch-size      : max parallel LDT jobs for step 1 (default: cpu count)
  --force           : ignore checkpoint and reprocess all hours
  --step            : run only step 1 (resampling) or step 2 (prediction); omit to run both

Examples:
  # Run both steps
  python wsf_opl.py 2024041300 2026041300

  # Run step 1 only (resampling)
  python wsf_opl.py 2024041300 2026041300 --step 1

  # Run step 2 only (prediction, resampled files already exist)
  python wsf_opl.py 2024050100 2025060100 --step 2

  # Reprocess all hours ignoring checkpoint
  python wsf_opl.py 2025090100 2025100100 --force

  # Use a config file at a non-default location
  python wsf_opl.py 2024041300 2026041300 --config /path/to/my_config.json

SLURM:
  sbatch --export=START_DT=2024041300,END_DT=2026041300 wsf_retrieval.slurm
  sbatch --export=START_DT=2024041300,END_DT=2026041300,STEP=1 wsf_retrieval.slurm
  sbatch --export=START_DT=2024041300,END_DT=2026041300,STEP=2 wsf_retrieval.slurm
"""

# --- standard library ---
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

V522_LAST  = "2025041323"   # v5.2.2 directory layout ends here
CHECKPOINT = Path("logs/resampling_checkpoint.txt")
DT_RE      = re.compile(r"WSFM_01_d(\d{8})_t(\d{2})")


def find_hours(base, start, end, lo="0"*10, hi="9"*10):
    """Scan a directory tree for WSF files; return YYYYMMDDHH strings in range."""
    if not Path(base).is_dir():
        return set()
    hours = set()
    for f in Path(base).rglob("*WSFM_01_d*_res_sdr.nc"):
        m = DT_RE.search(f.name)
        if m:
            dt = m.group(1) + m.group(2)
            if start <= dt <= end and lo <= dt <= hi:
                hours.add(dt)
    return hours


def src_dir(dt, v522_base, flat_base):
    """Return the raw SDR directory for a given YYYYMMDDHH."""
    if dt <= V522_LAST:
        return Path(v522_base) / dt[:4] / dt[4:6] / dt[6:8]
    return Path(flat_base) / dt[:6]


def write_config(template, dt, src, out):
    """Copy the LDT template and substitute the per-hour fields."""
    text = template.read_text(encoding='utf-8')
    for pattern, value in [
        (r"WSF valid date \(YYYYMMDDHH\):.*", f"WSF valid date (YYYYMMDDHH):            {dt}"),
        (r"WSF input directory:.*",            f"WSF input directory:                    {src}"),
        (r"WSF output directory:.*",           f"WSF output directory:                   {out}"),
        (r"LDT diagnostic file:.*",            f"LDT diagnostic file:                    logs/ldtlog_{dt}"),
        (r"WSF filelist suffix number:.*",     f"WSF filelist suffix number:             {dt}"),
    ]:
        text = re.sub(pattern, value, text)
    cfg = Path(f"ldt.config.wsf.{dt}")
    cfg.write_text(text, encoding='utf-8')
    return cfg


# =============================================================================
# Step 1 — OPL Resampling
# =============================================================================
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
def step1_resample(s1, program, template, out_base, start, end, batch_size, force):
    """Run LDT OPL resampling in parallel batches over all pending hours."""
    print("\n=== STEP 1: OPL Resampling ===")

    Path("logs").mkdir(exist_ok=True)
    CHECKPOINT.touch()

    hours  = find_hours(s1["v522_sdr_base"], start, end, hi=V522_LAST)
    hours |= find_hours(s1["raw_sdr_base"],  start, end, lo=str(int(V522_LAST) + 1))
    done    = set() if force else set(CHECKPOINT.read_text(encoding='utf-8').split())
    pending = sorted(hours - done)

    print(f"  Hours found: {len(hours)}  done: {len(done & hours)}  pending: {len(pending)}")

    if not pending:
        print("  All hours already processed.")
        return

    n_ok = n_fail = n_skip = 0
    n_batches = (len(pending) + batch_size - 1) // batch_size

    for b, i in enumerate(range(0, len(pending), batch_size)):
        batch = pending[i : i + batch_size]
        print(f"\n  Batch {b+1}/{n_batches} ({len(batch)} hours)")

        procs = []
        for dt in batch:
            sd = src_dir(dt, s1["v522_sdr_base"], s1["raw_sdr_base"])
            if not sd.is_dir():
                print(f"    SKIP {dt}: {sd} not found")
                n_skip += 1
                continue
            od = out_base / dt[:6]
            od.mkdir(parents=True, exist_ok=True)
            cfg_file = write_config(template, dt, sd, od)
            # pylint: disable=consider-using-with
            # Intentional: processes are launched as a batch and waited on
            # sequentially below — a context manager cannot span both loops.
            proc = subprocess.Popen(
                [str(program), str(cfg_file)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            # pylint: enable=consider-using-with
            procs.append((dt, proc, cfg_file))

        for dt, proc, cfg_file in procs:
            rc  = proc.wait()
            log = Path(f"logs/ldtlog_{dt}")
            if rc != 0:
                print(f"    FAIL {dt}: exit code {rc}")
                n_fail += 1
            elif log.exists() and "No WSF files found" in log.read_text(encoding='utf-8'):
                print(f"    WARN {dt}: no input files found")
                n_fail += 1
            else:
                with CHECKPOINT.open("a", encoding='utf-8') as fh:
                    fh.write(dt + "\n")
                n_ok += 1
            cfg_file.unlink(missing_ok=True)
            Path(f"WSF_filelist_{dt}.dat").unlink(missing_ok=True)

        print(f"    ok={n_ok} fail={n_fail} skip={n_skip}")

    n_out = sum(1 for _ in out_base.rglob("WSF_SDR_resampled_*.nc"))
    print(f"\n  Resampled files: {n_out}  ok: {n_ok}  fail: {n_fail}  skip: {n_skip}")
    if n_out == 0:
        sys.exit("ERROR: no resampled files produced — check logs/")
# pylint: enable=too-many-arguments,too-many-positional-arguments,too-many-locals


# =============================================================================
# Step 2 — ML Soil Moisture Prediction
# =============================================================================
# pylint: disable=too-many-locals
def step2_predict(s2, out_base, start, end, force=False):
    """Load the stratified RF predictor and run inference on all resampled files."""
    print("\n=== STEP 2: ML Prediction ===")

    script_path = s2["predict_script"]
    sys.path.insert(0, str(Path(script_path).parent))
    # pylint: disable=import-outside-toplevel
    # Dynamic import: script location comes from pipeline_config.json at runtime.
    from predict_sm_stratified import StratifiedSoilMoisturePredictor
    # pylint: enable=import-outside-toplevel

    predictor = StratifiedSoilMoisturePredictor(
        model_dir=s2["model_dir"],
        cluster_dir=s2["cluster_dir"],
        env_paths=s2["env_paths"],
        top_k=int(s2["top_k"]),
    )

    files = sorted(
        f for f in out_base.rglob("WSF_SDR_resampled_*.nc")
        if (m := re.search(r"WSF_SDR_resampled_(\d{8})_t(\d{4})_(ASC|DES)", f.name))
        and start <= m.group(1) + m.group(2)[:2] <= end
    )
    print(f"  Files to predict: {len(files)}")
    if not files:
        sys.exit("ERROR: no resampled files found for prediction")

    sm_out = Path(s2["sm_output_dir"])
    sm_out.mkdir(parents=True, exist_ok=True)

    n_ok = n_fail = n_skip = 0
    t0_all = time.time()

    for i, f in enumerate(files):
        # Skip if the output SM file already exists (unless --force was given).
        m = re.search(r"WSF_SDR_resampled_(\d{8})_t(\d{4})_(ASC|DES)", f.name)
        date, hhmm, orbit = m.group(1), m.group(2), m.group(3)
        if not force and any(sm_out.glob(f"ARFS_SM_WSFM_*{date}*{hhmm}*{orbit}*.nc")):
            n_skip += 1
            continue
        
        t0 = time.time()

        try:
            predictor.predict_from_hourly_nc(f, output_dir=sm_out)
            n_ok += 1
        except Exception as e:  # pylint: disable=broad-exception-caught
            # Catch-all is intentional: one bad file must not abort the full run.
            print(f"  ERROR {f.name}: {e}")
            n_fail += 1

        if (i + 1) % 10 == 0 or (i + 1) == len(files):
            print(f"  [{i+1}/{len(files)}] ok={n_ok} fail={n_fail} "
                  f"last={time.time()-t0:.1f}s total={time.time()-t0_all:.0f}s")

    print(f"\n  Done: {n_ok} ok, {n_fail} failed")
    if n_fail:
        sys.exit(1)
# pylint: enable=too-many-locals


# =============================================================================
# Entry point
# =============================================================================
def main():
    """Parse arguments, load config, and dispatch to step 1 and/or step 2."""
    parser = argparse.ArgumentParser()
    parser.add_argument("start_dt", nargs="?", default=os.environ.get("START_DT"))
    parser.add_argument("end_dt",   nargs="?", default=os.environ.get("END_DT"))
    parser.add_argument("--config",     default="./pipeline_config.json")
    parser.add_argument("--batch-size", type=int,
                        default=int(os.environ.get(
                            "MAX_PARALLEL",
                            max(10, (os.cpu_count() or 10) // 10 * 10)
                        )))
    parser.add_argument("--force", action="store_true",
                        help="Ignore checkpoint, reprocess all hours")
    parser.add_argument("--step", type=int, choices=[1, 2], default=None,
                        help="Run only step 1 or 2 (default: both)")
    args = parser.parse_args()

    if not args.start_dt or not args.end_dt:
        parser.error("START_DT and END_DT required")

    cfg      = json.loads(Path(args.config).read_text(encoding='utf-8'))
    s1       = cfg["step1_wsf_opl"]
    program  = Path(cfg["executables"]["ldt"]).resolve()
    template = Path(s1["ldt_config_template"]).resolve()
    out_base = Path(s1["resampled_base"])

    print(f"Period: {args.start_dt} → {args.end_dt}")

    if args.step != 2:
        step1_resample(s1, program, template, out_base,
                       args.start_dt, args.end_dt, args.batch_size, args.force)

    if args.step != 1:
        step2_predict(cfg["step2_prediction"], out_base, args.start_dt, args.end_dt, args.force)

    n_sm = sum(
        1 for _ in
        Path(cfg["step2_prediction"]["sm_output_dir"]).rglob("ARFS_SM_WSFM_*.nc")
    )
    print(f"\n=== DONE — SM files: {n_sm} ===")


if __name__ == "__main__":
    main()