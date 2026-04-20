#!/usr/bin/env python3
"""
submit_job.py  –  Unified SNIP Pipeline Manager for AMSR2 & WSF

Mode 1 (Submitter): Submit SNIP + LDT SLURM jobs for Discover or HPC11.
    Usage: python submit_job.py <YYYYMMDDHHMM> --input [AMSR2|WSF] [--system discover|hpc11] [--dry-run]

    Examples:
        python submit_job.py 202501200600 --input WSF
        python submit_job.py 202501200600 --input AMSR2 --system hpc11

Mode 2 (Worker): Run OPL Resampling (WSF ONLY - called automatically inside the SLURM job).
    Usage: python submit_job.py --input WSF --resample <START_DT> <END_DT>
"""

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# =============================================================================
# SUBMISSION CONSTANTS & PROFILES
# =============================================================================
VALID_HOURS = {"00", "06", "12", "18"}

# We now use a single, unified template for both AMSR2 and WSF
TEMPLATE_PATH = Path("./job_template.sh")

SYSTEMS = {
    "discover": {
        "sbatch_header": """\
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=1:00:00
#SBATCH --constraint="mil"
#SBATCH --mail-type=ALL
#SBATCH --account=s1189
#SBATCH --qos=debug""",
        "modules": """\
module purge
unset LD_LIBRARY_PATH
module use --append /home/emkemp/privatemodules/sles15
module load lisf_7.6_intel_2023.2.1_emk_aiml""",
    },
    "hpc11": {
        "sbatch_header": """\
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1:00:00
#SBATCH --cluster-constraint=blue
#SBATCH --account=NWP601
#SBATCH --exclusive
#SBATCH --mem=0""",
        "modules": """\
module use --append /ccs/home/emkemp/hpc11/privatemodules/
module load lisf_7.6_prgenv_cray_8.6.0_cpe_25.03_cce_19.0.0
module use --append /sw/afw_sw/modulefiles
module load afw-python/3.11-202511""",
    },
}

# =============================================================================
# RESAMPLING CONSTANTS (WSF ONLY)
# =============================================================================
V522_LAST  = "2025041323"
CHECKPOINT = Path("log/resampling_checkpoint.txt")
DT_RE      = re.compile(r"WSFM_01_d(\d{8})_t(\d{2})")


# =============================================================================
# MODE 1: SUBMISSION LOGIC
# =============================================================================
def validate_datetime(dt_str: str) -> str:
    """Return dt_str if valid, otherwise raise ValueError."""
    if not re.fullmatch(r"\d{12}", dt_str):
        raise ValueError(f"Datetime must be exactly 12 digits (YYYYMMDDHHMM), got: '{dt_str}'")
    hour = dt_str[8:10]
    if hour not in VALID_HOURS:
        raise ValueError(f"Hour must be one of {sorted(VALID_HOURS)}, got: '{hour}'")
    try:
        datetime.strptime(dt_str[:10], "%Y%m%d%H")
    except ValueError as exc:
        raise ValueError(f"Invalid calendar date in '{dt_str}': {exc}") from exc
    return dt_str

def detect_system() -> str:
    """Best-effort auto-detect of the current HPC system."""
    hostname = os.environ.get("HOSTNAME", "")
    if "hpc11" in hostname or "frontier" in hostname:
        return "hpc11"
    return "discover"

def build_script(dt: str, system: str, input_type: str) -> str:
    """Read the unified job_template and substitute placeholders."""
    if not TEMPLATE_PATH.exists():
        raise FileNotFoundError(f"ERROR: Template file not found: {TEMPLATE_PATH}")

    cfg = SYSTEMS[system]
    dt10 = dt[:10]  # YYYYMMDDHH – no minutes
    template = TEMPLATE_PATH.read_text()

    script = (template
        .replace("{dt}", dt)
        .replace("{dt10}", dt10)
        .replace("{input_type}", input_type)
        .replace("{system}", system)
        .replace("{sbatch_header}", cfg["sbatch_header"])
        .replace("{modules}", cfg["modules"])
    )
    return script

def submit(dt: str, system: str, input_type: str, dry_run: bool = False) -> None:
    log_dir = Path("./log")
    log_dir.mkdir(parents=True, exist_ok=True)

    script_content = build_script(dt, system, input_type)
    script_path = log_dir / f"job_{input_type}_{dt}.sh"
    script_path.write_text(script_content)
    script_path.chmod(0o755)

    print(f"[{datetime.now()}] System={system}  Input={input_type}  Datetime={dt}")
    print(f"  Script : {script_path}")
    print(f"  Log    : {log_dir / f'SNIP_LDT_{input_type}_{dt}.log'}")

    if dry_run:
        print("  DRY RUN – sbatch not called. Script written but not submitted.")
        print("\n--- Script preview (first 25 lines) ---")
        print("\n".join(script_content.splitlines()[:25]))
        print("...\n")
        return

    result = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True)
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        print(f"  Submitted job ID: {job_id}")
    else:
        print(f"  ERROR submitting job: {result.stderr.strip()}")
        sys.exit(1)

    print("\nCurrent queue:")
    subprocess.run(["squeue", "-u", os.environ.get("USER", os.environ.get("LOGNAME", "unknown"))])


# =============================================================================
# MODE 2: RESAMPLING LOGIC (WSF WORKER ONLY)
# =============================================================================
def find_hours(base, start, end, lo="0"*10, hi="9"*10):
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
    if dt <= V522_LAST:
        return Path(v522_base) / dt[:4] / dt[4:6] / dt[6:8]
    return Path(flat_base) / dt[:6]

def write_config(template, dt, src, out):
    text = template.read_text(encoding='utf-8')
    for pattern, value in [
        (r"WSF valid date \(YYYYMMDDHH\):.*", f"WSF valid date (YYYYMMDDHH):            {dt}"),
        (r"WSF input directory:.*",            f"WSF input directory:                    {src}"),
        (r"WSF output directory:.*",           f"WSF output directory:                   {out}"),
        (r"LDT diagnostic file:.*",            f"LDT diagnostic file:                    log/ldtlog_{dt}"),
        (r"WSF filelist suffix number:.*",     f"WSF filelist suffix number:             {dt}"),
    ]:
        text = re.sub(pattern, value, text)
    cfg = Path(f"ldt.config.wsf.{dt}")
    cfg.write_text(text, encoding='utf-8')
    return cfg

def run_resampling(cfg, program, template, out_base, start, end, batch_size, force):
    print("\n=== OPL Resampling (WSF) ===")
    Path("log").mkdir(exist_ok=True)
    CHECKPOINT.touch()

    hours  = find_hours(cfg["v522_sdr_base"], start, end, hi=V522_LAST)
    hours |= find_hours(cfg["raw_sdr_base"],  start, end, lo=str(int(V522_LAST) + 1))
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
            sd = src_dir(dt, cfg["v522_sdr_base"], cfg["raw_sdr_base"])
            if not sd.is_dir():
                print(f"    SKIP {dt}: {sd} not found")
                n_skip += 1
                continue
            od = out_base / dt[:6]
            od.mkdir(parents=True, exist_ok=True)
            cfg_file = write_config(template, dt, sd, od)

            proc = subprocess.Popen([str(program), str(cfg_file)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            procs.append((dt, proc, cfg_file))

        for dt, proc, cfg_file in procs:
            rc  = proc.wait()
            log = Path(f"log/ldtlog_{dt}")
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
        sys.exit("ERROR: no resampled files produced — check log/")


# =============================================================================
# ENTRY POINT (DISPATCHER)
# =============================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required Input selection
    parser.add_argument("--input", choices=["AMSR2", "WSF"], required=True, help="Input data type to process (AMSR2 or WSF)")

    # Mode 1: Submission Arguments
    parser.add_argument("datetime", nargs="?", default=None, help="Target datetime in YYYYMMDDHHMM format")
    parser.add_argument("--system", choices=list(SYSTEMS), default=None, help="HPC system profile (default: auto-detect)")
    parser.add_argument("--dry-run", action="store_true", help="Write the script but do not call sbatch")

    # Mode 2: Worker Arguments (WSF Only)
    parser.add_argument("--resample", nargs=2, metavar=("START_DT", "END_DT"), help="Worker mode: Run OPL resampling (WSF only)")
    parser.add_argument("--config", default="./SNIP_ops/config/SNIP_config.json", help="Path to config (for WSF resampling)")
    parser.add_argument("--batch-size", type=int, default=int(os.environ.get("MAX_PARALLEL", max(10, (os.cpu_count() or 10) // 10 * 10))), help="Max parallel LDT jobs")
    parser.add_argument("--force", action="store_true", help="Ignore checkpoint during resampling")

    args = parser.parse_args()

    # ---------------------------------------------------------
    # Dispatch to Mode 2: Worker (Resampling)
    # ---------------------------------------------------------
    if args.resample:
        if args.input != "WSF":
            parser.error("Resampling is only supported when --input WSF is selected.")

        config_path = Path(args.config)
        if not config_path.exists():
            sys.exit(f"ERROR: Configuration file '{args.config}' not found.")

        start_dt, end_dt = args.resample
        cfg = json.loads(config_path.read_text(encoding='utf-8'))
        program = Path(cfg["ldt"]).resolve()
        template = Path(cfg["ldt_config_template"]).resolve()
        out_base = Path(cfg["resampled_base"])

        run_resampling(cfg, program, template, out_base, start_dt, end_dt, args.batch_size, args.force)
        sys.exit(0)

    # ---------------------------------------------------------
    # Dispatch to Mode 1: Submitter
    # ---------------------------------------------------------
    if not args.datetime:
        parser.error("You must provide a target datetime for job submission (e.g., 202601010600).")

    try:
        dt = validate_datetime(args.datetime)
    except ValueError as exc:
        parser.error(str(exc))

    system = args.system or detect_system()
    print(f"Using system profile: '{system}'")

    submit(dt, system, args.input, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
