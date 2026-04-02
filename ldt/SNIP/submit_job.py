#!/usr/bin/env python3
"""
submit_job.py  –  Submit SNIP + LDT SLURM jobs for Discover or HPC11.

Usage:
    python submit_job.py <YYYYMMDDHHMM> [--system discover|hpc11] [--dry-run]

Examples:
    python submit_job.py 202501201200
    python submit_job.py 202501200600 --system hpc11
    python submit_job.py 202501200600 --dry-run
"""

import argparse
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# System profiles
# ──────────────────────────────────────────────────────────────────────────────
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

# ──────────────────────────────────────────────────────────────────────────────
# Validation
# ──────────────────────────────────────────────────────────────────────────────
VALID_HOURS = {"00", "06", "12", "18"}
TEMPLATE_PATH = Path("./job_template.sh")


def validate_datetime(dt_str: str) -> str:
    """Return dt_str if valid, otherwise raise ValueError."""
    if not re.fullmatch(r"\d{12}", dt_str):
        raise ValueError(
            f"Datetime must be exactly 12 digits (YYYYMMDDHHMM), got: '{dt_str}'"
        )
    hour = dt_str[8:10]
    if hour not in VALID_HOURS:
        raise ValueError(
            f"Hour must be one of {sorted(VALID_HOURS)}, got: '{hour}'"
        )
    try:
        datetime.strptime(dt_str[:10], "%Y%m%d%H")
    except ValueError as exc:
        raise ValueError(
            f"Invalid calendar date in '{dt_str}': {exc}"
        ) from exc
    return dt_str


def detect_system() -> str:
    """Best-effort auto-detect of the current HPC system."""
    hostname = os.environ.get("HOSTNAME", "")
    if "hpc11" in hostname or "frontier" in hostname:
        return "hpc11"
    return "discover"  # default


# ──────────────────────────────────────────────────────────────────────────────
# Script generation
# ──────────────────────────────────────────────────────────────────────────────
def build_script(dt: str, system: str) -> str:
    """Read job_template.sh and substitute placeholders."""
    if not TEMPLATE_PATH.exists():
        raise FileNotFoundError(f"Template not found: {TEMPLATE_PATH}")

    cfg = SYSTEMS[system]
    dt10 = dt[:10]  # YYYYMMDDHH – no minutes

    template = TEMPLATE_PATH.read_text()

    # Replace all placeholders
    script = (template
        .replace("{dt}", dt)
        .replace("{dt10}", dt10)
        .replace("{system}", system)
        .replace("{sbatch_header}", cfg["sbatch_header"])
        .replace("{modules}", cfg["modules"])
    )

    return script


# ──────────────────────────────────────────────────────────────────────────────
# Submission
# ──────────────────────────────────────────────────────────────────────────────
def submit(dt: str, system: str, dry_run: bool = False) -> None:
    log_dir = Path("./log")
    log_dir.mkdir(parents=True, exist_ok=True)

    script_content = build_script(dt, system)
    script_path = log_dir / f"job_{dt}.sh"
    script_path.write_text(script_content)
    script_path.chmod(0o755)

    print(f"[{datetime.now()}] System={system}  Datetime={dt}")
    print(f"  Script : {script_path}")
    print(f"  Log    : {log_dir / f'SNIP_LDT_{dt}.log'}")

    if dry_run:
        print("  DRY RUN – sbatch not called. Script written but not submitted.")
        print("\n--- Script preview (first 30 lines) ---")
        print("\n".join(script_content.splitlines()[:30]))
        print("...")
        return

    result = subprocess.run(
        ["sbatch", str(script_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        print(f"  Submitted job ID: {job_id}")
    else:
        print(f"  ERROR submitting job: {result.stderr.strip()}")
        sys.exit(1)

    # Show queue status
    print("\nCurrent queue:")
    subprocess.run([
        "squeue", "-u",
        os.environ.get("USER", os.environ.get("LOGNAME", "unknown"))
    ])


# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "datetime",
        help="Target datetime in YYYYMMDDHHMM format (hour must be 00/06/12/18)"
    )
    parser.add_argument(
        "--system", choices=list(SYSTEMS), default=None,
        help="HPC system profile (default: auto-detect from hostname)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Write the script but do not call sbatch"
    )
    args = parser.parse_args()

    try:
        dt = validate_datetime(args.datetime)
    except ValueError as exc:
        parser.error(str(exc))

    system = args.system or detect_system()
    print(f"Using system profile: '{system}'")

    submit(dt, system, dry_run=args.dry_run)


if __name__ == "__main__":
    main()