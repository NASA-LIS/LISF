#!/usr/bin/env python3
"""
WSF_data_reader.py  –  WSF OPL Resampling Worker

Usage: python WSF_data_reader.py --resample <START_DT> <END_DT>
[--config path/to/config.json]
   or: python WSF_data_reader.py --resample <END_DT>
"""
# pylint: disable=invalid-name
import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timedelta
from pathlib import Path

# =============================================================================
# RESAMPLING CONSTANTS
# =============================================================================
V522_LAST  = "2025041323"
DT_RE      = re.compile(r"WSFM_01_d(\d{8})_t(\d{2})")

def find_hours(base, start, end, lo="0"*10, hi="9"*10):
    """
        Finds unique 10-digit datetime hours from WSF-M NetCDF files
        Args:
            base (str or Path): The base directory to search recursively.
            start (str): The start datetime string
            end (str): The end datetime string
            lo (str, optional): An additional lower bound for the 10-digit dt.
            hi (str, optional): An additional upper bound for the 10-digit dt.

        Returns:
            set: A set of unique 10-character datetime strings (YYYYMMDDHH)
                 extracted from the matching filenames.
        """
    if not Path(base).is_dir():
        return set()
    start_10 = start[:10]
    end_10 = end[:10]
    # -----------------------------------------------------------------
    hours = set()
    for f in Path(base).rglob("*WSFM_01_d*_res_sdr.nc"):
        m = DT_RE.search(f.name)
        if m:
            dt = m.group(1) + m.group(2)  # This is 10 characters
            # Compare using the truncated 10-character boundaries
            if start_10 <= dt <= end_10 and lo <= dt <= hi:
                hours.add(dt)
    return hours

def src_dir(dt, v522_base, flat_base):
    """Return the correct source directory path based on the date."""
    if dt <= V522_LAST:
        return Path(v522_base) / dt[:4] / dt[4:6] / dt[6:8]
    return Path(flat_base) / dt[:6]

def write_config(template, dt, src, out):
    """Generate an LDT config file for a specific
    datetime by updating a template."""
    text = template.read_text(encoding='utf-8')
    for pattern, value in [
        (r"WSF valid date \(YYYYMMDDHH\):.*",
         f"WSF valid date (YYYYMMDDHH):            "
         f"{dt}"),
        (r"WSF input directory:.*",
         f"WSF input directory:                   "
         f" {src}"),
        (r"WSF output directory:.*",
         f"WSF output directory:                   "
         f"{out}"),
        (r"LDT diagnostic file:.*",
         f"LDT diagnostic file:                    "
         f"logs/01/ldtlog_{dt}"),
        (r"WSF filelist suffix number:.*",
         f"WSF filelist suffix number:             "
         f"{dt}"),
    ]:
        text = re.sub(pattern, value, text)
    cfg = Path(f"ldt.config.wsf.{dt}")
    cfg.write_text(text, encoding='utf-8')
    return cfg



def run_resampling(cfg, program, template, out_base, start, end, batch_size):
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    # pylint:disable=too-many-locals
    """
        Execute the resampling workflow for the provided
        datetimes in parallel batches.

        Args:
            cfg (dict): Configuration dictionary.
            program (Path): Path to the LDT executable.
            template (Path): Path to the LDT configuration template.
            out_base (Path): Base directory for the resampled output.
            start (str): Start datetime string.
            end (str): End datetime string.
            batch_size (int): Maximum number of parallel jobs to run.
        """
    print("\n=== OPL Resampling (WSF) ===")
    Path("logs/01").mkdir(exist_ok=True)

    hours  = find_hours(cfg["v522_sdr_base"], start, end, hi=V522_LAST)
    hours |= find_hours(cfg["raw_sdr_base"],  start, end,
                        lo=str(int(V522_LAST) + 1))
    pending = sorted(hours)

    print(f"  Hours found to process: {len(pending)}")
    if not pending:
        print("  No hours found in the specified timeframe.")
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

            # pylint: disable=consider-using-with
            proc = subprocess.Popen([str(program), str(cfg_file)],
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL)
            procs.append((dt, proc, cfg_file))

        for dt, proc, cfg_file in procs:
            rc  = proc.wait()
            log = Path(f"logs/01/ldtlog_{dt}")
            if rc != 0:
                print(f"    FAIL {dt}: exit code {rc}")
                n_fail += 1
            elif (log.exists() and "No WSF files found" in
                  log.read_text(encoding='utf-8')):
                print(f"    WARN {dt}: no input files found")
                n_fail += 1
            else:
                n_ok += 1
            cfg_file.unlink(missing_ok=True)
            Path(f"WSF_filelist_{dt}.dat").unlink(missing_ok=True)

        print(f"    ok={n_ok} fail={n_fail} skip={n_skip}")

    n_out = sum(1 for _ in out_base.rglob("WSF_SDR_resampled_*.nc"))
    print(f"\n  Resampled files: {n_out}  ok: {n_ok}  fail: {n_fail}  "
          f"skip: {n_skip}")
    if n_out == 0:
        sys.exit("ERROR: no resampled files produced — check logs/01/")

def main() -> None:
    """Parse arguments, load configuration,
    and execute the resampling workflow."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--resample", nargs='+',
                        required=True, metavar=("START_DT", "END_DT"),
                        help="Run OPL resampling (WSF only)."
                             " Provide either two args (START_DT END_DT) "
                             "or a single END_DT.")
    parser.add_argument("--config",
                        default="./SNIP_ops/config/SNIP_config.json",
                        help="Path to config")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=int(os.environ.get(
            "MAX_PARALLEL",
            max(10, (os.cpu_count() or 10) // 10 * 10))),
        help="Max parallel LDT jobs")

    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        sys.exit(f"ERROR: Configuration file '{args.config}' not found.")
    start_dt = None
    end_dt = None
    # Allow either: --resample START END  OR  --resample END
    if len(args.resample) == 2:
        start_dt, end_dt = args.resample
        cfg = json.loads(config_path.read_text(encoding='utf-8'))
    elif len(args.resample) == 1:
        end_dt = args.resample[0]
        cfg = json.loads(config_path.read_text(encoding='utf-8'))
        # Determine start offset (hours) from config; fall back to 6
        try:
            offset = int(cfg.get("start_offset_hours", 6))
        except (ValueError, TypeError):
            offset = 6
        try:
            dt_end = datetime.strptime(end_dt, "%Y%m%d%H%M")
        except ValueError as exc:
            sys.exit(f"ERROR: invalid end datetime '{end_dt}': {exc}")
        start_dt = (dt_end - timedelta(hours=offset)).strftime("%Y%m%d%H%M")
    else:
        parser.error("--resample requires 1 or 2 datetime arguments")
    program = Path(cfg["ldt"]).resolve()
    template = Path(cfg["ldt_config_template"]).resolve()
    out_base = Path(cfg["resampled_base"])

    run_resampling(cfg, program, template, out_base,
                   start_dt, end_dt, args.batch_size)

if __name__ == "__main__":
    main()
