#!/usr/bin/env python3
"""
WSF Soil Moisture Retrieval — Steps 1 & 2

  Step 1: OPL resampling  — batched parallel LDT
  Step 2: ML prediction   — sequential hours, all cores per hour

Usage:
  python wsf_opl.py START_DT END_DT [--config FILE] [--batch-size N] [--force]
                    [--step {1,2}] [--dry-run]

  START_DT / END_DT : period to process, format YYYYMMDDHH (e.g. 2024041300)
  --config          : path to the JSON config, see wsf_opl.json (default: ./pipeline_config.json)
  --batch-size      : max parallel LDT jobs for step 1 (default: cpu count)
  --force           : ignore checkpoint and reprocess all hours
  --step            : run only step 1 (resampling) or step 2 (prediction); omit to run both
  --dry-run         : show, month by month, which SDR files step 1 would read; run nothing

Step 1 input (step1_wsf_opl in the config):
  sdr_sources        ordered list of {"name", "path", "layout"}; layout is the
                     sub-folder for one day, with {yyyy} {mm} {dd} {yyyymm} {yyyymmdd}
                     (e.g. "{yyyy}/{mm}/{dd}" or "{yyyymm}")
  skip_file_prefixes file names starting with these are ignored (e.g. "temp")
  staging_dir        per-hour folders of links to the chosen files (LDT input),
                     named <YYYYMMDDHH>_<process id> and removed after each hour
  For every hour, each pass (a resampled SDR file, identified by its start
  time) is read in its newest processing version (_vNNNN_) found in any
  source; if several sources have that version, the first one listed wins.
  A file that runs past midnight (e.g. _t235900_e001059) is used for hour 23
  of its own day and hour 00 of the next day.
  Every chosen file is written to logs/sdr_selection.log.

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

  # See which SDR files would be read, without running anything
  python wsf_opl.py 2024041300 2026111200 --step 1 --dry-run

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
import shutil
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import NamedTuple

CHECKPOINT   = Path("logs/resampling_checkpoint.txt")
SELECTION_LOG = Path("logs/sdr_selection.log")
# Resampled SDR files the LDT resampler reads: *WSFM_01_d<YYYYMMDD>*_res_sdr.nc
SDR_RE = re.compile(
    r"WSFM_01_d(\d{8})_t(\d{6})_e(\d{6}).*_v(\d[\d.]*)_res_sdr\.nc$")
LAYOUT_KEYS = ("yyyy", "mm", "dd", "yyyymm", "yyyymmdd")


# =============================================================================
# Step 1 — choosing the SDR files for each hour
# =============================================================================
def version_key(tag):
    """Sortable processing version: v0602 -> (6, 0, 2), v051 -> (5, 1, 0), v6.1.0 -> (6, 1, 0)."""
    if "." in tag:
        return tuple(int(x) for x in tag.split(".") if x)
    if len(tag) == 4:
        return (int(tag[:2]), int(tag[2]), int(tag[3]))
    if len(tag) == 3:
        return (int(tag[:2]), int(tag[2]), 0)
    return (int(tag),)


class SdrFile(NamedTuple):
    """One resampled SDR file and where it was found."""
    rank: int        # position of its source in sdr_sources
    source: str
    path: Path
    day: str         # YYYYMMDD in the file name (_d)
    start: str       # HHMMSS (_t)
    end: str         # HHMMSS (_e)
    version: str     # processing version (_vNNNN_)


def belongs_to_hour(sdr, day, hour):
    """Whether a file has footprints in hour HH of day YYYYMMDD.

    A file that runs past midnight (e.g. _d20240601_t235900_e001059) covers
    hour 23 of its own day and hour 00 of the next day. The LDT resampler on
    its own would give it to hour 00 of its own day instead, almost a day
    early, so the driver decides here and LDT only sees the right files.
    """
    h0, h1 = int(sdr.start[:2]), int(sdr.end[:2])
    wraps = h1 < h0
    if sdr.day != day:                        # previous day, past midnight
        return wraps and hour <= h1
    return hour >= h0 if wraps else h0 <= hour <= h1


def link_name(sdr, day):
    """Name of the link LDT reads: a previous-day file gets this day's _d date,
    because LDT only lists *WSFM_01_d<day>*_res_sdr.nc."""
    if sdr.day == day:
        return sdr.path.name
    return sdr.path.name.replace(f"_d{sdr.day}_", f"_d{day}_", 1)


class SdrSources:
    """Ordered SDR folders (from the config) with cached directory listings."""

    def __init__(self, sources, skip_prefixes):
        self.sources = []
        for i, src in enumerate(sources):
            for key in ("name", "path", "layout"):
                if key not in src:
                    sys.exit(f"ERROR: sdr_sources[{i}] has no '{key}'")
            self.sources.append((src["name"], Path(src["path"]), src["layout"]))
        self.skip_prefixes = tuple(skip_prefixes)
        self._listing = {}

    def _list(self, folder):
        if folder not in self._listing:
            try:
                self._listing[folder] = os.listdir(folder)
            except FileNotFoundError:
                self._listing[folder] = []
        return self._listing[folder]

    def describe(self):
        """Lines describing the sources, in reading order."""
        return [f"{name}: {path}/{layout}" for name, path, layout in self.sources]

    def _files_named(self, day):
        """SdrFile entries whose name carries the date YYYYMMDD."""
        fields = {"yyyy": day[:4], "mm": day[4:6], "dd": day[6:8],
                  "yyyymm": day[:6], "yyyymmdd": day}
        found = []
        for rank, (name, path, layout) in enumerate(self.sources):
            folder = path / layout.format(**fields)
            for fname in self._list(folder):
                if fname.startswith(self.skip_prefixes):
                    continue
                m = SDR_RE.search(fname)
                if m and m.group(1) == day:
                    found.append(SdrFile(rank, name, folder / fname, day, m.group(2),
                                         m.group(3), m.group(4)))
        return found

    def files_for_day(self, day):
        """Files with footprints on YYYYMMDD: that day's files, plus the previous
        day's files that run past midnight."""
        prev = (datetime.strptime(day, "%Y%m%d") - timedelta(days=1)).strftime("%Y%m%d")
        late = [f for f in self._files_named(prev) if f.end[:2] < f.start[:2]]
        return self._files_named(day) + late


def select_hour(day_files, day, hour):
    """Files to resample for one hour: the newest version of each pass.

    A pass is identified by its date and start time (to the minute). Among
    the files of a pass, only the newest processing version is kept, from the
    first source (in config order) that has it; several copies of that
    version (e.g. other ground stations) are all kept, as LDT chooses
    between them.
    """
    passes = {}
    for sdr in day_files:
        if belongs_to_hour(sdr, day, hour):
            passes.setdefault((sdr.day, sdr.start[:4]), []).append(sdr)
    chosen, dropped = [], 0
    for files in passes.values():
        best = max(version_key(f.version) for f in files)
        newest = [f for f in files if version_key(f.version) == best]
        first_rank = min(f.rank for f in newest)
        keep = [f for f in newest if f.rank == first_rank]
        chosen.extend(keep)
        dropped += len(files) - len(keep)
    return sorted(chosen, key=lambda f: (f.day, f.path.name)), dropped


def plan_resampling(sources, start, end):
    """{YYYYMMDDHH: [chosen files]} for every hour with data, plus monthly stats."""
    plan, stats = {}, {}
    day = datetime.strptime(start[:8], "%Y%m%d")
    while day <= datetime.strptime(end[:8], "%Y%m%d"):
        yyyymmdd = day.strftime("%Y%m%d")
        day += timedelta(days=1)
        day_files = sources.files_for_day(yyyymmdd)
        month = stats.setdefault(yyyymmdd[:6], {"hours": 0, "files": 0, "dropped": 0,
                                                "versions": {}, "sources": {}})
        for hour in range(24):
            dt = f"{yyyymmdd}{hour:02d}"
            if not start <= dt <= end:
                continue
            chosen, dropped = select_hour(day_files, yyyymmdd, hour)
            if not chosen:
                continue
            plan[dt] = chosen
            month["hours"] += 1
            month["files"] += len(chosen)
            month["dropped"] += dropped
            for f in chosen:
                month["versions"][f"v{f.version}"] = month["versions"].get(f"v{f.version}", 0) + 1
                month["sources"][f.source] = month["sources"].get(f.source, 0) + 1
    return plan, stats


def print_plan(stats):
    """Monthly summary of the files chosen for resampling."""
    print(f"  {'month':<8}{'hours':>6}{'files':>8}{'older/duplicate dropped':>25}"
          "   versions | sources")
    for month, st in sorted(stats.items()):
        vers = ", ".join(f"{k} {v}" for k, v in sorted(st["versions"].items()))
        srcs = ", ".join(f"{k} {v}" for k, v in st["sources"].items())
        print(f"  {month:<8}{st['hours']:>6}{st['files']:>8}{st['dropped']:>25}"
              f"   {vers} | {srcs}")


def stage_hour(staging_dir, dt, chosen):
    """Folder of links to the chosen files: the LDT input directory for this hour.

    The folder name includes the process id, so two runs that share an hour
    (e.g. monthly jobs whose periods meet at 00Z of the 1st) never touch
    each other's folder.
    """
    folder = Path(staging_dir) / f"{dt}_{os.getpid()}"
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    for f in chosen:
        (folder / link_name(f, dt[:8])).symlink_to(f.path.resolve())
    return folder


def write_config(template, dt, src, out):
    """Copy the LDT template and substitute the per-hour fields."""
    text = template.read_text(encoding='utf-8')
    for pattern, value in [
        (r"WSF valid date \(YYYYMMDDHH\):.*", f"WSF valid date (YYYYMMDDHH):            {dt}"),
        (r"WSF input directory:.*",            f"WSF input directory:                    {src}"),
        (r"WSF output directory:.*",           f"WSF output directory:                   {out}"),
        (r"LDT diagnostic file:.*",
         f"LDT diagnostic file:                    logs/ldtlog_{dt}"),
        (r"WSF filelist suffix number:.*",     f"WSF filelist suffix number:             {dt}"),
    ]:
        text = re.sub(pattern, value, text)
    cfg = Path(f"ldt.config.wsf.{dt}")
    cfg.write_text(text, encoding='utf-8')
    return cfg


HINTS_SHOWN = set()   # failure hints already printed in this run


def check_period(parser, start, end):
    """Stop with a usage error unless START_DT and END_DT are valid YYYYMMDDHH.

    A date without the hour (YYYYMMDD) used to be accepted, and the whole end
    day was then left out of the period.
    """
    for name, value in (("START_DT", start), ("END_DT", end)):
        if not re.fullmatch(r"\d{10}", value):
            parser.error(f"{name} must be YYYYMMDDHH (10 digits), not {value!r}")
        try:
            datetime.strptime(value, "%Y%m%d%H")
        except ValueError:
            parser.error(f"{name} is not a valid date and hour: {value}")
    if start > end:
        parser.error(f"START_DT {start} is after END_DT {end}")


def ldt_stdout(dt):
    """File that receives LDT's screen output for hour dt (kept only on failure)."""
    return Path(f"logs/ldt_stdout_{dt}.txt")


def launch_ldt(program, cfg_file, dt):
    """Start LDT for one hour, with its screen output going to ldt_stdout(dt)."""
    with open(ldt_stdout(dt), "w", encoding="utf-8") as out:
        # The process is waited on later, with the rest of its batch; it keeps
        # its own copy of the output file after this one is closed.
        return subprocess.Popen(  # pylint: disable=consider-using-with
            [str(program), str(cfg_file)], stdout=out, stderr=subprocess.STDOUT)


def last_line(path, nbytes=4096):
    """Last non-empty line of a text file ('' if the file is missing or empty)."""
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            fh.seek(max(0, fh.tell() - nbytes))
            lines = fh.read().decode(errors="replace").splitlines()
    except FileNotFoundError:
        return ""
    return next((line.strip()[:200] for line in reversed(lines) if line.strip()), "")


def report_failure(dt, rc):
    """Print why LDT failed for hour dt; its screen output stays in logs/."""
    print(f"    FAIL {dt}: exit code {rc}")
    last = last_line(ldt_stdout(dt))
    if last:
        print(f"         {last}")
    if rc == 127 and rc not in HINTS_SHOWN:
        HINTS_SHOWN.add(rc)
        print("         LDT could not start (a library or command was not found):"
              " load the LDT modules, or submit through the batch script.")
    print(f"         LDT screen output: {ldt_stdout(dt)}")


# =============================================================================
# Step 1 — OPL Resampling
# =============================================================================
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements,too-many-branches
def step1_resample(s1, program, template, out_base, start, end, batch_size, force,
                   dry_run=False):
    """Run LDT OPL resampling in parallel batches over all pending hours."""
    print("\n=== STEP 1: OPL Resampling ===")

    sources = SdrSources(s1["sdr_sources"], s1["skip_file_prefixes"])
    print("  SDR sources (a pass is read in its newest version; ties go to the "
          "first source):")
    for line in sources.describe():
        print(f"    {line}")
    plan, stats = plan_resampling(sources, start, end)
    print_plan(stats)
    if not plan:
        sys.exit(f"ERROR: no SDR files for {start} to {end} in the SDR sources above")
    if dry_run:
        print("  Dry run: nothing resampled.")
        return

    Path("logs").mkdir(exist_ok=True)
    CHECKPOINT.touch()
    done    = set() if force else set(CHECKPOINT.read_text(encoding='utf-8').split())
    pending = sorted(set(plan) - done)

    print(f"  Hours with data: {len(plan)}  done: {len(done & set(plan))}  "
          f"pending: {len(pending)}")

    if not pending:
        print("  All hours already processed.")
        return

    n_ok = n_fail = 0
    n_batches = (len(pending) + batch_size - 1) // batch_size

    for b, i in enumerate(range(0, len(pending), batch_size)):
        batch = pending[i : i + batch_size]
        print(f"\n  Batch {b+1}/{n_batches} ({len(batch)} hours)")

        procs = []
        with SELECTION_LOG.open("a", encoding='utf-8') as log:
            for dt in batch:
                for f in plan[dt]:
                    note = " (previous day, past midnight)" if f.day != dt[:8] else ""
                    log.write(f"{dt} {f.source} v{f.version} {f.path}{note}\n")
                sd = stage_hour(s1["staging_dir"], dt, plan[dt])
                od = out_base / dt[:6]
                od.mkdir(parents=True, exist_ok=True)
                cfg_file = write_config(template, dt, sd.resolve(), od)
                proc = launch_ldt(program, cfg_file, dt)
                procs.append((dt, proc, cfg_file, sd))

        for dt, proc, cfg_file, sd in procs:
            rc  = proc.wait()
            log = Path(f"logs/ldtlog_{dt}")
            if rc != 0:
                report_failure(dt, rc)
                n_fail += 1
            elif log.exists() and "No WSF files found" in log.read_text(encoding='utf-8'):
                print(f"    WARN {dt}: no input files found")
                n_fail += 1
            else:
                with CHECKPOINT.open("a", encoding='utf-8') as fh:
                    fh.write(dt + "\n")
                n_ok += 1
            if rc == 0:
                ldt_stdout(dt).unlink(missing_ok=True)
            cfg_file.unlink(missing_ok=True)
            Path(f"WSF_filelist_{dt}.dat").unlink(missing_ok=True)
            shutil.rmtree(sd, ignore_errors=True)

        print(f"    ok={n_ok} fail={n_fail}")

    n_out = sum(1 for _ in out_base.rglob("WSF_SDR_resampled_*.nc"))
    print(f"\n  Resampled files: {n_out}  ok: {n_ok}  fail: {n_fail}")
    if n_ok == 0:
        sys.exit(f"ERROR: none of the {len(pending)} pending hours was resampled "
                 "(see the lines above and logs/)")
    if n_out == 0:
        sys.exit("ERROR: no resampled files produced — check logs/")
# pylint: enable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements,too-many-branches


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
    # Dynamic import: script location is predict_script in the config.
    from predict_sm_stratified import StratifiedSoilMoisturePredictor
    # pylint: enable=import-outside-toplevel

    required = ["model_dir", "cluster_dir", "env_paths", "blending", "top_k"]
    if s2.get("blending") == "distance":
        required.append("temperature")
    missing = [key for key in required if key not in s2]
    if missing:
        sys.exit("ERROR: step2_prediction in the config file has no "
                 + ", ".join(missing))

    predictor = StratifiedSoilMoisturePredictor(
        model_dir=s2["model_dir"],
        cluster_dir=s2["cluster_dir"],
        env_paths=s2["env_paths"],
        blending=s2["blending"],
        top_k=int(s2["top_k"]),
        temperature=float(s2["temperature"]) if "temperature" in s2 else None,
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
        # Exact output name (as written by predict_from_hourly_nc). A wildcard
        # match is not safe here: "*1000*" also matches the 01Z file
        # (...T010000...), so 10Z and 20Z were skipped once 01Z and 02Z existed.
        if not force and (sm_out / f"ARFS_SM_WSFM_{date}T{hhmm}00_{orbit}.nc").exists():
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

    print(f"\n  Done: {n_ok} ok, {n_fail} failed, {n_skip} already there (skipped)")
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
    parser.add_argument("--dry-run", action="store_true",
                        help="Show which SDR files step 1 would read; run nothing")
    args = parser.parse_args()

    if not args.start_dt or not args.end_dt:
        parser.error("START_DT and END_DT required")
    check_period(parser, args.start_dt, args.end_dt)

    cfg      = json.loads(Path(args.config).read_text(encoding='utf-8'))
    s1       = cfg["step1_wsf_opl"]
    if "sdr_sources" not in s1 and ("v522_sdr_base" in s1 or "raw_sdr_base" in s1):
        sys.exit("ERROR: step1_wsf_opl uses v522_sdr_base/raw_sdr_base, which are no "
                 "longer read. Replace them with sdr_sources, skip_file_prefixes and "
                 "staging_dir (see wsf_opl.json).")
    missing = [key for key in ("ldt_config_template", "resampled_base", "sdr_sources",
                               "skip_file_prefixes", "staging_dir") if key not in s1]
    if missing:
        sys.exit("ERROR: step1_wsf_opl in the config file has no " + ", ".join(missing))
    program  = Path(cfg["executables"]["ldt"]).resolve()
    template = Path(s1["ldt_config_template"]).resolve()
    out_base = Path(s1["resampled_base"])

    print(f"Period: {args.start_dt} → {args.end_dt}")

    if args.step != 2:
        step1_resample(s1, program, template, out_base,
                       args.start_dt, args.end_dt, args.batch_size, args.force,
                       dry_run=args.dry_run)
    if args.dry_run:
        return

    if args.step != 1:
        step2_predict(cfg["step2_prediction"], out_base, args.start_dt, args.end_dt, args.force)

    n_sm = sum(
        1 for _ in
        Path(cfg["step2_prediction"]["sm_output_dir"]).rglob("ARFS_SM_WSFM_*.nc")
    )
    print(f"\n=== DONE — SM files: {n_sm} ===")


if __name__ == "__main__":
    main()
