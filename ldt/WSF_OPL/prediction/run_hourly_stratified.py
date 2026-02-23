#!/usr/bin/env python
"""
Runner script for hourly NetCDF soil moisture prediction (stratified RF).

Processes hourly WSF NetCDF files organized in monthly subdirectories:
    INPUT_BASE/YYYYMM/*.nc

Output:
    OUTPUT_BASE/YYYYMM/ARFS_SM_WSFM_YYYYMMDDTHHMMSS.nc

Usage:
    # Process all months
    python run_hourly_stratified.py

    # Process specific month (for parallel SLURM jobs)
    python run_hourly_stratified.py --month 202404

    # Process a range of months
    python run_hourly_stratified.py --months 202401 202402 202403
"""

import os
import re
import glob
import sys
import argparse
import time as timer
from pathlib import Path
from predict_sm_stratified import StratifiedSoilMoisturePredictor

# =============================================================================
# Configuration
# =============================================================================
INPUT_BASE  = "/discover/nobackup/projects/usaf_lis/MET_FORCING/WSF/resampled"
OUTPUT_BASE = "/discover/nobackup/projects/usaf_lis/ejalilvand/data/WSF/hourly_sm"
MODEL_DIR   = "/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf/cluster_models"
CLUSTER_DIR = "/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf"

# File pattern to match input TB files within each monthly folder
FILE_PATTERN = "*.nc"

# Blending settings
TOP_K = 3
HARD  = False
TEMPERATURE = 1.0


# =============================================================================
# Helpers
# =============================================================================
def find_monthly_dirs(base_dir, months=None):
    """Find all YYYYMM subdirectories under base_dir."""
    all_dirs = sorted(glob.glob(os.path.join(base_dir, '[0-9]' * 6)))
    monthly = [(d, os.path.basename(d)) for d in all_dirs if os.path.isdir(d)]
    if months:
        monthly = [(d, name) for d, name in monthly if name in months]
    return monthly


def extract_file_key(path):
    """
    Extract a unique key (timestamp + overpass) from a filename.
    Returns a set of key strings, e.g., {'20240413T000000_ASC'}.
    """
    basename = os.path.basename(path).upper()
    timestamps = set()

    # Extract timestamp
    matches = re.findall(r'(\d{8})T(\d{6})', basename)
    for date, time in matches:
        timestamps.add(f'{date}T{time}')

    if not timestamps:
        matches = re.findall(r'(\d{8})T(\d{2})', basename)
        for date, hh in matches:
            timestamps.add(f'{date}T{hh}0000')

    if not timestamps:
        matches = re.findall(r'(\d{8})_T(\d{4})', basename)
        for date, hhmm in matches:
            timestamps.add(f'{date}T{hhmm}00')

    if not timestamps:
        matches = re.findall(r'(\d{8})_T(\d{2})', basename)
        for date, hh in matches:
            timestamps.add(f'{date}T{hh}0000')

    # Extract overpass direction
    if '_ASC' in basename or '.ASC' in basename:
        overpass = 'ASC'
    elif '_DES' in basename or '.DES' in basename:
        overpass = 'DES'
    else:
        overpass = None

    # Build keys: timestamp_overpass or just timestamp
    keys = set()
    for ts in timestamps:
        if overpass:
            keys.add(f'{ts}_{overpass}')
        else:
            keys.add(ts)
    return keys


def get_already_processed(output_dir):
    """Scan output dir for existing ARFS_SM_WSFM_*.nc file keys."""
    processed = set()
    pattern = os.path.join(output_dir, 'ARFS_SM_WSFM_*.nc')
    for f in glob.glob(pattern):
        keys = extract_file_key(f)
        processed.update(keys)
    return processed


# =============================================================================
# Main
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description='Run stratified RF soil moisture prediction on hourly WSF files')
    parser.add_argument('--month', type=str, default=None,
                        help='Single month to process, e.g., 202404')
    parser.add_argument('--months', nargs='+', type=str, default=None,
                        help='List of months to process, e.g., 202401 202402')
    args = parser.parse_args()

    # Determine which months to process
    if args.month:
        month_filter = [args.month]
    elif args.months:
        month_filter = args.months
    else:
        month_filter = None  # all

    # Load predictor once
    print("Initializing predictor...")
    init_t0 = timer.time()
    predictor = StratifiedSoilMoisturePredictor(
        model_dir=MODEL_DIR,
        cluster_dir=CLUSTER_DIR,
        top_k=TOP_K,
        temperature=TEMPERATURE,
        hard=HARD,
    )
    print(f"Predictor ready in {timer.time() - init_t0:.1f}s\n")

    # Find monthly directories
    monthly_dirs = find_monthly_dirs(INPUT_BASE, month_filter)
    if not monthly_dirs:
        print(f"No monthly directories found under {INPUT_BASE}")
        if month_filter:
            print(f"  Filter was: {month_filter}")
        sys.exit(1)

    print(f"Processing {len(monthly_dirs)} month(s): "
          f"{', '.join(name for _, name in monthly_dirs)}")

    grand_t0 = timer.time()
    total_processed = 0
    total_skipped = 0
    total_failed = 0

    for input_dir, month_name in monthly_dirs:
        print("\n" + "=" * 60)
        print(f"MONTH: {month_name}")
        print("=" * 60)

        output_dir = os.path.join(OUTPUT_BASE, month_name)
        os.makedirs(output_dir, exist_ok=True)

        # Find input files
        input_files = sorted(glob.glob(os.path.join(input_dir, FILE_PATTERN)))
        if not input_files:
            print(f"  No input files found, skipping")
            continue

        # Filter out already-processed
        already_done = get_already_processed(output_dir)
        todo = []
        skipped = 0
        for f in input_files:
            ts = extract_file_key(f)
            if ts and ts.issubset(already_done):
                skipped += 1
            else:
                todo.append(f)

        print(f"  Total files: {len(input_files)}, "
              f"to process: {len(todo)}, already done: {skipped}")

        if not todo:
            total_skipped += skipped
            continue

        # Process files sequentially — RF .predict() uses all cores internally
        month_t0 = timer.time()
        processed = 0
        failed = 0

        for i, input_path in enumerate(todo):
            input_stem = Path(input_path).stem
            t0 = timer.time()

            try:
                predictor.predict_from_hourly_nc(
                    input_path=input_path,
                    output_dir=output_dir,
                    apply_qc=True,
                )
                processed += 1
            except Exception as e:
                print(f'  [ERROR] {input_stem}: {e}')
                failed += 1

            # Progress every 10 files
            if (i + 1) % 10 == 0 or (i + 1) == len(todo):
                elapsed = timer.time() - month_t0
                rate = (i + 1) / elapsed * 60
                remaining = (len(todo) - i - 1) / rate if rate > 0 else 0
                print(f"  [{i+1:>4d}/{len(todo)}] "
                      f"{rate:.1f} files/min, "
                      f"~{remaining:.0f} min remaining")

        month_elapsed = timer.time() - month_t0
        rate = len(todo) / month_elapsed * 60 if month_elapsed > 0 else 0

        print(f"\n  {month_name}: processed={processed}, skipped={skipped}, "
              f"failed={failed}, time={month_elapsed/60:.1f} min, "
              f"rate={rate:.1f} files/min")

        total_processed += processed
        total_skipped += skipped
        total_failed += failed

    # Grand summary
    grand_elapsed = timer.time() - grand_t0
    print("\n" + "=" * 60)
    print("GRAND SUMMARY")
    print("=" * 60)
    print(f"  Months:    {len(monthly_dirs)}")
    print(f"  Processed: {total_processed}")
    print(f"  Skipped:   {total_skipped}")
    print(f"  Failed:    {total_failed}")
    print(f"  Total time: {grand_elapsed/60:.1f} min")
    if total_processed > 0:
        print(f"  Avg rate:  {total_processed / grand_elapsed * 60:.1f} files/min")
    print("=" * 60)


if __name__ == '__main__':
    main()
