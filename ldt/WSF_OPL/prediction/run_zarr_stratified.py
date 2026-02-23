#!/usr/bin/env python
"""
Runner script for monthly zarr soil moisture prediction (stratified RF).
Run: python run_zarr_stratified.py
"""

import os
import glob
from predict_sm_stratified import StratifiedSoilMoisturePredictor

# =============================================================================
# Configuration
# =============================================================================
INPUT_DIR  = "/discover/nobackup/projects/usaf_lis/MET_FORCING/WSF/resampled/daily_composite"
OUTPUT_DIR = "/discover/nobackup/projects/usaf_lis/ejalilvand/data/WSF/wsf_stratified_rf"
MODEL_DIR  = "/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf/cluster_models"
CLUSTER_DIR = "/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf"

YEARS  = [2024, 2025]
ORBITS = ['ASC', 'DES']

# Blending settings
TOP_K = 3           # Number of clusters to blend (set 1 or hard=True for no blending)
HARD  = False       # True = hard assignment, False = soft blending
TEMPERATURE = 1.0   # Softmax temperature (lower = sharper blending)


# =============================================================================
# Main
# =============================================================================
def main():
    # Load predictor once (model loading + env encoding is one-time cost)
    predictor = StratifiedSoilMoisturePredictor(
        model_dir=MODEL_DIR,
        cluster_dir=CLUSTER_DIR,
        top_k=TOP_K,
        temperature=TEMPERATURE,
        hard=HARD,
    )

    processed = 0
    skipped = 0
    failed = 0

    for year in YEARS:
        for month in range(1, 13):
            year_month = f'{year}{month:02d}'

            for orbit in ORBITS:
                input_name = f'{year_month}_{orbit}'
                input_path = f"{INPUT_DIR}/{input_name}.zarr"

                # Check if already processed
                existing = glob.glob(f'{OUTPUT_DIR}/wsf_sm_tb_only_*{input_name}*.zarr')
                if existing:
                    print(f'[SKIP] {input_name} already processed')
                    skipped += 1
                    continue

                # Check if input exists
                if not os.path.exists(input_path):
                    print(f'[SKIP] {input_name} input does not exist')
                    skipped += 1
                    continue

                # Process
                print(f'\n[PROCESSING] {input_name}...')
                try:
                    predictor.predict_from_zarr(
                        input_name=input_name,
                        input_dir=INPUT_DIR,
                        output_dir=OUTPUT_DIR,
                        apply_qc=True,
                    )
                    processed += 1
                except Exception as e:
                    print(f'[ERROR] {input_name} failed: {e}')
                    failed += 1

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Processed: {processed}")
    print(f"  Skipped:   {skipped}")
    print(f"  Failed:    {failed}")
    print("=" * 60)


if __name__ == '__main__':
    main()
