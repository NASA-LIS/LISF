#!/usr/bin/env python
# coding: utf-8
"""
Stratified RF Inference with Soft Blending
============================================
Predicts soil moisture from WSF TB using per-cluster RF models.
Uses GMM / distance-weighted soft assignment to eliminate boundary artifacts.

Usage:
    python predict_blended_sm.py --date 2025-04-15
    python predict_blended_sm.py --date-range 2025-04-01 2025-04-30
    python predict_blended_sm.py --date 2025-04-15 --top-k 3 --hard   # hard assignment (no blending)
    python predict_blended_sm.py --date 2025-04-15 --top-k 5          # blend top-5 clusters

Soft Blending Strategy:
    For each pixel, instead of picking ONE cluster model:
    1. Compute membership weights from GMM posterior (or KMeans centroid distances)
    2. Predict SM using the top-N most probable cluster models
    3. Blend: SM = Σ(weight_i * SM_pred_i) / Σ(weight_i)

    This smooths transitions at cluster boundaries while keeping
    the specialized cluster behavior in cluster cores.
"""

import os
import sys
import glob
import json
import time
import argparse
import warnings
import gc

import numpy as np
import pandas as pd
import xarray as xr
import joblib
from datetime import datetime, timedelta
from scipy.special import softmax

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION — must match train_stratified_rf.py
# =============================================================================
WSF_ZARR_DIR = '/discover/nobackup/projects/usaf_lis/MET_FORCING/WSF/resampled/daily_composite'
CLUSTER_DIR  = '/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf'
MODEL_DIR    = '/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf/cluster_models'
OUTPUT_DIR   = '/discover/nobackup/ejalilva/dev/wsf_ml_sm/predictions/stratified_rf'

ENV_PATHS = {
    'land_cover': '/discover/nobackup/ejalilva/data/LC/lc_arfs.nc',
    'climate':    '/discover/nobackup/ejalilva/data/LC/clm_usaf.nc',
    'clay':       '/discover/nobackup/ejalilva/data/soil/sg_clay_usaf.nc',
    'sand':       '/discover/nobackup/ejalilva/data/soil/sg_sand_usaf.nc',
    'silt':       '/discover/nobackup/ejalilva/data/soil/sg_silt_usaf.nc',
}

WSF_TB_VARS = [
    'TB_10V', 'TB_10H',
    'TB_18V', 'TB_18H',
    'TB_23V',
    'TB_36V', 'TB_36H',
    'TB_89V', 'TB_89H',
]

KOPPEN_TO_GROUP = {
    1: 1, 2: 1, 3: 1,
    4: 2, 5: 2,
    6: 3, 7: 3,
    8: 4, 9: 4, 10: 4,
    11: 5, 12: 5, 13: 5,
    14: 6, 15: 6, 16: 6,
    17: 7, 18: 7, 19: 7, 20: 7,
    21: 8, 22: 8, 23: 8, 24: 8,
    25: 9, 26: 9, 27: 9, 28: 9,
    29: 10, 30: 10,
}
LC_EXCLUDE = [0, 15, 17]
CLIMATE_EXCLUDE = [10]
CATEGORICAL_FEATURES = ['land_cover', 'climate_group']
CONTINUOUS_FEATURES = ['clay', 'sand', 'silt']


# =============================================================================
# ARGUMENT PARSING
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(description='Predict SM with stratified RF + soft blending')
    p.add_argument('--date', type=str, default=None, help='Single date YYYY-MM-DD')
    p.add_argument('--date-range', nargs=2, type=str, default=None,
                   help='Date range: start end (YYYY-MM-DD YYYY-MM-DD)')
    p.add_argument('--top-k', type=int, default=3,
                   help='Number of top clusters to blend (default: 3)')
    p.add_argument('--hard', action='store_true',
                   help='Use hard cluster assignment (no blending)')
    p.add_argument('--temperature', type=float, default=1.0,
                   help='Softmax temperature for distance-based blending (lower = sharper)')
    p.add_argument('--model-dir', type=str, default=MODEL_DIR)
    p.add_argument('--output-dir', type=str, default=OUTPUT_DIR)
    p.add_argument('--n-jobs', type=int, default=-1)
    return p.parse_args()


# =============================================================================
# LOAD MODELS + CLUSTER INFRASTRUCTURE
# =============================================================================
def load_all_models(model_dir):
    """Load all cluster RF models and the config."""
    print("\n" + "=" * 60)
    print("LOADING MODELS")
    print("=" * 60)

    config_path = os.path.join(model_dir, 'stratified_config.json')
    with open(config_path) as f:
        config = json.load(f)

    k = config['k']
    features = config['features']
    print(f"  K={k}, features={features}")

    models = {}

    # Global fallback
    global_path = os.path.join(model_dir, 'rf_global.pkl')
    models['global'] = joblib.load(global_path)
    print(f"  Loaded: global model")

    # Per-cluster
    for c in range(k):
        path = os.path.join(model_dir, f'rf_cluster_{c:02d}.pkl')
        if os.path.exists(path):
            models[c] = joblib.load(path)
        else:
            models[c] = None  # will use global
            print(f"  Cluster {c}: using global fallback")

    print(f"  Total models: {sum(1 for v in models.values() if v is not None)}")

    return models, config, k, features


def load_cluster_infrastructure(cluster_dir, model_dir, k):
    """
    Load everything needed for soft assignment:
    - Cluster grid (for hard assignment)
    - Clustering model (KMeans centroids for distance-based soft assignment)
    - GMM (if available, for posterior-based soft assignment)
    - Environmental layer encoders
    """
    print("\n" + "=" * 60)
    print("LOADING CLUSTER INFRASTRUCTURE")
    print("=" * 60)

    # Cluster grid
    cluster_da = xr.open_dataarray(
        os.path.join(cluster_dir, f'environmental_clusters_k{k}.nc')
    )
    print(f"  Cluster grid: {cluster_da.shape}")

    # Clustering model (KMeans + encoders)
    cluster_model = joblib.load(
        os.path.join(cluster_dir, f'clustering_model_k{k}.pkl')
    )
    print(f"  KMeans centroids: {cluster_model['kmeans'].cluster_centers_.shape}")

    # GMM (optional)
    gmm_path = os.path.join(model_dir, f'gmm_k{k}.pkl')
    if os.path.exists(gmm_path):
        gmm = joblib.load(gmm_path)
        print(f"  GMM loaded: {gmm_path}")
    else:
        gmm = None
        print(f"  No GMM found — will use distance-based soft assignment")

    return cluster_da, cluster_model, gmm


# =============================================================================
# ENCODE ENVIRONMENTAL FEATURES (for soft assignment)
# =============================================================================
def encode_pixel_features(layers, encoders, feature_names):
    """
    Encode the environmental layers to the same feature space used for clustering.
    Returns: X_encoded (n_pixels, n_features), valid_mask (n_pixels,)
    """
    print("  Encoding environmental features...")

    ref = list(layers.values())[0]
    n_lat, n_lon = ref.shape
    n_total = n_lat * n_lon

    # Stack raw features
    raw_arrays = []
    for name in feature_names:
        arr = layers[name].values
        if arr.ndim == 3:
            arr = arr[0]
        raw_arrays.append(arr.flatten())
    X_raw = np.column_stack(raw_arrays)

    # Encode
    encoded_parts = []
    for i, feat_name in enumerate(feature_names):
        col = X_raw[:, i].astype(np.float64)
        enc = encoders.get(feat_name, None)

        if enc is None:
            encoded_parts.append(col.reshape(-1, 1))
            continue

        if enc['type'] == 'onehot':
            cats = enc['categories']
            onehot = np.zeros((n_total, len(cats)), dtype=np.float64)
            for j, cat_val in enumerate(cats):
                onehot[:, j] = (col == cat_val).astype(np.float64)
            onehot[np.isnan(col), :] = np.nan
            encoded_parts.append(onehot)

        elif enc['type'] == 'continuous':
            scaled = np.full(n_total, np.nan, dtype=np.float64)
            valid = ~np.isnan(col) & (col >= 0) & (col <= 100)
            if enc['std'] > 0:
                scaled[valid] = (col[valid] - enc['mean']) / enc['std']
            encoded_parts.append(scaled.reshape(-1, 1))

    X_encoded = np.hstack(encoded_parts)
    valid_mask = ~np.isnan(X_encoded).any(axis=1)

    print(f"    Encoded shape: {X_encoded.shape}")
    print(f"    Valid pixels: {valid_mask.sum():,} / {n_total:,}")

    return X_encoded, valid_mask


def load_and_encode_env_layers(cluster_model):
    """Load environmental layers and encode them for cluster assignment."""

    layers = {}
    for name, path in ENV_PATHS.items():
        if not os.path.exists(path):
            continue
        ds = xr.open_dataset(path)
        var_name = list(ds.data_vars)[0]
        da = ds[var_name]
        if da.ndim == 3:
            da = da.squeeze()
        layers[name] = da

    # Group climate
    if 'climate' in layers:
        clm_arr = layers['climate'].values.copy()
        clm_grouped = np.full_like(clm_arr, np.nan, dtype=np.float64)
        for orig, grp in KOPPEN_TO_GROUP.items():
            clm_grouped[clm_arr == orig] = grp
        layers['climate_group'] = xr.DataArray(
            clm_grouped, dims=layers['climate'].dims, coords=layers['climate'].coords
        )
        del layers['climate']

    X_encoded, valid_mask = encode_pixel_features(
        layers, cluster_model['encoders'], cluster_model['feature_names']
    )

    ref = list(layers.values())[0]
    return X_encoded, valid_mask, ref.lat.values, ref.lon.values


# =============================================================================
# SOFT ASSIGNMENT: COMPUTE PIXEL WEIGHTS
# =============================================================================
def compute_soft_weights(X_encoded, valid_mask, cluster_model, gmm, top_k, temperature):
    """
    Compute soft cluster membership weights for each pixel.

    Returns:
        top_clusters: (n_valid, top_k) — cluster indices
        top_weights:  (n_valid, top_k) — normalized weights summing to 1
    """
    print("\n  Computing soft cluster weights...")
    t0 = time.time()

    X_valid = X_encoded[valid_mask]
    n_valid = X_valid.shape[0]

    if gmm is not None:
        # --- Method A: GMM posterior probabilities ---
        print("    Using GMM posterior probabilities")
        probs = gmm.predict_proba(X_valid)  # (n_valid, k)

    else:
        # --- Method B: Distance-based soft assignment ---
        print("    Using distance-based assignment (inverse squared distance)")
        kmeans = cluster_model['kmeans']
        centroids = kmeans.cluster_centers_  # (k, n_features)

        # Compute distances to all centroids
        # Chunked to avoid memory blow-up
        chunk_size = 100_000
        k = centroids.shape[0]
        probs = np.zeros((n_valid, k), dtype=np.float64)

        for start in range(0, n_valid, chunk_size):
            end = min(start + chunk_size, n_valid)
            X_chunk = X_valid[start:end]

            # Euclidean distance to each centroid
            dists = np.zeros((end - start, k))
            for c in range(k):
                dists[:, c] = np.sqrt(((X_chunk - centroids[c]) ** 2).sum(axis=1))

            # Convert to weights: softmax of negative distances
            probs[start:end] = softmax(-dists / temperature, axis=1)

    # Get top-K clusters and their weights
    top_clusters = np.argsort(probs, axis=1)[:, -top_k:][:, ::-1]  # descending
    top_weights = np.take_along_axis(probs, top_clusters, axis=1)

    # Renormalize weights to sum to 1
    weight_sums = top_weights.sum(axis=1, keepdims=True)
    weight_sums[weight_sums == 0] = 1.0
    top_weights = top_weights / weight_sums

    print(f"    Time: {time.time() - t0:.1f}s")
    print(f"    Top-{top_k} weight stats:")
    print(f"      Dominant cluster weight: mean={top_weights[:, 0].mean():.3f}, "
          f"min={top_weights[:, 0].min():.3f}")
    if top_k > 1:
        print(f"      2nd cluster weight:     mean={top_weights[:, 1].mean():.3f}")

    return top_clusters, top_weights


# =============================================================================
# PREDICT SOIL MOISTURE
# =============================================================================
def load_wsf_for_date(date_str):
    """Load WSF TB data for a specific date."""

    date = pd.Timestamp(date_str)
    year = date.year
    month = date.month

    # Find the zarr file for this month
    pattern = f'{WSF_ZARR_DIR}/{year}{month:02d}*_DES.zarr'
    files = sorted(glob.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No WSF file found for {pattern}")

    wsf_ds = xr.open_mfdataset(files, engine='zarr', chunks='auto',
                                combine='by_coords', decode_timedelta=True)
    wsf_ds['time'] = wsf_ds.time + pd.Timedelta(hours=6)

    # Select the specific date
    wsf_day = wsf_ds.sel(time=date.strftime('%Y-%m-%d'), method='nearest')

    # Quality mask
    if 'QUALITY_FLAG' in wsf_day.data_vars:
        ocean_mask = (wsf_day.QUALITY_FLAG.astype("int") & 1) > 0
        wsf_day = wsf_day.where(~ocean_mask)

    # Check available vars
    available_vars = [v for v in WSF_TB_VARS if v in wsf_day.data_vars]

    return wsf_day, available_vars


def predict_single_date(date_str, models, k, top_clusters, top_weights,
                        valid_mask, env_lats, env_lons, hard_mode=False):
    """
    Predict soil moisture for one date with soft blending.
    """
    print(f"\n  Predicting: {date_str}")
    t0 = time.time()

    # Load WSF TB
    wsf_day, available_vars = load_wsf_for_date(date_str)
    n_lat_wsf = len(wsf_day.lat)
    n_lon_wsf = len(wsf_day.lon)

    # Stack TB features
    tb_list = []
    for var in available_vars:
        arr = wsf_day[var].values
        if arr.ndim == 3:
            arr = arr[0]  # squeeze time if present
        tb_list.append(arr.flatten())
    X_tb = np.column_stack(tb_list)  # (n_pixels_wsf, n_features)

    # We need to match WSF pixels to the environmental grid
    # If grids match (same lat/lon), direct index mapping
    # If grids differ, use nearest-neighbor
    wsf_lats = wsf_day.lat.values
    wsf_lons = wsf_day.lon.values

    # Build mapping: for each WSF pixel, find nearest env grid index
    # This assumes both grids are regular — use searchsorted for speed
    lat_indices = np.searchsorted(env_lats, wsf_lats) - 1
    lon_indices = np.searchsorted(env_lons, wsf_lons) - 1
    lat_indices = np.clip(lat_indices, 0, len(env_lats) - 1)
    lon_indices = np.clip(lon_indices, 0, len(env_lons) - 1)

    n_env_lon = len(env_lons)

    # Map WSF pixel (i, j) → env flat index
    wsf_lat_grid, wsf_lon_grid = np.meshgrid(
        np.arange(n_lat_wsf), np.arange(n_lon_wsf), indexing='ij'
    )
    env_flat_indices = lat_indices[wsf_lat_grid.flatten()] * n_env_lon + lon_indices[wsf_lon_grid.flatten()]

    # Initialize output
    n_wsf_pixels = n_lat_wsf * n_lon_wsf
    sm_pred = np.full(n_wsf_pixels, np.nan)

    # Mask for valid TB data
    tb_valid = ~np.isnan(X_tb).any(axis=1)

    # Mask for valid env data
    env_valid_for_wsf = valid_mask[env_flat_indices]

    # Combined mask
    predict_mask = tb_valid & env_valid_for_wsf

    if predict_mask.sum() == 0:
        print(f"    No valid pixels for {date_str}")
        return None

    X_predict = X_tb[predict_mask]

    # Map to valid env indices (for weight lookup)
    # valid_mask cumsum gives the index into the "valid-only" arrays
    valid_cumsum = np.cumsum(valid_mask) - 1  # maps flat env index → valid index
    env_valid_idx = valid_cumsum[env_flat_indices[predict_mask]]

    if hard_mode:
        # ----- Hard assignment: single cluster model -----
        cluster_ids = top_clusters[env_valid_idx, 0].astype(int)
        unique_clusters = np.unique(cluster_ids)

        for c in unique_clusters:
            c_mask = cluster_ids == c
            model = models.get(c, None)
            if model is None:
                model = models['global']
            preds = model.predict(X_predict[c_mask])
            sm_pred_subset = np.full(predict_mask.sum(), np.nan)
            if 'sm_pred_subset' not in dir():
                sm_pred_subset = np.full(predict_mask.sum(), np.nan)
            sm_pred_subset = np.full(predict_mask.sum(), 0.0)

        # Redo hard mode more cleanly
        sm_pred_valid = np.zeros(predict_mask.sum())
        cluster_ids = top_clusters[env_valid_idx, 0].astype(int)
        for c in np.unique(cluster_ids):
            c_mask = cluster_ids == c
            model = models.get(int(c), None) or models['global']
            sm_pred_valid[c_mask] = model.predict(X_predict[c_mask])

    else:
        # ----- Soft blending: weighted average of top-K models -----
        n_predict = predict_mask.sum()
        top_k = top_clusters.shape[1]

        sm_pred_valid = np.zeros(n_predict)

        for rank in range(top_k):
            cluster_ids = top_clusters[env_valid_idx, rank].astype(int)
            weights = top_weights[env_valid_idx, rank]

            # Predict per-cluster
            preds_this_rank = np.zeros(n_predict)
            for c in np.unique(cluster_ids):
                c_mask = cluster_ids == c
                model = models.get(int(c), None) or models['global']
                preds_this_rank[c_mask] = model.predict(X_predict[c_mask])

            sm_pred_valid += weights * preds_this_rank

    # Clip and place back
    sm_pred_valid = np.clip(sm_pred_valid, 0.0, 0.6)
    sm_pred[predict_mask] = sm_pred_valid

    # Reshape to grid
    sm_grid = sm_pred.reshape(n_lat_wsf, n_lon_wsf)

    # Create DataArray
    sm_da = xr.DataArray(
        sm_grid,
        dims=['lat', 'lon'],
        coords={'lat': wsf_lats, 'lon': wsf_lons},
        name='soil_moisture',
        attrs={
            'units': 'm³/m³',
            'long_name': 'Soil Moisture (WSF stratified RF)',
            'date': date_str,
            'blending': 'soft' if not hard_mode else 'hard',
            'top_k': int(top_clusters.shape[1]),
        },
    )

    n_valid = np.isfinite(sm_grid).sum()
    print(f"    Valid pixels: {n_valid:,} / {n_wsf_pixels:,}")
    print(f"    SM range: [{np.nanmin(sm_grid):.3f}, {np.nanmax(sm_grid):.3f}]")
    print(f"    Time: {time.time() - t0:.1f}s")

    return sm_da


# =============================================================================
# MAIN
# =============================================================================
def main():
    args = parse_args()
    total_start = time.time()

    print("=" * 70)
    print("STRATIFIED RF PREDICTION WITH SOFT BLENDING")
    print("=" * 70)
    print(f"  top_k: {args.top_k}")
    print(f"  hard mode: {args.hard}")
    print(f"  temperature: {args.temperature}")

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # ----- Load models -----
    models, config, k, features = load_all_models(args.model_dir)

    # ----- Load cluster infrastructure -----
    cluster_da, cluster_model, gmm = load_cluster_infrastructure(
        CLUSTER_DIR, args.model_dir, k
    )

    # ----- Encode environmental features (one-time cost) -----
    print("\n  Loading and encoding environmental layers...")
    X_encoded, valid_mask, env_lats, env_lons = load_and_encode_env_layers(cluster_model)

    # ----- Compute soft weights (one-time cost) -----
    top_k = 1 if args.hard else args.top_k
    top_clusters, top_weights = compute_soft_weights(
        X_encoded, valid_mask, cluster_model, gmm, top_k, args.temperature
    )

    del X_encoded
    gc.collect()

    # ----- Determine dates to process -----
    if args.date:
        dates = [args.date]
    elif args.date_range:
        start = pd.Timestamp(args.date_range[0])
        end = pd.Timestamp(args.date_range[1])
        dates = pd.date_range(start, end, freq='D').strftime('%Y-%m-%d').tolist()
    else:
        print("ERROR: Specify --date or --date-range")
        sys.exit(1)

    print(f"\n  Processing {len(dates)} date(s)...")

    # ----- Process each date -----
    results = []
    for date_str in dates:
        try:
            sm_da = predict_single_date(
                date_str, models, k,
                top_clusters, top_weights, valid_mask,
                env_lats, env_lons, hard_mode=args.hard,
            )
            if sm_da is not None:
                results.append(sm_da)
        except Exception as e:
            print(f"    ERROR for {date_str}: {e}")
            continue

    if not results:
        print("No predictions generated.")
        return

    # ----- Save output -----
    if len(results) == 1:
        out_path = os.path.join(output_dir, f'wsf_sm_{dates[0]}.nc')
        results[0].to_netcdf(out_path)
    else:
        sm_combined = xr.concat(results, dim='time')
        sm_combined['time'] = pd.to_datetime(dates[:len(results)])
        out_path = os.path.join(output_dir,
                                f'wsf_sm_{dates[0]}_to_{dates[-1]}.nc')
        sm_combined.to_netcdf(out_path)

    print(f"\n  Saved: {out_path}")

    elapsed = time.time() - total_start
    print(f"\n{'=' * 70}")
    print(f"COMPLETE — {elapsed / 60:.1f} minutes")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
