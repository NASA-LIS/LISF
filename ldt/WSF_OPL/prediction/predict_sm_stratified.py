#!/usr/bin/env python
"""
Stratified Soil Moisture Prediction using Per-Cluster Random Forests
=====================================================================
Predicts soil moisture from WSF brightness temperatures using
per-cluster RF models with GMM soft blending at cluster boundaries.

Output files are LIS/LDT-compatible NetCDF:
    - Variable:  arfs_sm
    - Attribute: MAP_PROJECTION = 'EQUIDISTANT CYLINDRICAL'
    - Naming:    ARFS_SM_WSFM_YYYYMMDDTHHMMSS.nc
    - One timestep per file (time dimension length = 1)

Supports:
    1. Monthly zarr files (daily composites)
    2. Hourly NetCDF files (produces one output file per timestep)

Usage:
    # From zarr (monthly)
    python predict_sm_stratified.py zarr 202401_DES

    # From hourly NC
    python predict_sm_stratified.py hourly /path/to/wsf_hourly.nc --output-dir /path/to/out

    # Hard assignment (no blending, for comparison)
    python predict_sm_stratified.py zarr 202401_DES --hard

    # As module
    from predict_sm_stratified import StratifiedSoilMoisturePredictor
    predictor = StratifiedSoilMoisturePredictor()
    predictor.predict_from_hourly_nc('/path/to/input.nc', output_dir='/path/to/out')
"""

import numpy as np
import xarray as xr
import joblib
import json
import argparse
import os
import re
import sys
import gc
import time as timer
from pathlib import Path
from datetime import datetime
from scipy.special import softmax


class StratifiedSoilMoisturePredictor:
    """Soil moisture predictor using per-cluster RF models with soft blending."""

    # Default feature order (must match training)
    DEFAULT_FEATURES = [
        'TB_10V', 'TB_10H', 'TB_18V', 'TB_18H',
        'TB_23V', 'TB_36V', 'TB_36H', 'TB_89V', 'TB_89H'
    ]

    # Köppen grouping (must match clustering script)
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

    # Environmental layer paths
    ENV_PATHS = {
        'land_cover': '/discover/nobackup/ejalilva/data/LC/lc_arfs.nc',
        'climate':    '/discover/nobackup/ejalilva/data/LC/clm_usaf.nc',
        'clay':       '/discover/nobackup/ejalilva/data/soil/sg_clay_usaf.nc',
        'sand':       '/discover/nobackup/ejalilva/data/soil/sg_sand_usaf.nc',
        'silt':       '/discover/nobackup/ejalilva/data/soil/sg_silt_usaf.nc',
    }

    # Output variable name (must match LIS/LDT readers)
    OUTPUT_VAR_NAME = 'arfs_sm'

    # Output file prefix (ARFS_SM_WSFM_ matches the LIS/LDT glob pattern
    # ARFS_SM_*YYYYMMDDT*.nc)
    OUTPUT_FILE_PREFIX = 'ARFS_SM_WSFM'

    # Map projection string expected by LIS/LDT readers
    MAP_PROJECTION = 'EQUIDISTANT CYLINDRICAL'

    def __init__(self, model_dir=None, cluster_dir=None, config_path=None,
                 top_k=3, temperature=1.0, hard=False):
        """
        Load all cluster models and precompute soft-assignment weights.

        Parameters
        ----------
        model_dir : str, optional
            Directory containing rf_global.pkl, rf_cluster_XX.pkl, gmm_kXX.pkl
        cluster_dir : str, optional
            Directory containing environmental_clusters_kXX.nc, clustering_model_kXX.pkl
        config_path : str, optional
            Path to stratified_config.json (auto-detected from model_dir if None)
        top_k : int
            Number of top clusters to blend (default: 3)
        temperature : float
            Softmax temperature for distance-based weights (lower = sharper)
        hard : bool
            If True, use hard assignment (single cluster, no blending)
        """
        # Default paths
        if model_dir is None:
            model_dir = '/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf/cluster_models'
        if cluster_dir is None:
            cluster_dir = '/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf'
        if config_path is None:
            config_path = os.path.join(model_dir, 'stratified_config.json')

        self.model_dir = model_dir
        self.cluster_dir = cluster_dir
        self.top_k = 1 if hard else top_k
        self.temperature = temperature
        self.hard = hard

        # --- Load config ---
        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                self.config = json.load(f)
            self.features = self.config.get('features', self.DEFAULT_FEATURES)
            self.k = self.config.get('k')
        else:
            print(f"Config not found at {config_path}, using defaults")
            self.config = {}
            self.features = self.DEFAULT_FEATURES
            self.k = None

        # --- Load cluster models ---
        self.models = self._load_models()

        # --- Load cluster infrastructure + precompute weights ---
        self._setup_cluster_assignment()

        mode_str = "HARD" if self.hard else f"SOFT (top-{self.top_k})"
        print(f"Ready: K={self.k}, mode={mode_str}, features={self.features}")

    def _load_models(self):
        """Load global + per-cluster RF models."""
        print(f"Loading models from {self.model_dir}...")
        models = {}

        # Global fallback
        global_path = os.path.join(self.model_dir, 'rf_global.pkl')
        if os.path.exists(global_path):
            models['global'] = joblib.load(global_path)
            print(f"  Loaded: global model")
        else:
            raise FileNotFoundError(f"Global model not found: {global_path}")

        # Auto-detect K if not in config
        if self.k is None:
            cluster_files = sorted(Path(self.model_dir).glob('rf_cluster_*.pkl'))
            if cluster_files:
                self.k = len(cluster_files)
            else:
                raise FileNotFoundError("No cluster model files found")

        # Per-cluster
        for c in range(self.k):
            path = os.path.join(self.model_dir, f'rf_cluster_{c:02d}.pkl')
            if os.path.exists(path):
                models[c] = joblib.load(path)
            else:
                models[c] = None  # fallback to global
                print(f"  Cluster {c}: using global fallback")

        n_loaded = sum(1 for v in models.values() if v is not None)
        print(f"  Total models loaded: {n_loaded}")
        return models

    def _setup_cluster_assignment(self):
        """
        Load environmental layers, encode features, compute soft weights.
        This is a one-time cost at initialization.
        """
        print("Setting up cluster assignment...")
        t0 = timer.time()

        # Load clustering model (KMeans + encoders)
        model_path = os.path.join(self.cluster_dir, f'clustering_model_k{self.k}.pkl')
        self.cluster_model = joblib.load(model_path)
        print(f"  Loaded clustering model: {model_path}")

        # Load GMM if available
        gmm_path = os.path.join(self.model_dir, f'gmm_k{self.k}.pkl')
        if os.path.exists(gmm_path):
            self.gmm = joblib.load(gmm_path)
            print(f"  Loaded GMM: {gmm_path}")
        else:
            self.gmm = None
            print(f"  No GMM found, using distance-based soft assignment")

        # Load and encode environmental layers
        X_encoded, valid_mask, env_lats, env_lons = self._load_and_encode_env()
        self.env_lats = env_lats
        self.env_lons = env_lons
        self.env_valid_mask = valid_mask
        self.n_env_lon = len(env_lons)

        # Precompute cumulative sum for valid-index mapping
        self.env_valid_cumsum = np.cumsum(valid_mask) - 1

        # Compute soft weights for all valid pixels
        self.top_clusters, self.top_weights = self._compute_soft_weights(X_encoded, valid_mask)

        del X_encoded
        gc.collect()

        print(f"  Setup time: {timer.time() - t0:.1f}s")

    def _load_and_encode_env(self):
        """Load environmental layers and encode to clustering feature space."""
        print("  Loading environmental layers...")

        layers = {}
        for name, path in self.ENV_PATHS.items():
            if not os.path.exists(path):
                print(f"    SKIP {name}: not found")
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
            for orig, grp in self.KOPPEN_TO_GROUP.items():
                clm_grouped[clm_arr == orig] = grp
            layers['climate_group'] = xr.DataArray(
                clm_grouped, dims=layers['climate'].dims, coords=layers['climate'].coords
            )
            del layers['climate']

        # Get grid info
        ref = list(layers.values())[0]
        env_lats = ref.lat.values
        env_lons = ref.lon.values
        n_total = len(env_lats) * len(env_lons)

        # Stack raw features
        feature_names = self.cluster_model['feature_names']
        encoders = self.cluster_model['encoders']

        raw_arrays = []
        for name in feature_names:
            arr = layers[name].values
            if arr.ndim == 3:
                arr = arr[0]
            raw_arrays.append(arr.flatten())
        X_raw = np.column_stack(raw_arrays)

        # Encode (one-hot categorical, standardize continuous)
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

        print(f"    Encoded: {X_encoded.shape}, valid: {valid_mask.sum():,}/{n_total:,}")
        return X_encoded, valid_mask, env_lats, env_lons

    def _compute_soft_weights(self, X_encoded, valid_mask):
        """Compute soft cluster membership weights for all valid pixels."""
        print("  Computing soft weights...")

        X_valid = X_encoded[valid_mask]
        n_valid = X_valid.shape[0]

        if self.gmm is not None:
            # GMM posterior probabilities
            probs = self.gmm.predict_proba(X_valid)
        else:
            # Distance-based (inverse distance via softmax)
            centroids = self.cluster_model['kmeans'].cluster_centers_
            k = centroids.shape[0]
            chunk_size = 100_000
            probs = np.zeros((n_valid, k), dtype=np.float64)

            for start in range(0, n_valid, chunk_size):
                end = min(start + chunk_size, n_valid)
                X_chunk = X_valid[start:end]
                dists = np.zeros((end - start, k))
                for c in range(k):
                    dists[:, c] = np.sqrt(((X_chunk - centroids[c]) ** 2).sum(axis=1))
                probs[start:end] = softmax(-dists / self.temperature, axis=1)

        # Top-K clusters and weights
        top_clusters = np.argsort(probs, axis=1)[:, -self.top_k:][:, ::-1]
        top_weights = np.take_along_axis(probs, top_clusters, axis=1)

        # Renormalize
        weight_sums = top_weights.sum(axis=1, keepdims=True)
        weight_sums[weight_sums == 0] = 1.0
        top_weights = top_weights / weight_sums

        print(f"    Dominant weight: mean={top_weights[:, 0].mean():.3f}, "
              f"min={top_weights[:, 0].min():.3f}")

        return top_clusters, top_weights

    # =========================================================================
    # Grid mapping (WSF pixel → env pixel → cluster weights)
    # =========================================================================
    def _build_grid_mapping(self, wsf_lats, wsf_lons):
        """
        Map WSF grid pixels to environmental grid indices.
        Returns flat index array and per-pixel validity mask.
        """
        lat_idx = np.clip(np.searchsorted(self.env_lats, wsf_lats) - 1,
                          0, len(self.env_lats) - 1)
        lon_idx = np.clip(np.searchsorted(self.env_lons, wsf_lons) - 1,
                          0, len(self.env_lons) - 1)

        wsf_lat_grid, wsf_lon_grid = np.meshgrid(
            np.arange(len(wsf_lats)), np.arange(len(wsf_lons)), indexing='ij'
        )
        env_flat_idx = lat_idx[wsf_lat_grid.flatten()] * self.n_env_lon \
                       + lon_idx[wsf_lon_grid.flatten()]

        env_valid_for_wsf = self.env_valid_mask[env_flat_idx]
        return env_flat_idx, env_valid_for_wsf

    # =========================================================================
    # Core prediction (single 2D slice)
    # =========================================================================
    def _predict_2d(self, X_tb, env_flat_idx, env_valid_for_wsf):
        """
        Predict soil moisture for one 2D field (one timestep).

        Parameters
        ----------
        X_tb : np.ndarray, shape (n_pixels, n_features)
        env_flat_idx : np.ndarray, shape (n_pixels,)
        env_valid_for_wsf : np.ndarray, bool, shape (n_pixels,)

        Returns
        -------
        sm : np.ndarray, shape (n_pixels,), with NaN for invalid
        """
        n_pixels = X_tb.shape[0]
        sm = np.full(n_pixels, np.nan, dtype=np.float32)

        # Combined mask: valid TB + valid env
        tb_valid = ~np.isnan(X_tb).any(axis=1)
        predict_mask = tb_valid & env_valid_for_wsf

        if predict_mask.sum() == 0:
            return sm

        X_pred = X_tb[predict_mask]
        env_valid_idx = self.env_valid_cumsum[env_flat_idx[predict_mask]]

        # Predict with blending
        sm_valid = np.zeros(predict_mask.sum(), dtype=np.float64)

        for rank in range(self.top_k):
            cluster_ids = self.top_clusters[env_valid_idx, rank].astype(int)
            weights = self.top_weights[env_valid_idx, rank]

            preds = np.zeros(predict_mask.sum(), dtype=np.float64)
            for c in np.unique(cluster_ids):
                c_mask = cluster_ids == c
                model = self.models.get(int(c), None) or self.models['global']
                preds[c_mask] = model.predict(X_pred[c_mask])

            sm_valid += weights * preds

        sm[predict_mask] = np.clip(sm_valid, 0.0, 0.6).astype(np.float32)
        return sm

    # =========================================================================
    # Helper: parse timestamp from WSF filename
    # =========================================================================
    @staticmethod
    def _parse_time_from_filename(path):
        """
        Extract datetime from WSF filename.

        Supported patterns:
            - YYYYMMDDTHHMMSS  (e.g., 20240413T120000)
            - YYYYMMDDTHH      (e.g., 20240413T12)
            - YYYYMMDD_tHHMM   (e.g., 20240413_t0000)  — resampled WSF files
            - YYYYMMDD_tHH     (e.g., 20240413_t00)

        Parameters
        ----------
        path : str
            File path (only basename is parsed)

        Returns
        -------
        datetime
        """
        basename = os.path.basename(path)

        # Pattern 1: YYYYMMDDTHHMMSS (uppercase T, no separator)
        match = re.search(r'(\d{8})T(\d{6})', basename)
        if match:
            return datetime.strptime(match.group(1) + match.group(2),
                                     '%Y%m%d%H%M%S')

        # Pattern 2: YYYYMMDDTHH (uppercase T, 2-digit hour)
        match = re.search(r'(\d{8})T(\d{2})', basename)
        if match:
            return datetime.strptime(match.group(1) + match.group(2),
                                     '%Y%m%d%H')

        # Pattern 3: YYYYMMDD_tHHMM (underscore + lowercase t, resampled files)
        match = re.search(r'(\d{8})_t(\d{4})', basename)
        if match:
            return datetime.strptime(match.group(1) + match.group(2),
                                     '%Y%m%d%H%M')

        # Pattern 4: YYYYMMDD_tHH (underscore + lowercase t, 2-digit hour)
        match = re.search(r'(\d{8})_t(\d{2})', basename)
        if match:
            return datetime.strptime(match.group(1) + match.group(2),
                                     '%Y%m%d%H')

        raise ValueError(f"Cannot parse timestamp from filename: {basename}")

    @staticmethod
    def _parse_overpass_from_filename(path):
        """
        Extract overpass direction (ASC/DES) from filename if present.

        Returns
        -------
        str or None
            'ASC', 'DES', or None if not found
        """
        basename = os.path.basename(path).upper()
        if '_ASC' in basename or '.ASC' in basename:
            return 'ASC'
        if '_DES' in basename or '.DES' in basename:
            return 'DES'
        return None

    # =========================================================================
    # Helper: write single-timestep LIS/LDT-compatible NetCDF
    # =========================================================================
    def _write_lis_netcdf(self, sm_2d, time_val, wsf_lats, wsf_lons,
                          output_dir, overpass=None, input_ds=None):
        """
        Write one single-timestep NetCDF file compatible with LIS/LDT readers.

        Parameters
        ----------
        sm_2d : np.ndarray, shape (lat, lon)
            Soil moisture field for this timestep
        time_val : datetime or np.datetime64
            Timestamp for this field
        wsf_lats : np.ndarray
            Latitude coordinate array
        wsf_lons : np.ndarray
            Longitude coordinate array
        output_dir : str
            Directory to write the file into
        overpass : str or None
            'ASC' or 'DES' — appended to filename if provided
        input_ds : xr.Dataset, optional
            Input dataset to copy coordinate attributes from

        Returns
        -------
        output_path : str
            Path to the written file
        """
        import pandas as pd

        # Convert time to python datetime for filename formatting
        if isinstance(time_val, np.datetime64):
            dt_obj = pd.Timestamp(time_val).to_pydatetime()
        elif isinstance(time_val, datetime):
            dt_obj = time_val
        else:
            dt_obj = pd.Timestamp(time_val).to_pydatetime()

        # Build filename: ARFS_SM_WSFM_YYYYMMDDTHHMMSS_ASC.nc (or _DES)
        # The trailing _ASC/_DES is absorbed by the * in the LIS/LDT glob
        # pattern ARFS_SM_*YYYYMMDDT*.nc
        time_str = dt_obj.strftime('%Y%m%dT%H%M%S')
        if overpass:
            out_fname = f"{self.OUTPUT_FILE_PREFIX}_{time_str}_{overpass}.nc"
        else:
            out_fname = f"{self.OUTPUT_FILE_PREFIX}_{time_str}.nc"
        output_path = os.path.join(output_dir, out_fname)

        # Reshape to (1, lat, lon) — LIS/LDT readers expect time dim of length 1
        sm_3d = sm_2d.reshape(1, len(wsf_lats), len(wsf_lons))

        out_ds = xr.Dataset(
            {
                self.OUTPUT_VAR_NAME: xr.DataArray(
                    sm_3d,
                    dims=['time', 'lat', 'lon'],
                    coords={
                        'time': [time_val],
                        'lat': wsf_lats,
                        'lon': wsf_lons,
                    },
                    attrs={
                        'units': 'm3/m3',
                        'long_name': 'Volumetric Soil Moisture',
                        'valid_range': [0.0, 0.6],
                    }
                )
            },
            attrs={
                'MAP_PROJECTION': self.MAP_PROJECTION,
                'model': f'Stratified Random Forest (K={self.k})',
                'blending': 'hard' if self.hard else f'soft_top{self.top_k}',
                'features': str(self.features),
                'created': datetime.now().isoformat(),
            }
        )

        # Copy coordinate attributes from input dataset if available
        if input_ds is not None:
            for coord in ['lat', 'lon']:
                if coord in input_ds.coords and coord in out_ds.coords:
                    out_ds[coord].attrs = input_ds[coord].attrs

        out_ds.to_netcdf(output_path)
        return output_path

    # =========================================================================
    # Method 1: Predict from Monthly Zarr (Daily Composites)
    # =========================================================================
    def predict_from_zarr(self, input_name,
                          input_dir='/discover/nobackup/projects/usaf_lis/MET_FORCING/WSF/resampled/daily_composite',
                          output_dir='/discover/nobackup/projects/usaf_lis/ejalilvand/data/WSF/wsf_stratified_rf',
                          apply_qc=True, chunks=None):
        """
        Predict soil moisture from monthly zarr file.

        Parameters
        ----------
        input_name : str
            Input identifier, e.g., '202401_DES' or '202401_ASC'
        input_dir : str
            Directory containing input zarr files
        output_dir : str
            Directory for output zarr files
        apply_qc : bool
            Apply ocean mask from QUALITY_FLAG
        chunks : dict, optional
            Chunk sizes for dask

        Returns
        -------
        output_path : str
            Path to saved output zarr
        """
        input_path = f"{input_dir}/{input_name}.zarr"
        output_path = f"{output_dir}/wsf_sm_tb_only_{input_name}.zarr"

        print("=" * 60)
        print(f"ZARR PREDICTION: {input_name}")
        print("=" * 60)
        print(f"  Input:  {input_path}")
        print(f"  Output: {output_path}")

        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input not found: {input_path}")

        # Load data
        print("\nLoading data...")
        if chunks:
            ds = xr.open_zarr(input_path, chunks=chunks)
        else:
            ds = xr.open_zarr(input_path)

        print(f"  Time range: {ds.time.values[0]} to {ds.time.values[-1]}")
        print(f"  Shape: lat={len(ds.lat)}, lon={len(ds.lon)}, time={len(ds.time)}")

        # Apply QC mask
        if apply_qc and 'QUALITY_FLAG' in ds.data_vars:
            print("\nApplying quality mask (ocean)...")
            ocean_mask = (ds.QUALITY_FLAG.astype(int) & 1) > 0
            ds = ds.where(~ocean_mask)

        # Check features
        available = [v for v in self.features if v in ds.data_vars]
        missing = [v for v in self.features if v not in ds.data_vars]
        if missing:
            raise ValueError(f"Missing features in dataset: {missing}")

        # Build grid mapping (one-time for this file)
        wsf_lats = ds.lat.values
        wsf_lons = ds.lon.values
        n_lat = len(wsf_lats)
        n_lon = len(wsf_lons)
        env_flat_idx, env_valid_for_wsf = self._build_grid_mapping(wsf_lats, wsf_lons)

        # Predict all timesteps
        n_times = len(ds.time)
        sm_all = np.full((n_times, n_lat, n_lon), np.nan, dtype=np.float32)

        print(f"\nPredicting {n_times} timesteps...")
        total_valid = 0

        for t in range(n_times):
            t0 = timer.time()

            # Stack TB for this timestep
            tb_list = [ds[var].isel(time=t).values.flatten() for var in available]
            X_tb = np.column_stack(tb_list)

            # Predict
            sm_flat = self._predict_2d(X_tb, env_flat_idx, env_valid_for_wsf)
            sm_all[t] = sm_flat.reshape(n_lat, n_lon)

            n_valid = np.sum(~np.isnan(sm_all[t]))
            total_valid += n_valid

            if t % 5 == 0 or t == n_times - 1:
                print(f"  [{t+1:>3d}/{n_times}] valid={n_valid:,}, "
                      f"SM=[{np.nanmin(sm_all[t]):.3f}, {np.nanmax(sm_all[t]):.3f}], "
                      f"{timer.time()-t0:.1f}s")

        print(f"\n  Total valid predictions: {total_valid:,}")

        # Create output DataArray
        sm_da = xr.DataArray(
            sm_all,
            coords={'time': ds.time.values, 'lat': wsf_lats, 'lon': wsf_lons},
            dims=['time', 'lat', 'lon'],
            name=self.OUTPUT_VAR_NAME,
            attrs={
                'units': 'm3/m3',
                'long_name': 'Volumetric Soil Moisture',
                'valid_range': [0.0, 0.6],
                'model': f'Stratified Random Forest (K={self.k})',
                'blending': 'hard' if self.hard else f'soft_top{self.top_k}',
                'features': self.features,
                'created': datetime.now().isoformat(),
            }
        )

        out_ds = sm_da.to_dataset()
        out_ds.attrs['MAP_PROJECTION'] = self.MAP_PROJECTION

        # Copy coordinate attributes
        for coord in ['lat', 'lon', 'time']:
            if coord in ds.coords:
                out_ds[coord].attrs = ds[coord].attrs

        # Save
        print(f"\nSaving to {output_path}...")
        os.makedirs(output_dir, exist_ok=True)
        if os.path.exists(output_path):
            import shutil
            shutil.rmtree(output_path)
        out_ds.to_zarr(output_path)
        print("Done!")

        return output_path

    # =========================================================================
    # Method 2: Predict from Hourly NetCDF
    # =========================================================================
    def predict_from_hourly_nc(self, input_path, output_dir=None,
                                apply_qc=True, time_var='time'):
        """
        Predict soil moisture from hourly NetCDF file.
        Produces one LIS/LDT-compatible NetCDF file per timestep.

        Output files follow the naming convention:
            ARFS_SM_WSFM_YYYYMMDDTHHMMSS.nc

        Each file contains:
            - Variable 'arfs_sm' with dims (time=1, lat, lon)
            - Global attribute MAP_PROJECTION = 'EQUIDISTANT CYLINDRICAL'

        Parameters
        ----------
        input_path : str
            Path to input hourly NetCDF file
        output_dir : str, optional
            Directory for output files. Defaults to same directory as input.
        apply_qc : bool
            Apply ocean mask if QUALITY_FLAG exists
        time_var : str
            Name of time dimension

        Returns
        -------
        output_paths : list of str
            Paths to all written output files
        """
        print("=" * 60)
        print(f"HOURLY NC PREDICTION")
        print("=" * 60)
        print(f"  Input: {input_path}")

        if output_dir is None:
            output_dir = str(Path(input_path).parent)
        os.makedirs(output_dir, exist_ok=True)

        print(f"  Output dir: {output_dir}")

        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input not found: {input_path}")

        # Load data
        print("\nLoading data...")
        ds = xr.open_dataset(input_path, decode_timedelta=False)

        if time_var in ds.dims:
            print(f"  Time steps: {len(ds[time_var])}")
            print(f"  Time range: {ds[time_var].values[0]} to {ds[time_var].values[-1]}")

        # Apply QC mask
        # Squeeze time dim from ocean_mask — QC flag is static per file
        # but may carry a length-1 time dim that won't broadcast to 2D TB vars
        if apply_qc and 'QUALITY_FLAG' in ds.data_vars:
            print("\nApplying quality mask (ocean)...")
            ocean_mask = (ds.QUALITY_FLAG.astype(int) & 1) > 0
            if 'time' in ocean_mask.dims:
                ocean_mask = ocean_mask.squeeze('time', drop=True)
            ds = ds.where(~ocean_mask)

        # Check features
        available = [v for v in self.features if v in ds.data_vars]
        missing = [v for v in self.features if v not in ds.data_vars]
        if missing:
            raise ValueError(f"Missing features in dataset: {missing}")

        # Build grid mapping (one-time for this file)
        wsf_lats = ds.lat.values
        wsf_lons = ds.lon.values
        n_lat = len(wsf_lats)
        n_lon = len(wsf_lons)
        env_flat_idx, env_valid_for_wsf = self._build_grid_mapping(wsf_lats, wsf_lons)

        # Check if TB variables actually have a time dimension
        # (the dataset may have 'time' from other variables like 'hour',
        #  but TB fields may be 2D only)
        template = ds[available[0]]
        has_time = time_var in template.dims

        # Determine timesteps
        if has_time:
            n_times = len(ds[time_var])
        else:
            n_times = 1

        output_paths = []
        total_valid = 0

        # Extract overpass direction (ASC/DES) from input filename
        overpass = self._parse_overpass_from_filename(input_path)
        if overpass:
            print(f"  Overpass: {overpass}")

        print(f"\nPredicting {n_times} timestep(s)...")

        for t in range(n_times):
            t0 = timer.time()

            # Always parse timestamp from filename — the file's time
            # dimension only carries the date, the actual hour is
            # encoded in the filename (e.g., _t0000_ = 00:00 UTC)
            dt_obj = self._parse_time_from_filename(input_path)
            time_val = np.datetime64(dt_obj)

            # Extract TB features
            if has_time:
                tb_list = [ds[var].isel(**{time_var: t}).values.flatten()
                           for var in available]
            else:
                tb_list = [ds[var].values.flatten() for var in available]

            X_tb = np.column_stack(tb_list)

            # Predict
            sm_flat = self._predict_2d(X_tb, env_flat_idx, env_valid_for_wsf)
            sm_2d = sm_flat.reshape(n_lat, n_lon)

            n_valid = np.sum(~np.isnan(sm_2d))
            total_valid += n_valid

            # Write single-timestep LIS/LDT-compatible NetCDF
            out_path = self._write_lis_netcdf(
                sm_2d=sm_2d,
                time_val=time_val,
                wsf_lats=wsf_lats,
                wsf_lons=wsf_lons,
                output_dir=output_dir,
                overpass=overpass,
                input_ds=ds,
            )
            output_paths.append(out_path)

            elapsed = timer.time() - t0
            out_fname = os.path.basename(out_path)
            if t % 10 == 0 or t == n_times - 1:
                print(f"  [{t+1:>3d}/{n_times}] valid={n_valid:,}, "
                      f"{elapsed:.1f}s -> {out_fname}")

        print(f"\n  Total valid predictions: {total_valid:,}")
        print(f"  Files written: {len(output_paths)}")
        print("Done!")

        return output_paths


# =============================================================================
# Convenience Functions (for direct import)
# =============================================================================

_predictor = None


def get_predictor(**kwargs):
    """Get or create singleton predictor instance."""
    global _predictor
    if _predictor is None:
        _predictor = StratifiedSoilMoisturePredictor(**kwargs)
    return _predictor


def predict_zarr(input_name, **kwargs):
    """Convenience function for zarr prediction."""
    predictor = get_predictor()
    return predictor.predict_from_zarr(input_name, **kwargs)


def predict_hourly(input_path, output_dir=None, **kwargs):
    """Convenience function for hourly NC prediction."""
    predictor = get_predictor()
    return predictor.predict_from_hourly_nc(input_path, output_dir, **kwargs)


# =============================================================================
# Command Line Interface
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Predict soil moisture using stratified RF with soft blending',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Monthly zarr (daily composite)
  python predict_sm_stratified.py zarr 202504_DES
  python predict_sm_stratified.py zarr 202504_DES --hard

  # Hourly NetCDF (writes one ARFS_SM_WSFM_*.nc per timestep)
  python predict_sm_stratified.py hourly wsf_2024010112.nc
  python predict_sm_stratified.py hourly wsf_2024010112.nc --output-dir /path/to/out

  # Custom top-K blending
  python predict_sm_stratified.py zarr 202504_DES --top-k 5
        """
    )

    # Global options
    parser.add_argument('--model-dir', default=None, help='Directory with cluster RF models')
    parser.add_argument('--cluster-dir', default=None, help='Directory with clustering output')
    parser.add_argument('--top-k', type=int, default=3, help='Top-K clusters to blend (default: 3)')
    parser.add_argument('--temperature', type=float, default=1.0, help='Softmax temperature')
    parser.add_argument('--hard', action='store_true', help='Use hard assignment (no blending)')
    parser.add_argument('--no-qc', action='store_true', help='Skip QC masking')

    subparsers = parser.add_subparsers(dest='mode', help='Prediction mode')

    # --- Zarr mode ---
    zarr_parser = subparsers.add_parser('zarr', help='Predict from monthly zarr')
    zarr_parser.add_argument('input_name', help='Input name, e.g., 202401_DES')
    zarr_parser.add_argument('--input-dir',
                             default='/discover/nobackup/projects/usaf_lis/MET_FORCING/WSF/resampled/daily_composite',
                             help='Input directory')
    zarr_parser.add_argument('--output-dir',
                             default='/discover/nobackup/projects/usaf_lis/ejalilvand/data/WSF/wsf_trained_rf',
                             help='Output directory')

    # --- Hourly NC mode ---
    hourly_parser = subparsers.add_parser('hourly', help='Predict from hourly NetCDF')
    hourly_parser.add_argument('input_path', help='Input NetCDF path')
    hourly_parser.add_argument('--output-dir', '-o', default=None,
                               help='Output directory (default: same as input)')

    args = parser.parse_args()

    if args.mode is None:
        parser.print_help()
        sys.exit(1)

    # Create predictor
    predictor = StratifiedSoilMoisturePredictor(
        model_dir=args.model_dir,
        cluster_dir=args.cluster_dir,
        top_k=args.top_k,
        temperature=args.temperature,
        hard=args.hard,
    )

    if args.mode == 'zarr':
        predictor.predict_from_zarr(
            input_name=args.input_name,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            apply_qc=not args.no_qc,
        )

    elif args.mode == 'hourly':
        predictor.predict_from_hourly_nc(
            input_path=args.input_path,
            output_dir=args.output_dir,
            apply_qc=not args.no_qc,
        )


if __name__ == '__main__':
    main()
