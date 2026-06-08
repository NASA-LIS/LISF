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
    # From hourly NC (all paths required — no defaults)
    python predict_sm_stratified.py hourly /path/to/wsf_hourly.nc
        --model-dir   ./models/cluster_models
        --cluster-dir ./models
        --env-land-cover ./static_data/LC/lc_arfs.nc
        --env-climate    ./static_data/LC/clm_usaf.nc
        --env-clay       ./static_data/soil/sg_clay_usaf.nc
        --env-sand       ./static_data/soil/sg_sand_usaf.nc
        --env-silt       ./static_data/soil/sg_silt_usaf.nc
        --output-dir /path/to/out

    # From zarr (monthly)
    python predict_sm_stratified.py zarr 202401_DES
        --input-dir  ./wsf_sdr_data/resampled/daily_composite
        --output-dir ./output/zarr_sm
        [same --model-dir / --cluster-dir / --env-* flags as above]

    # Hard assignment (no blending)
    python predict_sm_stratified.py hourly ... --hard

    # As module (called from wsf_opl.py — all paths supplied by caller)
    from predict_sm_stratified import StratifiedSoilMoisturePredictor
    predictor = StratifiedSoilMoisturePredictor(
        model_dir='./models/cluster_models',
        cluster_dir='./models',
        env_paths={
            'land_cover': './static_data/LC/lc_arfs.nc',
            'climate':    './static_data/LC/clm_usaf.nc',
            'clay':       './static_data/soil/sg_clay_usaf.nc',
            'sand':       './static_data/soil/sg_sand_usaf.nc',
            'silt':       './static_data/soil/sg_silt_usaf.nc',
        },
    )
    predictor.predict_from_hourly_nc('/path/to/input.nc', output_dir='/path/to/out')
"""

# --- standard library ---
import argparse
import gc
import json
import os
import re
import sys
import time as timer
from datetime import datetime
from pathlib import Path

# --- third-party ---
import numpy as np
import pandas as pd
import xarray as xr
import joblib
from scipy.special import softmax


# X_*, X_tb, X_pred etc. are standard matrix-variable names in ML/scientific Python.
# pylint: disable=invalid-name


# pylint: disable=too-many-instance-attributes
# The predictor is a stateful object that caches models, grid mappings,
# and precomputed cluster weights — all attributes are necessary.
class StratifiedSoilMoisturePredictor:
    """Soil moisture predictor using per-cluster RF models with GMM soft blending.

    Each land pixel is assigned to one of K environmental clusters (defined by
    land cover, Köppen climate, and soil texture). A cluster-specific Random
    Forest predicts soil moisture from WSF TB channels. Soft blending smooths
    transitions at cluster boundaries by weighted-averaging the top-K models.

    Parameters
    ----------
    model_dir : str
        Directory containing rf_global.pkl, rf_cluster_XX.pkl, gmm_kXX.pkl,
        and stratified_config.json.
    cluster_dir : str
        Directory containing environmental_clusters_kXX.nc and
        clustering_model_kXX.pkl.
    env_paths : dict
        Paths to static environmental NetCDF layers keyed by:
        'land_cover', 'climate', 'clay', 'sand', 'silt'.
    config_path : str, optional
        Path to stratified_config.json. Auto-detected from model_dir if None.
    top_k : int
        Number of top clusters to blend (default: 3).
    temperature : float
        Softmax temperature for distance-based weights (lower = sharper).
    hard : bool
        If True, use hard single-cluster assignment (no blending).
    """

    # Default feature order — must match the order used during training.
    DEFAULT_FEATURES = [
        'TB_10V', 'TB_10H', 'TB_18V', 'TB_18H',
        'TB_23V', 'TB_36V', 'TB_36H', 'TB_89V', 'TB_89H',
    ]

    # Köppen climate grouping — must match the clustering script.
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

    # Output conventions — must match LIS/LDT Fortran readers.
    OUTPUT_VAR_NAME    = 'arfs_sm'
    OUTPUT_FILE_PREFIX = 'ARFS_SM_WSFM'
    MAP_PROJECTION     = 'EQUIDISTANT CYLINDRICAL'

    # pylint: disable=too-many-arguments,too-many-positional-arguments
    def __init__(self, model_dir, cluster_dir, env_paths,
                 config_path=None, top_k=3, temperature=1.0, hard=False):
        """Load all cluster models and precompute soft-assignment weights."""
        self.model_dir   = model_dir
        self.cluster_dir = cluster_dir
        self.env_paths   = env_paths
        self.top_k       = 1 if hard else top_k
        self.temperature = temperature
        self.hard        = hard

        if config_path is None:
            config_path = os.path.join(model_dir, 'stratified_config.json')

        if os.path.exists(config_path):
            with open(config_path, 'r', encoding='utf-8') as fh:
                self.config = json.load(fh)
            self.features = self.config.get('features', self.DEFAULT_FEATURES)
            self.k = self.config.get('k')
        else:
            print(f"Config not found at {config_path}, using defaults")
            self.config   = {}
            self.features = self.DEFAULT_FEATURES
            self.k        = None

        self.models = self._load_models()
        self._setup_cluster_assignment()

        mode_str = "HARD" if self.hard else f"SOFT (top-{self.top_k})"
        print(f"Ready: K={self.k}, mode={mode_str}, features={self.features}")
    # pylint: enable=too-many-arguments,too-many-positional-arguments

    # =========================================================================
    # Model loading
    # =========================================================================

    def _load_models(self):
        """Load the global fallback RF and all per-cluster RF models."""
        print(f"Loading models from {self.model_dir}...")
        models = {}

        global_path = os.path.join(self.model_dir, 'rf_global.pkl')
        if os.path.exists(global_path):
            models['global'] = joblib.load(global_path)
            print("  Loaded: global model")
        else:
            raise FileNotFoundError(f"Global model not found: {global_path}")

        if self.k is None:
            cluster_files = sorted(Path(self.model_dir).glob('rf_cluster_*.pkl'))
            if cluster_files:
                self.k = len(cluster_files)
            else:
                raise FileNotFoundError("No cluster model files found")

        for c in range(self.k):
            path = os.path.join(self.model_dir, f'rf_cluster_{c:02d}.pkl')
            if os.path.exists(path):
                models[c] = joblib.load(path)
            else:
                models[c] = None   # falls back to global at prediction time
                print(f"  Cluster {c}: using global fallback")

        n_loaded = sum(1 for v in models.values() if v is not None)
        print(f"  Total models loaded: {n_loaded}")
        return models

    # =========================================================================
    # Cluster assignment setup (one-time cost at init)
    # =========================================================================

    def _setup_cluster_assignment(self):
        """Load env layers, encode features, compute soft weights."""
        print("Setting up cluster assignment...")
        t0 = timer.time()

        model_path = os.path.join(self.cluster_dir, f'clustering_model_k{self.k}.pkl')
        self.cluster_model = joblib.load(model_path)
        print(f"  Loaded clustering model: {model_path}")

        gmm_path = os.path.join(self.model_dir, f'gmm_k{self.k}.pkl')
        if os.path.exists(gmm_path):
            self.gmm = joblib.load(gmm_path)
            print(f"  Loaded GMM: {gmm_path}")
        else:
            self.gmm = None
            print("  No GMM found, using distance-based soft assignment")

        X_encoded, valid_mask, env_lats, env_lons = self._load_and_encode_env()
        self.env_lats         = env_lats
        self.env_lons         = env_lons
        self.env_valid_mask   = valid_mask
        self.n_env_lon        = len(env_lons)
        self.env_valid_cumsum = np.cumsum(valid_mask) - 1

        self.top_clusters, self.top_weights = self._compute_soft_weights(
            X_encoded, valid_mask
        )

        del X_encoded
        gc.collect()

        print(f"  Setup time: {timer.time() - t0:.1f}s")

    # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    def _load_and_encode_env(self):
        """Load static env NetCDF layers and encode to the clustering feature space."""
        print("  Loading environmental layers...")

        layers = {}
        for name, path in self.env_paths.items():
            if not os.path.exists(path):
                print(f"    SKIP {name}: not found at {path}")
                continue
            ds = xr.open_dataset(path)
            var_name = list(ds.data_vars)[0]
            da = ds[var_name]
            if da.ndim == 3:
                da = da.squeeze()
            layers[name] = da

        if not layers:
            raise FileNotFoundError(
                "No environmental layers could be loaded. "
                "Check env_paths in pipeline_config.json."
            )

        # Remap raw Köppen codes to broad climate groups.
        if 'climate' in layers:
            clm_arr     = layers['climate'].values.copy()
            clm_grouped = np.full_like(clm_arr, np.nan, dtype=np.float64)
            for orig, grp in self.KOPPEN_TO_GROUP.items():
                clm_grouped[clm_arr == orig] = grp
            layers['climate_group'] = xr.DataArray(
                clm_grouped,
                dims=layers['climate'].dims,
                coords=layers['climate'].coords,
            )
            del layers['climate']

        ref      = list(layers.values())[0]
        env_lats = ref.lat.values
        env_lons = ref.lon.values
        n_total  = len(env_lats) * len(env_lons)

        feature_names = self.cluster_model['feature_names']
        encoders      = self.cluster_model['encoders']

        raw_arrays = []
        for name in feature_names:
            arr = layers[name].values
            if arr.ndim == 3:
                arr = arr[0]
            raw_arrays.append(arr.flatten())
        X_raw = np.column_stack(raw_arrays)

        encoded_parts = []
        for i, feat_name in enumerate(feature_names):
            col = X_raw[:, i].astype(np.float64)
            enc = encoders.get(feat_name, None)

            if enc is None:
                encoded_parts.append(col.reshape(-1, 1))
                continue

            if enc['type'] == 'onehot':
                cats   = enc['categories']
                onehot = np.zeros((n_total, len(cats)), dtype=np.float64)
                for j, cat_val in enumerate(cats):
                    onehot[:, j] = (col == cat_val).astype(np.float64)
                onehot[np.isnan(col), :] = np.nan
                encoded_parts.append(onehot)

            elif enc['type'] == 'continuous':
                scaled = np.full(n_total, np.nan, dtype=np.float64)
                valid  = ~np.isnan(col) & (col >= 0) & (col <= 100)
                if enc['std'] > 0:
                    scaled[valid] = (col[valid] - enc['mean']) / enc['std']
                encoded_parts.append(scaled.reshape(-1, 1))

        X_encoded  = np.hstack(encoded_parts)
        valid_mask = ~np.isnan(X_encoded).any(axis=1)

        print(f"    Encoded: {X_encoded.shape}, valid: {valid_mask.sum():,}/{n_total:,}")
        return X_encoded, valid_mask, env_lats, env_lons
    # pylint: enable=too-many-locals,too-many-branches,too-many-statements

    # pylint: disable=too-many-locals
    def _compute_soft_weights(self, X_encoded, valid_mask):
        """Compute soft cluster membership weights for all valid env pixels."""
        print("  Computing soft weights...")

        X_valid = X_encoded[valid_mask]
        n_valid = X_valid.shape[0]

        if self.gmm is not None:
            probs = self.gmm.predict_proba(X_valid)
        else:
            centroids  = self.cluster_model['kmeans'].cluster_centers_
            k          = centroids.shape[0]
            chunk_size = 100_000
            probs      = np.zeros((n_valid, k), dtype=np.float64)

            for start in range(0, n_valid, chunk_size):
                end     = min(start + chunk_size, n_valid)
                X_chunk = X_valid[start:end]
                dists   = np.zeros((end - start, k))
                for c in range(k):
                    dists[:, c] = np.sqrt(((X_chunk - centroids[c]) ** 2).sum(axis=1))
                probs[start:end] = softmax(-dists / self.temperature, axis=1)

        top_clusters = np.argsort(probs, axis=1)[:, -self.top_k:][:, ::-1]
        top_weights  = np.take_along_axis(probs, top_clusters, axis=1)

        weight_sums = top_weights.sum(axis=1, keepdims=True)
        weight_sums[weight_sums == 0] = 1.0
        top_weights = top_weights / weight_sums

        print(f"    Dominant weight: mean={top_weights[:, 0].mean():.3f}, "
              f"min={top_weights[:, 0].min():.3f}")

        return top_clusters, top_weights

    # =========================================================================
    # Grid mapping
    # =========================================================================

    def _build_grid_mapping(self, wsf_lats, wsf_lons):
        """Map WSF grid pixels to environmental grid flat indices."""
        lat_idx = np.clip(
            np.searchsorted(self.env_lats, wsf_lats) - 1,
            0, len(self.env_lats) - 1,
        )
        lon_idx = np.clip(
            np.searchsorted(self.env_lons, wsf_lons) - 1,
            0, len(self.env_lons) - 1,
        )

        wsf_lat_grid, wsf_lon_grid = np.meshgrid(
            np.arange(len(wsf_lats)), np.arange(len(wsf_lons)), indexing='ij'
        )
        env_flat_idx = (
            lat_idx[wsf_lat_grid.flatten()] * self.n_env_lon
            + lon_idx[wsf_lon_grid.flatten()]
        )

        env_valid_for_wsf = self.env_valid_mask[env_flat_idx]
        return env_flat_idx, env_valid_for_wsf

    # =========================================================================
    # Core prediction
    # =========================================================================

    # pylint: disable=too-many-locals
    def _predict_2d(self, X_tb, env_flat_idx, env_valid_for_wsf):
        """Predict soil moisture for one 2D timestep using soft-blended RF models.

        Parameters
        ----------
        X_tb : np.ndarray, shape (n_pixels, n_tb_features)
        env_flat_idx : np.ndarray, shape (n_pixels,)
        env_valid_for_wsf : np.ndarray, bool, shape (n_pixels,)

        Returns
        -------
        sm : np.ndarray, shape (n_pixels,), NaN for invalid pixels
        """
        n_pixels = X_tb.shape[0]
        sm = np.full(n_pixels, np.nan, dtype=np.float32)

        tb_valid     = ~np.isnan(X_tb).any(axis=1)
        predict_mask = tb_valid & env_valid_for_wsf

        if predict_mask.sum() == 0:
            return sm

        X_pred        = X_tb[predict_mask]
        env_valid_idx = self.env_valid_cumsum[env_flat_idx[predict_mask]]

        sm_valid = np.zeros(predict_mask.sum(), dtype=np.float64)

        for rank in range(self.top_k):
            cluster_ids = self.top_clusters[env_valid_idx, rank].astype(int)
            weights     = self.top_weights[env_valid_idx, rank]

            preds = np.zeros(predict_mask.sum(), dtype=np.float64)
            for c in np.unique(cluster_ids):
                c_mask = cluster_ids == c
                model  = self.models.get(int(c), None) or self.models['global']
                preds[c_mask] = model.predict(X_pred[c_mask])

            sm_valid += weights * preds

        sm[predict_mask] = np.clip(sm_valid, 0.0, 0.6).astype(np.float32)
        return sm
    # pylint: enable=too-many-locals

    # =========================================================================
    # Filename helpers
    # =========================================================================

    @staticmethod
    def _parse_time_from_filename(path):
        """Extract datetime from a WSF filename.

        Supported patterns:
            YYYYMMDDTHHMMSS  e.g. 20240413T120000
            YYYYMMDDTHH      e.g. 20240413T12
            YYYYMMDD_tHHMM   e.g. 20240413_t0000  (resampled files)
            YYYYMMDD_tHH     e.g. 20240413_t00
        """
        basename = os.path.basename(path)

        match = re.search(r'(\d{8})T(\d{6})', basename)
        if match:
            return datetime.strptime(match.group(1) + match.group(2), '%Y%m%d%H%M%S')

        match = re.search(r'(\d{8})T(\d{2})', basename)
        if match:
            return datetime.strptime(match.group(1) + match.group(2), '%Y%m%d%H')

        match = re.search(r'(\d{8})_t(\d{4})', basename)
        if match:
            return datetime.strptime(match.group(1) + match.group(2), '%Y%m%d%H%M')

        match = re.search(r'(\d{8})_t(\d{2})', basename)
        if match:
            return datetime.strptime(match.group(1) + match.group(2), '%Y%m%d%H')

        raise ValueError(f"Cannot parse timestamp from filename: {basename}")

    @staticmethod
    def _parse_overpass_from_filename(path):
        """Return 'ASC', 'DES', or None based on the filename."""
        basename = os.path.basename(path).upper()
        if '_ASC' in basename or '.ASC' in basename:
            return 'ASC'
        if '_DES' in basename or '.DES' in basename:
            return 'DES'
        return None

    # =========================================================================
    # NetCDF writer
    # =========================================================================

    # pylint: disable=too-many-arguments,too-many-positional-arguments
    def _write_lis_netcdf(self, sm_2d, time_val, wsf_lats, wsf_lons,
                          output_dir, overpass=None, input_ds=None):
        """Write one single-timestep LIS/LDT-compatible NetCDF file.

        Parameters
        ----------
        sm_2d : np.ndarray, shape (lat, lon)
        time_val : datetime or np.datetime64
        wsf_lats : np.ndarray
        wsf_lons : np.ndarray
        output_dir : str
        overpass : str or None — 'ASC' or 'DES'
        input_ds : xr.Dataset, optional — source for coordinate attributes

        Returns
        -------
        output_path : str
        """
        if isinstance(time_val, np.datetime64):
            dt_obj = pd.Timestamp(time_val).to_pydatetime()
        elif isinstance(time_val, datetime):
            dt_obj = time_val
        else:
            dt_obj = pd.Timestamp(time_val).to_pydatetime()

        time_str = dt_obj.strftime('%Y%m%dT%H%M%S')
        out_fname = (
            f"{self.OUTPUT_FILE_PREFIX}_{time_str}_{overpass}.nc"
            if overpass
            else f"{self.OUTPUT_FILE_PREFIX}_{time_str}.nc"
        )
        output_path = os.path.join(output_dir, out_fname)

        sm_3d  = sm_2d.reshape(1, len(wsf_lats), len(wsf_lons))
        out_ds = xr.Dataset(
            {
                self.OUTPUT_VAR_NAME: xr.DataArray(
                    sm_3d,
                    dims=['time', 'lat', 'lon'],
                    coords={'time': [time_val], 'lat': wsf_lats, 'lon': wsf_lons},
                    attrs={
                        'units':       'm3/m3',
                        'long_name':   'Volumetric Soil Moisture',
                        'valid_range': [0.0, 0.6],
                    },
                )
            },
            attrs={
                'MAP_PROJECTION': self.MAP_PROJECTION,
                'model':          f'Stratified Random Forest (K={self.k})',
                'blending':       'hard' if self.hard else f'soft_top{self.top_k}',
                'features':       str(self.features),
                'created':        datetime.now().isoformat(),
            },
        )

        if input_ds is not None:
            for coord in ['lat', 'lon']:
                if coord in input_ds.coords and coord in out_ds.coords:
                    out_ds[coord].attrs = input_ds[coord].attrs

        out_ds.to_netcdf(output_path)
        return output_path
    # pylint: enable=too-many-arguments,too-many-positional-arguments

    # =========================================================================
    # Public prediction methods
    # =========================================================================

    # pylint: disable=too-many-locals,too-many-arguments,too-many-positional-arguments
    def predict_from_zarr(self, input_name, input_dir, output_dir,
                          apply_qc=True, chunks=None):
        """Predict soil moisture from a monthly zarr composite file.

        Parameters
        ----------
        input_name : str — e.g. '202401_DES'
        input_dir : str — directory containing input zarr files
        output_dir : str — directory for output zarr files
        apply_qc : bool — apply ocean mask from QUALITY_FLAG
        chunks : dict, optional — dask chunk sizes
        """
        input_path  = f"{input_dir}/{input_name}.zarr"
        output_path = f"{output_dir}/wsf_sm_tb_only_{input_name}.zarr"

        print("=" * 60)
        print(f"ZARR PREDICTION: {input_name}")
        print(f"  Input:  {input_path}")
        print(f"  Output: {output_path}")

        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input not found: {input_path}")

        ds = xr.open_zarr(input_path, chunks=chunks) if chunks else xr.open_zarr(input_path)

        if apply_qc and 'QUALITY_FLAG' in ds.data_vars:
            ocean_mask = (ds.QUALITY_FLAG.astype(int) & 1) > 0
            ds = ds.where(~ocean_mask)

        available = [v for v in self.features if v in ds.data_vars]
        missing   = [v for v in self.features if v not in ds.data_vars]
        if missing:
            raise ValueError(f"Missing features in dataset: {missing}")

        wsf_lats = ds.lat.values
        wsf_lons = ds.lon.values
        n_lat    = len(wsf_lats)
        n_lon    = len(wsf_lons)
        n_times  = len(ds.time)

        env_flat_idx, env_valid_for_wsf = self._build_grid_mapping(wsf_lats, wsf_lons)
        sm_all = np.full((n_times, n_lat, n_lon), np.nan, dtype=np.float32)

        print(f"\nPredicting {n_times} timesteps...")
        for t in range(n_times):
            t0     = timer.time()
            tb_list = [ds[var].isel(time=t).values.flatten() for var in available]
            X_tb   = np.column_stack(tb_list)
            sm_flat = self._predict_2d(X_tb, env_flat_idx, env_valid_for_wsf)
            sm_all[t] = sm_flat.reshape(n_lat, n_lon)

            if t % 5 == 0 or t == n_times - 1:
                n_valid = np.sum(~np.isnan(sm_all[t]))
                print(f"  [{t+1:>3d}/{n_times}] valid={n_valid:,}, "
                      f"SM=[{np.nanmin(sm_all[t]):.3f}, {np.nanmax(sm_all[t]):.3f}], "
                      f"{timer.time()-t0:.1f}s")

        sm_da = xr.DataArray(
            sm_all,
            coords={'time': ds.time.values, 'lat': wsf_lats, 'lon': wsf_lons},
            dims=['time', 'lat', 'lon'],
            name=self.OUTPUT_VAR_NAME,
            attrs={
                'units':       'm3/m3',
                'long_name':   'Volumetric Soil Moisture',
                'valid_range': [0.0, 0.6],
                'model':       f'Stratified Random Forest (K={self.k})',
                'blending':    'hard' if self.hard else f'soft_top{self.top_k}',
                'features':    self.features,
                'created':     datetime.now().isoformat(),
            },
        )
        out_ds = sm_da.to_dataset()
        out_ds.attrs['MAP_PROJECTION'] = self.MAP_PROJECTION
        for coord in ['lat', 'lon', 'time']:
            if coord in ds.coords:
                out_ds[coord].attrs = ds[coord].attrs

        os.makedirs(output_dir, exist_ok=True)
        if os.path.exists(output_path):
            import shutil  # pylint: disable=import-outside-toplevel
            shutil.rmtree(output_path)
        out_ds.to_zarr(output_path)
        print(f"\n  Saved: {output_path}")
        return output_path
    # pylint: enable=too-many-locals,too-many-arguments,too-many-positional-arguments

    # pylint: disable=too-many-locals
    def predict_from_hourly_nc(self, input_path, output_dir=None,
                               apply_qc=True, time_var='time'):
        """Predict soil moisture from a single hourly resampled NetCDF file.

        Produces one LIS/LDT-compatible NetCDF per timestep named
        ARFS_SM_WSFM_YYYYMMDDTHHMMSS[_ASC|_DES].nc.

        Parameters
        ----------
        input_path : str
        output_dir : str, optional — defaults to same directory as input
        apply_qc : bool — apply ocean mask if QUALITY_FLAG present
        time_var : str — name of the time dimension

        Returns
        -------
        output_paths : list of str
        """
        dt_obj = self._parse_time_from_filename(input_path)
        print(f"Predicting: {dt_obj.strftime('%Y-%m-%d %H:%M UTC')}")

        if output_dir is None:
            output_dir = str(Path(input_path).parent)
        os.makedirs(output_dir, exist_ok=True)

        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input not found: {input_path}")

        ds = xr.open_dataset(input_path, decode_timedelta=False)

        if apply_qc and 'QUALITY_FLAG' in ds.data_vars:
            ocean_mask = (ds.QUALITY_FLAG.astype(int) & 1) > 0
            if 'time' in ocean_mask.dims:
                ocean_mask = ocean_mask.squeeze('time', drop=True)
            ds = ds.where(~ocean_mask)

        available = [v for v in self.features if v in ds.data_vars]
        missing   = [v for v in self.features if v not in ds.data_vars]
        if missing:
            raise ValueError(f"Missing features in dataset: {missing}")

        wsf_lats = ds.lat.values
        wsf_lons = ds.lon.values
        n_lat    = len(wsf_lats)
        n_lon    = len(wsf_lons)

        env_flat_idx, env_valid_for_wsf = self._build_grid_mapping(wsf_lats, wsf_lons)

        has_time = time_var in ds[available[0]].dims
        n_times  = len(ds[time_var]) if has_time else 1
        overpass = self._parse_overpass_from_filename(input_path)

        output_paths = []

        for t in range(n_times):
            t0       = timer.time()
            time_val = np.datetime64(self._parse_time_from_filename(input_path))

            tb_list = (
                [ds[var].isel(**{time_var: t}).values.flatten() for var in available]
                if has_time
                else [ds[var].values.flatten() for var in available]
            )
            X_tb    = np.column_stack(tb_list)
            sm_flat = self._predict_2d(X_tb, env_flat_idx, env_valid_for_wsf)
            sm_2d   = sm_flat.reshape(n_lat, n_lon)

            out_path = self._write_lis_netcdf(
                sm_2d=sm_2d, time_val=time_val,
                wsf_lats=wsf_lats, wsf_lons=wsf_lons,
                output_dir=output_dir, overpass=overpass, input_ds=ds,
            )
            output_paths.append(out_path)

            if t % 10 == 0 or t == n_times - 1:
                n_valid = np.sum(~np.isnan(sm_2d))
                print(f"  [{t+1:>3d}/{n_times}] valid={n_valid:,}, "
                      f"{timer.time()-t0:.1f}s -> {os.path.basename(out_path)}")

        print(f"  Files written: {len(output_paths)}")
        return output_paths
    # pylint: enable=too-many-locals


# =============================================================================
# Command Line Interface
# =============================================================================

def main():
    """Parse CLI arguments and run prediction."""
    parser = argparse.ArgumentParser(
        description='Predict soil moisture using stratified RF with soft blending',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required path arguments — no defaults, no hardcoded paths.
    parser.add_argument('--model-dir',      required=True, help='Directory with cluster RF models')
    parser.add_argument('--cluster-dir',    required=True, help='Directory with clustering output')
    parser.add_argument('--env-land-cover', required=True, help='Land cover .nc path')
    parser.add_argument('--env-climate',    required=True, help='Climate .nc path')
    parser.add_argument('--env-clay',       required=True, help='Clay fraction .nc path')
    parser.add_argument('--env-sand',       required=True, help='Sand fraction .nc path')
    parser.add_argument('--env-silt',       required=True, help='Silt fraction .nc path')

    # Optional arguments.
    parser.add_argument('--top-k',       type=int,   default=3,   help='Top-K clusters to blend')
    parser.add_argument('--temperature', type=float, default=1.0, help='Softmax temperature')
    parser.add_argument('--hard',        action='store_true',     help='Use hard assignment')
    parser.add_argument('--no-qc',       action='store_true',     help='Skip QC masking')

    subparsers = parser.add_subparsers(dest='mode', help='Prediction mode')

    zarr_parser = subparsers.add_parser('zarr', help='Predict from monthly zarr')
    zarr_parser.add_argument('input_name')
    zarr_parser.add_argument('--input-dir',  required=True, help='Input zarr directory')
    zarr_parser.add_argument('--output-dir', required=True, help='Output directory')

    hourly_parser = subparsers.add_parser('hourly', help='Predict from hourly NetCDF')
    hourly_parser.add_argument('input_path')
    hourly_parser.add_argument('--output-dir', '-o', required=True, help='Output directory')

    args = parser.parse_args()

    if args.mode is None:
        parser.print_help()
        sys.exit(1)

    env_paths = {
        'land_cover': args.env_land_cover,
        'climate':    args.env_climate,
        'clay':       args.env_clay,
        'sand':       args.env_sand,
        'silt':       args.env_silt,
    }

    predictor = StratifiedSoilMoisturePredictor(
        model_dir=args.model_dir,
        cluster_dir=args.cluster_dir,
        env_paths=env_paths,
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
