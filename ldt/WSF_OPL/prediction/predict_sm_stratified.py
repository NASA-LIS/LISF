#!/usr/bin/env python

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.8
#
# Copyright (c) 2026 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
Stratified Soil Moisture Prediction using Per-Cluster Random Forests
=====================================================================
Predicts soil moisture from WSF brightness temperatures using
per-cluster RF models, blended across clusters at their boundaries.

Each pixel gets weights w_k for its top-K environmental clusters, and its
soil moisture is the weighted sum of those clusters' Random Forests:

    SM = sum_k w_k * RF_k(TB)        (k over the top-K clusters, clipped 0-0.6)

The weights are chosen with --blending (required, no default):
    distance — from the distance between the pixel's encoded static features
               x (land-cover and climate flags, standardised clay/sand/silt)
               and each KMeans cluster centre c_k:
                   d_k = ||x - c_k||                       (Euclidean)
                   w_k = exp(-d_k / T) / sum_j exp(-d_j / T)
               with j over the same top-K (nearest) clusters and T the
               --temperature (required). Small T: the nearest cluster
               dominates and only pixels between clusters blend (T -> 0 is
               'hard'); large T: the weights approach 1/K everywhere.
               Example, d = 1.0 and 1.5: T = 1 gives 0.62/0.38,
               T = 0.2 gives 0.92/0.08, T = 0.1 gives 0.99/0.01.
    hard     — the nearest cluster centre only (weight 1): the cluster the
               pixel was trained in.
    gmm      — posterior probabilities of a Gaussian mixture fitted to the
               encoded static layers (gmm_kXX.pkl in the model directory).
               Its most probable component can differ from the cluster the
               pixel was trained in.

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
        --blending distance --temperature 0.1
        --output-dir /path/to/out

    # From zarr (monthly)
    python predict_sm_stratified.py zarr 202401_DES
        --input-dir  ./wsf_sdr_data/resampled/daily_composite
        --output-dir ./output/zarr_sm
        [same --model-dir / --cluster-dir / --env-* flags as above]

    # Distance-based softmax, or hard assignment (no blending)
    python predict_sm_stratified.py hourly ... --blending distance --temperature 0.1
    python predict_sm_stratified.py hourly ... --blending hard

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
        blending='distance',
        temperature=0.1,
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

# How the cluster weights of each pixel are computed (see module docstring).
BLENDING_MODES = ('gmm', 'distance', 'hard')

# Static layers, keyed as in env_paths.
ENV_LAYER_NAMES = ('land_cover', 'climate', 'clay', 'sand', 'silt')


def load_env_layers(env_paths):
    """Read the five static layers and group the Köppen classes.

    Every layer must exist, have dims (lat, lon) after dropping length-1
    dims, and share the same lat/lon. Returns a dict of 2-D DataArrays
    keyed by clustering feature name ('climate' becomes 'climate_group').
    """
    missing = [name for name in ENV_LAYER_NAMES if name not in env_paths]
    if missing:
        raise KeyError(f"env_paths has no entry for: {', '.join(missing)}")

    layers = {}
    for name in ENV_LAYER_NAMES:
        path = env_paths[name]
        if not os.path.exists(path):
            raise FileNotFoundError(f"Static layer '{name}' not found: {path}")
        with xr.open_dataset(path) as ds:
            da = ds[list(ds.data_vars)[0]].squeeze(drop=True).load()
        if da.dims != ('lat', 'lon'):
            raise ValueError(f"Static layer '{name}' must have dims (lat, lon), "
                             f"found {da.dims} in {path}")
        layers[name] = da

    ref = layers['land_cover']
    for name, da in layers.items():
        if (da.shape != ref.shape
                or not np.allclose(da['lat'].values, ref['lat'].values, atol=1e-4)
                or not np.allclose(da['lon'].values, ref['lon'].values, atol=1e-4)):
            raise ValueError(f"Static layer '{name}' ({env_paths[name]}) is not on "
                             "the same lat/lon grid as 'land_cover'")

    climate = layers.pop('climate')
    grouped = np.full(climate.shape, np.nan, dtype=np.float64)
    codes = climate.values
    for orig, grp in KOPPEN_TO_GROUP.items():
        grouped[codes == orig] = grp
    layers['climate_group'] = climate.copy(data=grouped)
    return layers


def encode_env_features(layers, cluster_model):
    """Encode the static layers into the clustering feature space.

    Uses the encoders stored with the KMeans model. A cell is invalid
    (NaN row) where any layer is missing, where a soil fraction is outside
    0-100, or where its class was not used in the clustering (water,
    snow/ice and polar climate were excluded there, so they are excluded
    here too). Rows are in (lat, lon) row-major order.

    Returns (X_encoded, valid_mask).
    """
    parts = []
    for name in cluster_model['feature_names']:
        if name not in layers:
            raise KeyError(f"Clustering feature '{name}' has no static layer")
        col = np.asarray(layers[name].values, dtype=np.float64).ravel()
        enc = cluster_model['encoders'].get(name)
        if enc is None:
            parts.append(col.reshape(-1, 1))
        elif enc['type'] == 'onehot':
            cats = np.asarray(enc['categories'], dtype=np.float64)
            onehot = (col[:, None] == cats[None, :]).astype(np.float64)
            onehot[~onehot.any(axis=1), :] = np.nan   # missing or unused class
            parts.append(onehot)
        elif enc['type'] == 'continuous':
            scaled = np.full(col.shape, np.nan, dtype=np.float64)
            ok = np.isfinite(col) & (col >= 0) & (col <= 100)
            if enc['std'] > 0:
                scaled[ok] = (col[ok] - enc['mean']) / enc['std']
            parts.append(scaled.reshape(-1, 1))
        else:
            raise ValueError(f"Unknown encoder type {enc['type']!r} for '{name}'")
    X_encoded = np.hstack(parts)
    return X_encoded, ~np.isnan(X_encoded).any(axis=1)


# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
def cluster_weights(X_valid, cluster_model, gmm, blending, top_k=3,
                    temperature=1.0, chunk_size=100_000):
    """Top-k cluster indices and their normalised weights for each cell.

    blending is one of BLENDING_MODES; gmm is the fitted GaussianMixture
    (used only for 'gmm'). 'hard' returns the nearest KMeans centroid with
    weight 1, i.e. the label KMeans gives the cell. 'distance' uses
        d_k = ||x - c_k||,   w_k = exp(-d_k / T) / sum_j exp(-d_j / T)
    over the top_k nearest centroids c_k (softmax of -d/T over all clusters,
    then the top_k weights renormalised to sum to 1, which is the same).
    """
    if blending not in BLENDING_MODES:
        raise ValueError(f"blending must be one of {BLENDING_MODES}, got {blending!r}")
    if blending == 'gmm' and gmm is None:
        raise ValueError("blending 'gmm' needs a fitted GMM")
    if blending == 'distance' and temperature <= 0:
        raise ValueError(f"temperature must be > 0, got {temperature}")

    centroids = cluster_model['kmeans'].cluster_centers_
    k = centroids.shape[0]
    top_k = 1 if blending == 'hard' else min(top_k, k)
    n_valid = X_valid.shape[0]
    top_clusters = np.zeros((n_valid, top_k), dtype=np.int64)
    top_weights = np.zeros((n_valid, top_k), dtype=np.float64)

    for start in range(0, n_valid, chunk_size):
        end = min(start + chunk_size, n_valid)
        X_chunk = X_valid[start:end]
        if blending == 'gmm':
            probs = gmm.predict_proba(X_chunk)
        else:
            dists = np.zeros((end - start, k))
            for c in range(k):
                dists[:, c] = np.sqrt(((X_chunk - centroids[c]) ** 2).sum(axis=1))
            # softmax(-d / T): w_k = exp(-d_k / T) / sum_j exp(-d_j / T)
            probs = softmax(-dists / (temperature if blending == 'distance' else 1.0),
                            axis=1)
        idx = np.argsort(probs, axis=1)[:, ::-1][:, :top_k]
        weights = np.take_along_axis(probs, idx, axis=1)
        sums = weights.sum(axis=1, keepdims=True)
        sums[sums == 0] = 1.0
        top_clusters[start:end] = idx
        top_weights[start:end] = weights / sums
    return top_clusters, top_weights
# pylint: enable=too-many-arguments,too-many-positional-arguments,too-many-locals


# pylint: disable=too-many-instance-attributes
# The predictor is a stateful object that caches models, grid mappings,
# and precomputed cluster weights — all attributes are necessary.
class StratifiedSoilMoisturePredictor:
    """Soil moisture predictor using per-cluster RF models with cluster blending.

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
        Directory containing clustering_model_kXX.pkl (KMeans centres and
        the encoders of the static layers). The cluster map
        environmental_clusters_kXX.nc may be there too, but it is not read:
        the cluster weights and the valid-cell mask are computed from the
        static layers (env_paths) when the predictor is set up.
    env_paths : dict
        Paths to static environmental NetCDF layers keyed by:
        'land_cover', 'climate', 'clay', 'sand', 'silt'.
    config_path : str, optional
        Path to stratified_config.json. Auto-detected from model_dir if None.
    blending : str
        'gmm', 'distance' or 'hard' (required; see module docstring).
    top_k : int
        Number of top clusters to blend (default: 3; 1 for 'hard').
    temperature : float
        T in w_k = exp(-d_k / T) / sum_j exp(-d_j / T) (required for
        'distance'; lower = sharper, e.g. 0.1). See the module docstring.
    """

    # Output conventions — must match LIS/LDT Fortran readers.
    OUTPUT_VAR_NAME    = 'arfs_sm'
    OUTPUT_FILE_PREFIX = 'ARFS_SM_WSFM'
    MAP_PROJECTION     = 'EQUIDISTANT CYLINDRICAL'

    # pylint: disable=too-many-arguments
    def __init__(self, model_dir, cluster_dir, env_paths, *, blending,
                 config_path=None, top_k=3, temperature=None):
        """Load all cluster models and precompute the cluster weights."""
        if blending not in BLENDING_MODES:
            raise ValueError(f"blending must be one of {BLENDING_MODES}, got {blending!r}")
        if blending == 'distance' and (temperature is None or temperature <= 0):
            raise ValueError("blending 'distance' needs a temperature > 0 "
                             f"(got {temperature!r}), e.g. 0.1")
        self.model_dir   = model_dir
        self.cluster_dir = cluster_dir
        self.env_paths   = env_paths
        self.blending    = blending
        self.top_k       = 1 if blending == 'hard' else top_k
        self.temperature = temperature

        if config_path is None:
            config_path = os.path.join(model_dir, 'stratified_config.json')
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Model config not found: {config_path}")
        with open(config_path, 'r', encoding='utf-8') as fh:
            self.config = json.load(fh)
        for key in ('features', 'k'):
            if key not in self.config:
                raise KeyError(f"'{key}' missing from {config_path}")
        self.features = self.config['features']
        self.k = self.config['k']

        self.models = self._load_models()
        self._setup_cluster_assignment()

        mode_str = self.blending if self.blending == 'hard' else \
            f"{self.blending} (top-{self.top_k})"
        print(f"Ready: K={self.k}, blending={mode_str}, features={self.features}")
    # pylint: enable=too-many-arguments

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
        """Load the clustering model (and GMM), encode env layers, compute weights."""
        print("Setting up cluster assignment...")
        t0 = timer.time()

        model_path = os.path.join(self.cluster_dir, f'clustering_model_k{self.k}.pkl')
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Clustering model not found: {model_path}")
        self.cluster_model = joblib.load(model_path)
        print(f"  Loaded clustering model: {model_path}")

        self.gmm = None
        if self.blending == 'gmm':
            gmm_path = os.path.join(self.model_dir, f'gmm_k{self.k}.pkl')
            if not os.path.exists(gmm_path):
                raise FileNotFoundError(
                    f"blending 'gmm' needs {gmm_path}; fit it with "
                    "build_cluster_map.py, or use blending 'distance' or 'hard'")
            self.gmm = joblib.load(gmm_path)
            print(f"  Loaded GMM: {gmm_path}")

        print("  Loading environmental layers...")
        layers = load_env_layers(self.env_paths)
        X_encoded, valid_mask = encode_env_features(layers, self.cluster_model)
        print(f"    Encoded: {X_encoded.shape}, valid: {valid_mask.sum():,}/{valid_mask.size:,}")

        self.env_lats         = layers['land_cover']['lat'].values
        self.env_lons         = layers['land_cover']['lon'].values
        self.env_valid_mask   = valid_mask
        self.n_env_lon        = len(self.env_lons)
        self.env_valid_cumsum = np.cumsum(valid_mask) - 1

        print(f"  Computing cluster weights ({self.blending})...")
        self.top_clusters, self.top_weights = cluster_weights(
            X_encoded[valid_mask], self.cluster_model, self.gmm, self.blending,
            top_k=self.top_k, temperature=self.temperature,
        )
        print(f"    Dominant weight: mean={self.top_weights[:, 0].mean():.3f}, "
              f"min={self.top_weights[:, 0].min():.3f}")

        del X_encoded, layers
        gc.collect()
        print(f"  Setup time: {timer.time() - t0:.1f}s")

    # =========================================================================
    # Grid mapping
    # =========================================================================

    def _blending_label(self):
        """Blending description written to the output attributes."""
        if self.blending == 'hard':
            return 'hard'
        label = f'{self.blending}_top{self.top_k}'
        if self.blending == 'distance':
            label += f'_T{self.temperature:g}'
        return label

    @staticmethod
    def _nearest_index(env_coord, target):
        """Index of the nearest env coordinate for each target coordinate.

        Works whether env_coord is ascending or descending, and for any
        WSF grid (the WSF grid need not match the env grid).
        """
        env = np.asarray(env_coord, dtype=np.float64)
        tgt = np.asarray(target, dtype=np.float64)
        order = np.argsort(env)
        env_sorted = env[order]
        if env_sorted.size == 1:
            return np.zeros(tgt.shape, dtype=np.int64)
        pos = np.clip(np.searchsorted(env_sorted, tgt), 1, env_sorted.size - 1)
        closer_left = (tgt - env_sorted[pos - 1]) <= (env_sorted[pos] - tgt)
        return order[pos - closer_left.astype(np.int64)]

    def _build_grid_mapping(self, wsf_lats, wsf_lons):
        """Map WSF grid pixels to environmental grid flat indices (nearest cell)."""
        lat_idx = self._nearest_index(self.env_lats, wsf_lats)
        lon_idx = self._nearest_index(self.env_lons, wsf_lons)

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
                'blending':       self._blending_label(),
                'features':       str(self.features),
                'created':        datetime.now().isoformat(),
            },
        )

        if input_ds is not None:
            for coord in ['lat', 'lon']:
                if coord in input_ds.coords and coord in out_ds.coords:
                    out_ds[coord].attrs = input_ds[coord].attrs

        # Compressed: most cells of an hourly swath are NaN, so this cuts the
        # file from ~26 MB to a small fraction; readers see no difference.
        out_ds.to_netcdf(output_path, encoding={
            self.OUTPUT_VAR_NAME: {'zlib': True, 'complevel': 4}})
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
            # No-data cells (fill value, read as NaN) count as flag 0: NaN
            # cannot be cast to an integer ("invalid value encountered in cast").
            ocean_mask = (ds.QUALITY_FLAG.fillna(0).astype(int) & 1) > 0
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
                'blending':    self._blending_label(),
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
            # No-data cells (fill value, read as NaN) count as flag 0: NaN
            # cannot be cast to an integer ("invalid value encountered in cast").
            ocean_mask = (ds.QUALITY_FLAG.fillna(0).astype(int) & 1) > 0
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
        description='Predict soil moisture using stratified RF with cluster blending',
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
    parser.add_argument('--blending',       required=True, choices=BLENDING_MODES,
                        help='How the top-K cluster RFs are combined: distance '
                             '(weights from the distance to the cluster centres, '
                             'see --temperature), hard (own cluster only) or gmm')

    # Optional arguments.
    parser.add_argument('--top-k',       type=int,   default=3,   help='Top-K clusters to blend')
    parser.add_argument('--temperature', type=float, default=None,
                        help='T in w_k = exp(-d_k/T) / sum_j exp(-d_j/T), d_k = distance '
                             'to cluster centre k (required with --blending distance; '
                             'smaller = sharper, e.g. 0.1)')
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
    if args.blending == 'distance' and args.temperature is None:
        parser.error('--blending distance needs --temperature (e.g. 0.1)')

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
        blending=args.blending,
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
