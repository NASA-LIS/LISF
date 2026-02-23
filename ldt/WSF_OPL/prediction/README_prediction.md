# WSF Soil Moisture Prediction — Stratified Random Forest

Predicts volumetric soil moisture (m³/m³) from WSF brightness temperatures using **per-cluster Random Forest models** (K=30). Supports hard assignment from a pre-computed cluster grid or GMM soft blending at cluster boundaries. Output is LIS/LDT-compatible NetCDF.

---

## How It Works

Each land pixel is assigned to one of 30 environmental clusters (defined by land cover, Köppen climate, and soil texture). A cluster-specific RF model predicts soil moisture from 9 TB channels.

**Two assignment modes:**

| Mode | What happens | Required artifacts |
|------|-------------|-------------------|
| **Hard** (`--hard`) | Reads cluster ID directly from `environmental_clusters_k30.nc` | Cluster grid + RF models |
| **Soft** (default) | Computes GMM posterior weights, blends top-K cluster predictions | Cluster grid + GMM + environmental layers + encoders + RF models |

Soft blending eliminates boundary artifacts but requires re-encoding environmental layers and computing membership weights at startup (~30s one-time cost).

---

## Directory Layout

```
/discover/nobackup/ejalilva/dev/wsf_ml_sm/training/stratified_rf/
├── environmental_clusters_k30.nc     # Pre-computed cluster map (2560×1920)
├── clustering_model_k30.pkl          # KMeans + encoders (needed for soft mode)
└── cluster_models/
    ├── gmm_k30.pkl                   # GMM for soft assignment weights
    ├── rf_global.pkl                 # Global fallback RF (~16 GB)
    ├── rf_cluster_00.pkl             # Per-cluster RFs (00–29, ~54 GB total)
    │   ...
    └── rf_cluster_29.pkl
    └── stratified_config.json            # K, features, training metadata

```

---

## Scripts

### `predict_sm_stratified.py` — Primary (LIS/LDT Integration)

Class-based (`StratifiedSoilMoisturePredictor`). Two input modes:

- **`hourly`** — reads a single resampled TB NetCDF, writes one `ARFS_SM_WSFM_YYYYMMDDTHHMMSS.nc` per timestep (operational use).
- **`zarr`** — reads a monthly zarr archive, writes a zarr dataset (batch/research).

### `predict_blended_sm.py` — Standalone Research

Functional style, processes date ranges from zarr archives. Same blending logic.

---

## Usage

```bash
# Hourly (operational)
python predict_sm_stratified.py hourly /path/to/WSF_SDR_resampled_20240413_t1200.nc --output-dir /out

# Monthly zarr
python predict_sm_stratified.py zarr 202401_DES

# Hard assignment (faster, uses cluster grid directly)
python predict_sm_stratified.py zarr 202401_DES --hard

# Custom blending
python predict_sm_stratified.py hourly input.nc --top-k 5 --temperature 0.5

# Standalone: single date or range
python predict_blended_sm.py --date 2025-04-15
python predict_blended_sm.py --date-range 2025-04-01 2025-04-30 --hard
```

```python
# As module
from predict_sm_stratified import StratifiedSoilMoisturePredictor

predictor = StratifiedSoilMoisturePredictor(top_k=3, hard=False)
predictor.predict_from_hourly_nc('/path/to/input.nc', output_dir='/path/to/out')
```

---

## Input Features

**RF predictors (9 TB channels):** `TB_10V`, `TB_10H`, `TB_18V`, `TB_18H`, `TB_23V`, `TB_36V`, `TB_36H`, `TB_89V`, `TB_89H`

**Clustering features (not used by RF, only for cluster assignment):** land cover (one-hot), Köppen climate group (30→10 groups, one-hot), clay/sand/silt fraction (standardized). Only needed in soft blending mode.

---

## Output Conventions

Must match the LIS Fortran reader (`read_WSFsm.F90`):

- Variable: `arfs_sm`
- Dims: `(time=1, lat, lon)` — one timestep per file
- Grid: 2560 × 1920 (~0.14° × 0.09°)
- Attribute: `MAP_PROJECTION = 'EQUIDISTANT CYLINDRICAL'`
- Naming: `ARFS_SM_WSFM_YYYYMMDDTHHMMSS.nc`
- Valid range: [0.0, 0.6] m³/m³

---

## LDT/LIS Integration

1. **LDT ARFS resampling** → `WSF_SDR_resampled_YYYYMMDD_tHH00.nc` (hourly TB on ARFS grid)
2. **`predict_sm_stratified.py hourly`** → `ARFS_SM_WSFM_YYYYMMDDTHHMMSS.nc` (soil moisture)
3. **LIS `read_WSFsm`** → ingests `arfs_sm` for data assimilation

---

## QC & Fallbacks

- Ocean pixels masked via bit 0 of `QUALITY_FLAG`
- Pixels with any NaN TB channel are skipped
- Pixels outside environmental grid coverage are skipped
- Missing per-cluster model → global RF fallback

---

## Dependencies

`numpy`, `xarray`, `scipy`, `joblib`, `scikit-learn` (must match training version), `zarr`, `pandas` (standalone script only), `dask` (optional for chunked loading)
