#!/usr/bin/env python3

# -----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# -------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
SCRIPT: run_prediction_WSF.py

Script for invoking the Two-Stage AI/ML algorithm to retrieve snow depth
from WSF data.

REVISION HISTORY:
12 Apr 2026: Kehan Yang, Initial specification
"""
# pylint: disable=import-error
import glob
import logging
import os
import re
import time
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import xarray as xr
import xgboost as xgb

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('WSFSnowDepthPredictor')

# Encoding used for all NetCDF outputs
_ENCODING = {
    'snowdepth': {
        'zlib': True,
        'complevel': 4,
        'shuffle': True,
        '_FillValue': -9999.0,
        'dtype': 'float32'
    },
    'lat': {'dtype': 'float32'},
    'lon': {'dtype': 'float32'}
}


# pylint: disable=too-many-instance-attributes, invalid-name,
# too-many-locals, too-many-statements, too-many-branches
class WSFSnowDepthPredictor:
    """
    Predicts snow depth using Two-Stage machine learning models with WSF
    passive microwave data.
    """

    def __init__(self, config=None):
        self.config = config
        self.target_datetime = self.config.target_datetime

        self.models = {}
        self.model_feature_names = None

        # WSF specific channels
        self.model_channels = [
            'TB_10H', 'TB_10V', 'TB_18H', 'TB_18V',
            'TB_23V', 'TB_36H', 'TB_36V', 'TB_89H', 'TB_89V'
        ]

    TB_DIFF_PAIRS = {
        '10_V-18_V': ('TB_10V', 'TB_18V'),
        '10_V-23_V': ('TB_10V', 'TB_23V'),
        '10_V-36_V': ('TB_10V', 'TB_36V'),
        '10_V-89_V': ('TB_10V', 'TB_89V'),
        '18_V-23_V': ('TB_18V', 'TB_23V'),
        '18_V-36_V': ('TB_18V', 'TB_36V'),
        '18_V-89_V': ('TB_18V', 'TB_89V'),
        '23_V-36_V': ('TB_23V', 'TB_36V'),
        '23_V-89_V': ('TB_23V', 'TB_89V'),
        '36_V-89_V': ('TB_36V', 'TB_89V'),
    }

    def _add_tb_differences(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add brightness temperature difference features."""
        for diff_name, (ch1, ch2) in self.TB_DIFF_PAIRS.items():
            if ch1 in df.columns and ch2 in df.columns:
                df[diff_name] = df[ch1] - df[ch2]
            else:
                logger.warning('Skipping diff %s: missing channels',
                               diff_name)
        return df

    def load_model(self) -> None:
        """Load the pre-trained Two-Stage ML model from disk."""
        try:
            start_time = time.time()
            base_model_path = str(
                self.config.project_path / self.config.model_path)

            base_model_path = base_model_path + '_WSF'
            # Load Classifier
            clf_path = f"{base_model_path}_classifier.json"
            if os.path.exists(clf_path):
                clf = xgb.XGBClassifier()
                clf.load_model(clf_path)
                self.models['classifier'] = clf
                self.model_feature_names = clf.get_booster().feature_names
            else:
                raise FileNotFoundError(f"Classifier not found at {clf_path}")

            # Load Regressors
            for reg_name in ['reg_shallow', 'reg_other']:
                reg_path = f"{base_model_path}_{reg_name}.json"
                if os.path.exists(reg_path):
                    reg = xgb.XGBRegressor()
                    reg.load_model(reg_path)
                    self.models[reg_name] = reg

            logger.info(
                "Two-Stage WSF Models loaded successfully in %.2f seconds",
                time.time() - start_time)

        except Exception as e:
            logger.error("Error loading WSF two-stage model: %s", e)
            raise

    def get_final_file_path(self) -> str:
        """Generate the exact file path for LDT bash script compatibility."""
        target_datetime_str = self.target_datetime.strftime("%Y%m%d%H%M")
        dir_out = self.config.project_path / self.config.output_dir
        os.makedirs(dir_out, exist_ok=True)
        return os.path.join(dir_out,
                            f'wsf_snip_0p1deg_{target_datetime_str}_AFgrid.nc')

    def get_hourly_file_path(self, pmw_file: str) -> str:
        """Generate a temporary filename for the hourly prediction."""
        dir_out = self.config.project_path / self.config.output_dir
        os.makedirs(dir_out, exist_ok=True)

        # Extract date and time from the PMW file name
        # (e.g., WSF_SDR_resampled_allfp_20250120_t0600_DES.nc)

        # Matches an 8-digit date followed by _t and a 4-digit time
        match = re.search(r'(\d{8})_t(\d{4})', os.path.basename(pmw_file))

        if match:
            return os.path.join(dir_out,
                                f"wsf_snip_0p1deg_{match.group(1)}_"
                                f"t{match.group(2)}.nc")

        return os.path.join(dir_out, f"temp_wsf_pred_{int(time.time())}.nc")

    def apply_flags(self, tb_data, land_water_frac, n, m):
        """Generate snow and precipitation flags based on TB values."""

        def get_tb(key):
            if key in tb_data:
                val = tb_data[key].squeeze().values if (
                    hasattr(tb_data[key], 'values')) else np.array(tb_data[key])
                val = val.astype(np.float32).copy()
                val[val <= -9999.0] = np.nan
                return val
            return np.full((n, m), np.nan, dtype=np.float32)

        tb_18v, tb_18h, tb_23v, tb_36v, tb_89v = get_tb('TB_18V'), get_tb(
            'TB_18H'), get_tb('TB_23V'), get_tb('TB_36V'), get_tb('TB_89V')

        lw = land_water_frac.squeeze().values if hasattr(
            land_water_frac,
            'values') else np.array(
            land_water_frac)
        lw = lw.astype(np.float32).copy()
        lw[lw <= -9999.0] = np.nan

        rain, cold_deserts, frozen_ground = np.zeros((n, m),
                                                     dtype=np.int32), np.zeros(
            (n, m), dtype=np.int32), np.zeros((n, m), dtype=np.int32)
        glacier, snow, precip = np.zeros((n, m), dtype=np.int32), np.zeros(
            (n, m), dtype=np.int32), np.zeros((n, m), dtype=np.int32)

        land_mask = lw >= 0.9
        valid_mask = (np.isfinite(tb_18v) & np.isfinite(tb_18h) & np.isfinite(
            tb_23v) & np.isfinite(tb_36v) & np.isfinite(tb_89v))
        base_mask = land_mask & valid_mask

        scat = np.maximum(tb_23v - tb_89v, tb_18v - tb_36v)
        sc37 = tb_18v - tb_36v
        pd19 = tb_18v - tb_18h
        scx = tb_36v - tb_89v

        tt = 165.0 + 0.49 * tb_89v
        sil = (451.88 - 0.44 * tb_18v - 1.775 * tb_23v +
               0.00574 * tb_23v ** 2 - tb_89v)

        rain[base_mask & (
                ((tb_23v >= 254.0) & (scat <= 2.0)) | (tb_23v >= 258.0) | (
                tb_23v >= tt))] = 1
        cold_deserts[
            base_mask & ((pd19 >= 18.0) & (sc37 <= 10.0) & (scx <= 10.0))] = 1
        frozen_ground[base_mask & ((scat <= 6.0) & (pd19 >= 8.0))] = 1
        glacier[base_mask & ((tb_23v <= 210.0) | (
                (tb_23v <= 229.0) & (pd19 >= 23.0)))] = 1

        sil_mask = base_mask & (sil > 10)
        snow_cond = ((tb_23v <= 264.0) & (tb_23v <= (175.0 + 0.49 * tb_89v)))

        snow[sil_mask & snow_cond] = 1
        snow[sil_mask & snow_cond & (
                (pd19 >= 18) & ((tb_18v - tb_36v) <= 10) & (
                (tb_36v - tb_89v) <= 10))] = 0
        snow[sil_mask & snow_cond & ((pd19 >= 8) & ((tb_18v - tb_36v) <= 2) & (
                (tb_23v - tb_89v) <= 6))] = 1

        precip_cond = ~snow_cond
        precip[sil_mask & precip_cond] = 1
        precip[sil_mask & precip_cond & (pd19 > 20)] = 0
        precip[sil_mask & precip_cond & ((tb_89v > 253) & (pd19 > 7))] = 0

        return rain, cold_deserts, frozen_ground, glacier, snow, precip

    def apply_filter(self, data):
        """Apply traditional WSF filter thresholds."""
        is_coldsnow = (data['TB_36H'] < 245) & (data['TB_36V'] < 255)
        is_medium_deep = (data['TB_10V'] > data['TB_36V']) | (
                data['TB_10V'] > data['TB_36H'])
        T = 58.08 - 0.39 * data['TB_18V'] + 1.21 * data['TB_23V'] - 0.37 * data[
            'TB_36H'] + 0.36 * data['TB_89V']
        is_shallow = ((data['TB_89V'] <= 255) & (data['TB_89H'] <= 265) & (
                data['TB_23V'] > data['TB_89V']) & (T < 267))

        depth_flag = np.full_like(data['TB_36H'].values, np.nan, dtype=float)
        valid_mask = ~np.isnan(data['TB_36H'].values)

        depth_flag[valid_mask] = 0
        depth_flag = np.where(np.array(is_shallow) & valid_mask, 2, depth_flag)
        depth_flag = np.where(np.array(is_medium_deep) & valid_mask, 1,
                              depth_flag)
        depth_flag = np.where(np.array(is_coldsnow) & valid_mask, 1, depth_flag)

        return depth_flag

    def _extract_features(self, data: xr.Dataset, target_datetime: str) -> \
            Tuple[pd.DataFrame, np.ndarray]:
        """Extract WSF features and apply combined masks."""
        df = pd.DataFrame()
        nlat, nlon = data['lat'].shape[0], data['lon'].shape[0]

        land_frac = data['LAND_FRAC'].squeeze().values.flatten()
        land_frac[land_frac <= -9999.0] = np.nan
        land_mask = (land_frac < 0.9) | np.isnan(land_frac)

        _, cold_deserts, frozen_ground, glacier, _, precip = self.apply_flags(
            data, data['LAND_FRAC'], nlat, nlon)

        combined_mask = (land_mask | (cold_deserts.squeeze().flatten() == 1) | (
                frozen_ground.squeeze().flatten() == 1) |
                         (glacier.squeeze().flatten() == 1) | (
                                 precip.squeeze().flatten() == 1))

        lon_coords, lat_coords = data['lon'].values, data['lat'].values
        Lon, Lat = np.meshgrid(lon_coords, lat_coords)
        df['lat'], df['lon'] = Lat.flatten(), Lon.flatten()
        df['Date'] = target_datetime
        df['date'] = pd.to_datetime(df['Date'], format="%Y%m%d%H")

        for channel in self.model_channels:
            if channel in data:
                values = data[channel].squeeze().values.flatten()
                values[values <= -9999.0] = np.nan
                values[combined_mask] = np.nan
                df[channel] = values

        return self._add_tb_differences(df), combined_mask

    def _apply_two_stage_model(self, df: pd.DataFrame) -> Optional[np.ndarray]:
        """Apply the Two-Stage ML model to make predictions."""
        if not self.models:
            self.load_model()

        missing = [f for f in self.model_feature_names if f not in df.columns]
        if missing:
            logger.error(f"Missing features required by model: {missing}")
            return None

        df_sub = df[self.model_feature_names]
        nan_mask = np.isnan(df_sub).any(axis=1)
        x_test_no_nan = df_sub[~nan_mask]

        if len(x_test_no_nan) < 100:
            return None

        pred_class = self.models['classifier'].predict(x_test_no_nan)
        y_pred_no_nan = np.zeros(len(x_test_no_nan))

        mask_shallow = pred_class == 1
        if mask_shallow.sum() > 0 and 'reg_shallow' in self.models:
            y_pred_no_nan[mask_shallow] = self.models['reg_shallow'].predict(
                x_test_no_nan[mask_shallow])

        mask_other = pred_class == 2
        if mask_other.sum() > 0 and 'reg_other' in self.models:
            y_pred_no_nan[mask_other] = self.models['reg_other'].predict(
                x_test_no_nan[mask_other])

        y_pred_no_nan = np.clip(y_pred_no_nan, 0, None)
        y_pred = np.full(df.shape[0], np.nan)
        y_pred[~nan_mask] = y_pred_no_nan

        return y_pred

    def _format_output(self, predictions: np.ndarray, ds_pmw: xr.Dataset,
                       combined_mask: np.ndarray) -> xr.DataArray:
        """Format predictions, apply post-filters, and convert to DataArray."""
        predictions[combined_mask] = np.nan
        predictions[predictions <= 0.025] = 0

        nlat, nlon = ds_pmw['lat'].shape[0], ds_pmw['lon'].shape[0]
        output_np = predictions.reshape(nlat, nlon)

        # Apply Physical WSF Filter
        # depth_flag = self.apply_filter(ds_pmw)
        # output_np[depth_flag == 0] = 0
        # output_np[depth_flag == 2] = 0.050

        output = xr.DataArray(
            output_np,
            coords={'lat': ds_pmw['lat'].values, 'lon': ds_pmw['lon'].values},
            dims=['lat', 'lon']
        )
        return output.rio.write_crs(self.config.proj, inplace=True)

    def run_hourly_prediction(self, pmw_file: str) -> Optional[str]:
        """Generate prediction for a single file and save directly to disk."""
        hourly_out_path = self.get_hourly_file_path(pmw_file)

        try:
            with xr.open_dataset(pmw_file, decode_timedelta=False) as ds_pmw:
                ds_pmw = ds_pmw.rio.write_crs("EPSG:4326", inplace=True)

                # Extract features
                df, combined_mask = self._extract_features(
                    ds_pmw,
                    self.target_datetime.strftime(
                        "%Y%m%d%H"))

                # Predict
                y_pred = self._apply_two_stage_model(df)
                if y_pred is None:
                    return None

                # Format output
                da_output = self._format_output(y_pred, ds_pmw, combined_mask)

                # Save to disk immediately to clear memory
                # Create a properly formatted Dataset so the projection is recognized
                ds_out = xr.Dataset(
                    data_vars={
                        'snowdepth': (
                            ['lat', 'lon'],
                            da_output.values,
                            {
                                'units': 'meters',
                                'long_name': 'Snow Depth',
                                'standard_name': 'snow_depth',
                                'valid_range': [0.0, 10.0]
                            }
                        )
                    },
                    coords={
                        'lat': (['lat'], da_output.lat.values,
                                {'units': 'degrees_north',
                                 'long_name': 'Latitude'}),
                        'lon': (['lon'], da_output.lon.values,
                                {'units': 'degrees_east',
                                 'long_name': 'Longitude'})
                    },
                    attrs={
                        'map_projection': 'EQUIDISTANT CYLINDRICAL',
                        'title': 'WSF Snow Depth Hourly',
                        'source': 'Two-Stage XGBoost model prediction',
                        'Conventions': 'CF-1.8',
                        'institution': 'NASA GSFC'
                    }
                )

                if ds_out.lat.values[0] > ds_out.lat.values[-1]:
                    ds_out = ds_out.sortby(
                        'lat')  # Sorts from smallest (-90) to largest (90)

                ds_out = ds_out.rio.write_crs(self.config.proj, inplace=True)

                # Use _ENCODING here to apply the proper _FillValue and compression
                ds_out.to_netcdf(hourly_out_path, encoding=_ENCODING)

                return hourly_out_path

        except Exception as e:
            logger.error("Error processing WSF PMW data %s: %s",
                         pmw_file, e)
            return None

    def save_composite_to_netcdf(self, output: xr.DataArray) -> bool:
        """Save the final composite snow depth output to a NetCDF file."""
        if output is None:
            return False

        output_file = self.get_final_file_path()

        try:

            ds_out = xr.Dataset(
                data_vars={
                    'snowdepth': (
                        ['lat', 'lon'],
                        output.values,
                        {
                            'units': 'meters',
                            'long_name': 'Snow Depth',
                            'standard_name': 'snow_depth',
                            'valid_range': [0.0, 10.0]
                        }
                    )
                },
                coords={
                    'lat': (['lat'], output.lat.values,
                            {'units': 'degrees_north',
                             'long_name': 'Latitude'}),
                    'lon': (['lon'], output.lon.values,
                            {'units': 'degrees_east', 'long_name': 'Longitude'})
                },
                attrs={
                    'map_projection': 'EQUIDISTANT CYLINDRICAL',
                    'title': 'WSF Snow Depth',
                    'source': 'Two-Stage XGBoost model composite',
                    'Conventions': 'CF-1.8',
                    'institution': 'NASA GSFC'
                }
            )

            if ds_out.lat.values[0] > ds_out.lat.values[-1]:
                ds_out = ds_out.sortby(
                    'lat')  # Sorts from smallest (-90) to largest (90)

            ds_out = ds_out.rio.write_crs(self.config.proj, inplace=True)
            ds_out.to_netcdf(output_file, encoding=_ENCODING)
            logger.info("Final 6-hour composite saved to %s",
                        output_file)
            return True

        except Exception as e:
            logger.error("Error saving WSF output: %s", e)
            return False


class WSFSnowWorkflow:
    """Workflow manager for the WSF Snow Depth prediction pipeline."""

    def __init__(self, config):
        self.config = config
        self.predictor = WSFSnowDepthPredictor(config)
        self.target_datetime = self.config.target_datetime

        self.time_windows = {
            '0600': ['t0000', 't0100', 't0200', 't0300', 't0400', 't0500'],
            '1200': ['t0600', 't0700', 't0800', 't0900', 't1000', 't1100'],
            '1800': ['t1200', 't1300', 't1400', 't1500', 't1600', 't1700'],
            '0000': ['t1800', 't1900', 't2000', 't2100', 't2200', 't2300'],
        }

    def _find_pmw_files(self) -> list:
        """Find WSF PMW data files for the specific 6-hour window."""
        window_end = self.target_datetime.strftime("%H%M")
        if window_end not in self.time_windows:
            logger.error(
                "Invalid target hour %s. Must be 0000, "
                "0600, 1200, or 1800.",
                window_end)
            return []

        expected_hours = self.time_windows[window_end]

        data_date = self.target_datetime - pd.Timedelta(hours=1)
        date_str = data_date.strftime("%Y%m%d")
        year_month_str = data_date.strftime("%Y%m")
        # Directly access resampled_base from the merged config
        if not hasattr(self.config,
                       'resampled_base') or not self.config.resampled_base:
            logger.error(
                "Configuration missing 'resampled_base'. "
                "Cannot find resampled WSF files.")
            return []

        wsf_resample_dir = self.config.resampled_base
        logger.info("Using resampled_base from config: %s",
                    wsf_resample_dir)

        search_pattern = os.path.join(
            self.config.project_path, '..', wsf_resample_dir, year_month_str,
            f"*{date_str}*_t*DES.nc"
        )

        all_files = sorted(glob.glob(search_pattern))

        window_files = [f for f in all_files if
                        any(h in os.path.basename(f) for h in expected_hours)]

        return window_files

    def run_workflow(self) -> bool:
        """Execute the full prediction workflow,
        composite the results, and save."""
        logger.info("Starting WSF workflow for %s",
                    self.target_datetime)

        final_output_file = self.predictor.get_final_file_path()
        if os.path.exists(final_output_file):
            logger.info(
                "Final output file %s already exists, skipping workflow.",
                final_output_file)
            return True

        try:
            pmw_files = self._find_pmw_files()
            if not pmw_files:
                logger.error("Required WSF input data not found. Aborting.")
                return False

            logger.info("Found %d input files to process.",
                        len(pmw_files))

            # 1. Predict and save each hour to disk
            hourly_files = []
            for pmw_file in pmw_files:
                logger.info("Predicting file: %s",
                            os.path.basename(pmw_file))
                hourly_path = self.predictor.run_hourly_prediction(pmw_file)
                if hourly_path:
                    hourly_files.append(hourly_path)

            if not hourly_files:
                logger.error(
                    "All predictions failed or returned no valid data.")
                return False

            # 2. Open hourly files from disk and composite them
            logger.info(
                "Compositing %d predicted swaths "
                "into a single 6-hour map...",
                len(hourly_files))

            # Start with the most recent file as the base composite
            with (xr.open_dataset(hourly_files[-1], decode_timedelta=False)
                  as ds):
                composite_da = ds['snowdepth'].load()

            # Layer the older files underneath (latest data takes priority)
            for f in reversed(hourly_files[:-1]):
                with xr.open_dataset(f, decode_timedelta=False) as ds:
                    da = ds['snowdepth'].load()
                    composite_da = composite_da.combine_first(da)

            # 3. Save final composite
            success = self.predictor.save_composite_to_netcdf(composite_da)

            # 4. Clean up temporary hourly files
            if success:
                logger.info("Cleaning up temporary hourly prediction files...")
                for f in hourly_files:
                    os.remove(f)
                logger.info("WSF composite workflow completed successfully.")

            return success

        except Exception as e:
            logger.error("Unexpected error in WSF workflow: %s", e)
            return False
