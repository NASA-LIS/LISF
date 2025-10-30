#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
SCRIPT: run_prediction.py

Script for invoking AI/ML algorithms to retrieve snow depth from AMSR2
data.

REVISION HISTORY:
15 Aug 2025: Kehan Yang, Initial specification
18 Aug 2025: Eric Kemp, Code cleanup.
30 Oct 2025: Eric Kemp, More code cleanup.
"""

# Standard modules
from datetime import datetime
import logging
import os
import time
from typing import Optional, Tuple

# Third party modules
import numpy as np
import pandas as pd
import xarray as xr
import xgboost as xgb

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('SnowDepthPredictor')

# Pylint flags catching general exceptions, which is excessive.  We disable
# that here.
# pylint: disable=W0718

class SnowDepthPredictor:
    """
    Predicts snow depth using machine learning models with passive
    microwave data.
    """

    def __init__(self, config=None):
        """
        Initialize the snow depth predictor.

        Args:
            config: Configuration dictionary containing paths and parameters
            target_datetime: The datetime for which to predict snow depth
        """
        self.config = config
        self.target_datetime = self.config.target_datetime
        self.model = None
        self.pmw_file = None
        self.model_feature_names = None
        self.model_channels = [
            '6.9 GHz H', '6.9 GHz V',
            '7.3 GHz H', '7.3 GHz V',
            '10.7 GHz H', '10.7 GHz V',
            '18.7 GHz H', '18.7 GHz V',
            '23.8 GHz H', '23.8 GHz V',
            '36.6 GHz H', '36.6 GHz V',
            '89.0 GHz H', '89.0 GHz V'
        ]

        self.tb_channels = ['tb_6h', 'tb_6v',
                            'tb_7h', 'tb_7v',
                            'tb_10h','tb_10v',
                            'tb_18h', 'tb_18v',
                            'tb_23h', 'tb_23v',
                            'tb_36h','tb_36v',
                            'tb_89h', 'tb_89v']

        self.channel_mapping = {
            'tb_6h': '6.9 GHz H',
            'tb_6v': '6.9 GHz V',
            'tb_7h': '7.3 GHz H',
            'tb_7v': '7.3 GHz V',
            'tb_10h': '10.7 GHz H',
            'tb_10v': '10.7 GHz V',
            'tb_18h': '18.7 GHz H',
            'tb_18v': '18.7 GHz V',
            'tb_23h': '23.8 GHz H',
            'tb_23v': '23.8 GHz V',
            'tb_36h': '36.6 GHz H',
            'tb_36v': '36.6 GHz V',
            'tb_89h': '89.0 GHz H',
            'tb_89v': '89.0 GHz V'
        }

    def load_model(self) -> None:
        """
        Load the pre-trained ML model from disk.

        Raises:
            FileNotFoundError: If the model file doesn't exist
            Exception: For other errors during model loading
        """
        try:
            start_time = time.time()
            model_path = self.config.project_path / self.config.model_path
            if not os.path.exists(model_path):
                raise FileNotFoundError(f"Model file not found: {model_path}")

            with open(model_path, 'rb') as f:
                model = xgb.XGBRegressor()
                model.load_model(model_path)
                self.model = model
                self.model_feature_names = model.feature_names_in_.tolist()
            txt = f"Model loaded successfully from {model_path}" + \
                f" in {time.time() - start_time:.2f} seconds"
            logger.info(txt)
        except FileNotFoundError as e:
            logger.error("Model file not found: %s", e)
            raise
        except Exception as e:
            logger.error("Error loading model: %s", e)
            raise

    def get_file_paths(self) -> Tuple[str, str]:
        """
        Generate file paths for input and output data based on target date.

        Returns:
            Tuple containing output_file path, PMW data path, and target_datetime string
        """
        target_datetime = self.target_datetime.strftime("%Y%m%d%H")

        dir_out =  self.config.project_path / self.config.output_dir
        os.makedirs(dir_out, exist_ok=True)

        # Update output file path to use subfolder (keeping original
        # filename)
        output_file = os.path.join(dir_out, \
                                f'amsr2_snip_0p1deg_{target_datetime}.nc')

        return output_file, target_datetime

    def predict_snow_depth(self, pmw_file) -> Optional[xr.DataArray]:
        """
        Run the full ML prediction pipeline to generate snow depth
        estimates.

        Returns:
            xr.DataArray: Predicted snow depth data array, or None if
                          processing fails
        """
        output_file, target_datetime = self.get_file_paths()
        self.pmw_file = pmw_file
        logger.info("Processing snow depth for %s", target_datetime)

        # Check if output already exists
        if os.path.exists(output_file):
            logger.info("Output file %s already exists, skipping processing", output_file)
            return None

        # Find input data file
        if not pmw_file:
            logger.warning("No PMW data found for %s", target_datetime)
            return None

        try:
            # Process PMW data
            logger.info('Reading TB data from %s', pmw_file)
            return self._process_pmw_data(pmw_file, target_datetime)
        except Exception as e:
            logger.error("Error processing PMW data: %s", e)
            return None

    def _process_pmw_data(self, pmw_file: str, \
                          target_datetime: str) -> Optional[xr.DataArray]:
        """
        Process passive microwave data and generate snow depth prediction.

        Args:
            pmw_file: Path to the PMW data file
            target_datetime: Target datetime string

        Returns:
            xr.DataArray: Snow depth predictions
        """
        # Open PMW dataset
        with xr.open_dataset(pmw_file) as ds_pmw:
            logger.info('%s file opened', pmw_file)

            # Extract data for model input
            df = self._extract_features(ds_pmw, target_datetime)

            # Apply model to make predictions
            y_pred = self._apply_model(df)
            if y_pred is None:
                return None

            # Reshape predictions to original dimensions
            lat = ds_pmw.squeeze().lat.values
            lon = ds_pmw.squeeze().lon.values
            output_shape = (len(lat), len(lon))

            # File will be automatically closed when exiting the 'with'
            # block
            return self._format_output(y_pred, output_shape, lat, lon)

    def add_days_since_wy(self, df):
        """Add days from start of water year, and sine and cosine
        of that value"""
        df = df.copy()
        df['date'] = pd.to_datetime(df['Date'], format="%Y%m%d%H")

        # Northern Hemisphere
        nh_mask = df['lat'] >= 0
        nh_oct_plus = (df['date'].dt.month >= 10) & nh_mask
        nh_jan_sep = (df['date'].dt.month < 10) & nh_mask

        # Southern Hemisphere
        sh_mask = df['lat'] < 0
        sh_apr_plus = (df['date'].dt.month >= 4) & sh_mask
        sh_jan_mar = (df['date'].dt.month < 4) & sh_mask

        # Calculate water year starts
        wy_starts = pd.Series(index=df.index, dtype='datetime64[ns]')

        wy_starts[nh_oct_plus] = pd.to_datetime(
            df.loc[nh_oct_plus, 'date'].dt.year.astype(str) + '-10-01')
        wy_starts[nh_jan_sep] = pd.to_datetime(
            (df.loc[nh_jan_sep, 'date'].dt.year - 1).astype(str) + '-10-01')
        wy_starts[sh_apr_plus] = pd.to_datetime(
            df.loc[sh_apr_plus, 'date'].dt.year.astype(str) + '-04-01')
        wy_starts[sh_jan_mar] = pd.to_datetime(
            (df.loc[sh_jan_mar, 'date'].dt.year - 1).astype(str) + '-04-01')

        df['days_since_wy'] = (df['date'] - wy_starts).dt.days
        # Add cyclical features
        df['sin_period'] = np.sin(2 * np.pi * df['days_since_wy'] / \
                                  365.25)
        df['cos_period'] = np.cos(2 * np.pi * df['days_since_wy'] / \
                                  365.25)

        return df
    def _extract_features(self, data: xr.Dataset, \
                          target_datetime: str) -> pd.DataFrame:
        """
               Extract features from xarray data for model input.

               Args:
                   data: Input xarray data
                   target_datetime: Target datetime string

               Returns:
                   pd.DataFrame: DataFrame with features for model
        """

        df = pd.DataFrame()
        for tb_channel in self.tb_channels:
            # Get corresponding model channel
            model_channel = self.channel_mapping[tb_channel]
            if hasattr(data, tb_channel):
                # Access data directly and squeeze to remove singleton
                # time dimension
                tb_data = getattr(data, tb_channel).squeeze()
                df[model_channel] = tb_data.values.flatten()
            else:
                print(f"Warning: {tb_channel} not found in dataset")

        df['Date'] = target_datetime
        lat_coords = data['lat'].values
        lon_coords = data['lon'].values
        #Lon, Lat = np.meshgrid(lon_coords, lat_coords)
        #df['lat'] = Lat.flatten()
        df['lat'] = np.meshgrid(lon_coords, lat_coords)[1].flatten()

        df = self.add_days_since_wy(df)

        return df

    def _apply_model(self, df: pd.DataFrame) -> Optional[np.ndarray]:
        """
        Apply the ML model to make predictions.

        Args:
            df: Input features DataFrame

        Returns:
            np.ndarray: Array of predictions, or None if insufficient data
        """
        # load model
        self.load_model()

        # Select only the features used by the model
        df = df[self.model_feature_names]
        # Identify rows with NaN values
        nan_mask = np.isnan(df).any(axis=1)
        x_test_no_nan = df[~nan_mask]

        # Select only the features used by the model
        # x_test_no_nan = x_test_no_nan[self.model_feature_names]

        # Check if we have enough data to make predictions
        if len(x_test_no_nan) < 100:
            txt = "Insufficient valid data points for prediction " + \
                "(less than 100)"
            logger.warning(txt)
            return None

        # Make predictions
        logger.info("Making predictions on %s valid data points", \
                    len(x_test_no_nan))
        y_pred_no_nan = self.model.predict(x_test_no_nan)

        # Reconstruct full prediction array with NaNs for invalid inputs
        y_pred = np.full(df.shape[0], np.nan)
        y_pred[~nan_mask] = y_pred_no_nan

        # Apply threshold to remove very small values
        threshold = 0.01
        y_pred[y_pred <= threshold] = 0

        return y_pred

    def _format_output(self, predictions: np.ndarray, shape: Tuple, \
                       y_coords: np.ndarray,
                       x_coords: np.ndarray) -> xr.DataArray:
        """
        Format predictions into an xarray DataArray with proper
        coordinates.

        Args:
            predictions: Array of predictions
            shape: Shape to reshape predictions to
            y_coords: Y coordinates
            x_coords: X coordinates

        Returns:
            xr.DataArray: Formatted output
        """
        # Reshape predictions to the target grid shape
        reshaped_predictions = predictions.reshape(shape)

        # Handle NaN values and invalid predictions
        # Set negative values to 0 or NaN if appropriate for snow depth
        reshaped_predictions = np.where(reshaped_predictions < 0, \
                                        np.nan, reshaped_predictions)

        # Create DataArray with coordinates
        output = xr.DataArray(
            reshaped_predictions,
            coords={
                'y': ('y', y_coords, {
                    'long_name': 'Latitude',
                    'standard_name': 'latitude',
                    'units': 'degrees_north'
                }),
                'x': ('x', x_coords, {
                    'long_name': 'Longitude',
                    'standard_name': 'longitude',
                    'units': 'degrees_east'
                })
            },
            dims=["y", "x"],
            attrs={
                'long_name': 'Snow Depth',
                'units': 'meters',
                'valid_range': [0.0, 10.0]
            }
        )

        # Set coordinate reference system (CRS)
        output = output.rio.write_crs(self.config.proj, inplace=True)

        return output



    def save_to_netcdf(self, output: xr.DataArray, pmw_file: str) -> bool:
        """
           Save the snow depth output to a NetCDF file.

           Args:
               output: Snow depth data array to save
               pmw_file: Path to the PMW file for masking operations

           Returns:
               bool: True if successful, False otherwise
           """
        if output is None:
            logger.warning("No output data to save")
            return False

        output_file, _ = self.get_file_paths()

        try:
            # Create dataset with metadata
            ds_out = output.to_dataset(name="snow_depth")

            # Verify we have the expected coordinates
            if 'y' not in output.coords or 'x' not in output.coords:
                txt = "Expected coordinates 'y' and 'x' in output " + \
                    "DataArray"
                raise ValueError(txt)

            # Update coordinate attributes directly
            ds_out.coords['y'].attrs.update({
                'long_name': 'Latitude',
                'standard_name': 'latitude',
                'units': 'degrees_north',
                'valid_range': [-90., 90.],
                'axis': 'Y'
            })

            ds_out.coords['x'].attrs.update({
                'long_name': 'Longitude',
                'standard_name': 'longitude',
                'units': 'degrees_east',
                'valid_range': [-180., 180.],
                'axis': 'X'
            })

            # Add global attributes
            ds_out.attrs.update({
                'title': 'AMSR2 Snow Depth',
                'description': \
                   'Snow depth derived from AMSR2 passive microwave data',
                'creation_date': \
                   datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'source': 'XGBoost model prediction',
                'Conventions': 'CF-1.8',
                'institution': 'NASA GSFC',
                'geospatial_lat_min': float(ds_out.y.min()),
                'geospatial_lat_max': float(ds_out.y.max()),
                'geospatial_lon_min': float(ds_out.x.min()),
                'geospatial_lon_max': float(ds_out.x.max()),
                'geospatial_lat_resolution': \
                    float(abs(ds_out.y.diff('y').mean())) \
                    if len(ds_out.y) > 1 else 0.0,
                'geospatial_lon_resolution': \
                    float(abs(ds_out.x.diff('x').mean())) \
                    if len(ds_out.x) > 1 else 0.0
            })

            # Add variable attributes
            ds_out['snow_depth'].attrs.update({
                'units': 'meters',
                'long_name': 'Snow Depth',
                'standard_name': 'surface_snow_thickness',
                'valid_range': [0.0, 10.0],
                'grid_mapping': 'crs'
            })

            # Set encoding for better compression and type handling
            encoding = {
                'snow_depth': {
                    'zlib': True,
                    'complevel': 4,
                    'shuffle': True,
                    '_FillValue': -9999.0,
                    'dtype': 'float32'
                },
                'y': {'dtype': 'float32'},
                'x': {'dtype': 'float32'}
            }

            # Create parent directories if they don't exist
            os.makedirs(os.path.dirname(output_file), exist_ok=True)

            # open pmw file
            ds_pmw = xr.open_dataset(pmw_file).squeeze()

            # Apply land ocean frac threshold - only predict land
            if self.config.land_frac_th is not None:
                txt = "Apply land water fraction flag: " + \
                    f"{self.config.land_frac_th}%"
                logger.info(txt)
                flag_layer = 'land_water_frac'
                # Ensure proper alignment and dimensions
                land_frac = ds_pmw[flag_layer].squeeze()

                # Align coordinates if necessary
                if land_frac.dims == ('lat', 'lon'):
                    land_frac = land_frac.rename({'lat': 'y', 'lon': 'x'})

                land_mask = land_frac >= self.config.land_frac_th

                # Apply mask maintaining original dimensions
                land_mask = land_mask.squeeze()
                ds_out['snow_depth'] = ds_out['snow_depth']. \
                    where(land_mask)

            # apply flag
            if self.config.flag_cold:
                #  Flag rain, cold deserts, frozen ground, and glacier
                txt = "Apply cold deserts, frozen ground, and " + \
                    "glacier flag"
                logger.info(txt)

                ds_pmw = xr.open_dataset(pmw_file)
                flag_layer = 'pixel_qual_flag'

                # Set bits for each condition (1 means condition is
                # present)
                # Bit 0: rain
                # Bit 1: cold deserts
                # Bit 2: frozen ground
                # Bit 3: glacier
                # Bit 4: snow cover
                # Bit 5: precip
                # Bit 6: rfi_h
                # Bit 7: rfi_v

                # Create mask for pixels with any of the specified flag
                # values
                flag_values_to_mask = [1, 2, 3, 5]  # Skip bit 4 (snow cover)
                # Add rain flag if configured
                if self.config.flag_rain:
                    logger.info("Apply rain flag")
                    flag_values_to_mask.append(0)  # Bit 0: rain

                # Add RFI flags if configured
                if self.config.flag_rfi:
                    logger.info("Apply RFI h and v flag for C bands")
                    flag_values_to_mask.extend([6, 7])  # Bit 6: rfi_h, Bit 7: rfi_v

                flag_cleaned = ds_pmw[flag_layer].fillna(255)
                flag_uint8 = flag_cleaned.astype('uint8')

                flag_mask = xr.zeros_like(flag_uint8, dtype=bool)
                for bit in flag_values_to_mask:
                    bit_mask = (flag_uint8 & (1 << bit)) != 0
                    flag_mask = flag_mask | bit_mask

                txt = f"{flag_layer}: {flag_mask.sum().values} " + \
                    "pixels flagged"
                logger.info(txt)
                txt = f"Total masked pixels: {flag_mask.sum().values}"
                logger.info(txt)
                txt = f"Total valid pixels: {(~flag_mask).sum().values}"
                logger.info(txt)

                flag_mask = flag_mask.squeeze()
                if flag_mask.dims == ('lat', 'lon'):
                    flag_mask = flag_mask.rename({'lat': 'y', 'lon': 'x'})

                # Apply mask to your output dataset
                ds_out = ds_out.where(~flag_mask.squeeze())

            # Save to NetCDF
            ds_out.to_netcdf(output_file, encoding=encoding, \
                             format='NETCDF4')
            logger.info("Output saved to %s", output_file)
            return True

        except Exception as e:
            logger.error("Error saving output: %s", e)
            return False

    def reproject_to_usaf(self) -> bool:
        """
        reproject snow depth output based on a usaf nc file.

        Returns:
            bool: True if successful, False otherwise
        """
        output_file, _ = self.get_file_paths()
        # save reprojected data to netcdf
        base_name = os.path.splitext(os.path.basename(output_file))[0]
        filename = f"{base_name}00_AFgrid.nc"
        full_path = os.path.join(os.path.dirname(output_file), filename)
        if not os.path.exists(full_path):
            try:
                template_path = self.config.project_path / \
                    self.config.template_path
                template_data = xr.open_dataset(template_path)

                # open snow depth data
                ds_sd = xr.open_dataset(output_file)

                # reproject snow depth data based on template
                # Ensure SCA has CRS
                template_data = template_data.rio. \
                    write_crs(self.config.proj)

                # Assuming output is in WGS84 lat/lon based on your
                # coordinates
                ds_sd = ds_sd.rio.write_crs(self.config.proj)

                ds_sd_reprojected = ds_sd.rio. \
                    reproject_match(template_data)

                ds_out = xr.Dataset(
                    data_vars={
                        'snowdepth': (
                            ['lat', 'lon'],
                            ds_sd_reprojected['snow_depth'].data,
                            {
                                'units': 'meters',  # or 'cm', adjust as needed
                                'long_name': 'Snow Depth',
                                'standard_name': 'snow_depth',
                                '_FillValue': -9999.0,
                                'valid_range': [0.0, 10.0]  # adjust range as appropriate
                            }
                        )
                    },
                    coords={
                        'lat': (
                            ['lat'],
                            ds_sd_reprojected.y.data,
                            {
                                'units': 'degrees_north',
                                'long_name': 'Latitude',
                                'standard_name': 'latitude'
                            }
                        ),
                        'lon': (
                            ['lon'],
                            ds_sd_reprojected.x.data,
                            {
                                'units': 'degrees_east',
                                'long_name': 'Longitude',
                                'standard_name': 'longitude'
                            }
                        )
                    },
                    attrs={
                        'map_projection': 'EQUIDISTANT CYLINDRICAL',
                        'title': 'AMSR2 Snow Depth',
                        'description': 'Snow depth derived from AMSR2 passive microwave data',
                        'creation_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        'DX': 0.140625,
                        'DY': 0.09375,
                        'source': 'XGBoost model prediction',
                        'Conventions': 'CF-1.8',
                        'institution': 'NASA GSFC Hydrological Sciences Laboratory',
                    }
                )

                # Save to NetCDF
                ds_out.to_netcdf(full_path)
                logger.info("Reprojected file saved to %s", full_path)
                return True

            except Exception as e:
                logger.error("Error in reproject_to_usaf: %s", e)
                return False
        else:
            logger.info("Reprojected file exist: %s", full_path)
            return True


    def run_pipeline(self, pmw_file) -> bool:
        """
        Run the complete snow depth prediction pipeline.

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Check if output already exists
            output_file, _ = self.get_file_paths()
            if os.path.exists(output_file):
                txt = f"Output file {output_file} already exists, " + \
                      "skipping processing"
                logger.info(txt)

                # check if AFgrid data exists
                base_name = os.path.splitext(os.path.basename(output_file))[0]
                filename = f"{base_name}_AFgrid.nc"
                full_path = os.path.join(os.path.dirname(output_file), filename)
                if os.path.exists(full_path):
                    txt = f"AF grid output file {full_path} already exists, " + \
                          "skipping processing"
                    logger.info(txt)
                else:
                    logger.info("Output at AF grid not exists;"
                                "Reproject to USAF grid")
                    self.reproject_to_usaf()
                return True

            # Generate initial prediction
            start_time = time.time()
            output = self.predict_snow_depth(pmw_file=pmw_file)
            if output is None:
                return False

            # Save results
            success = self.save_to_netcdf(output=output,
                                          pmw_file=pmw_file)

            if not success:
                txt = "Failed to save output to NetCDF. Pipeline aborted."
                logger.error(txt)
                return False

            # Reproject to USAF grid if configured
            if self.config.reproject_USAF:
                reproject_success = self.reproject_to_usaf()

                if reproject_success:
                    logger.info(
                        "USAF reprojection completed successfully")
                else:
                    txt = "USAF reprojection failed, but " + \
                          "primary output was saved successfully"
                    logger.warning(txt)
            txt = "Pipeline completed in " + \
                f"{time.time() - start_time:.2f} seconds"
            logger.info(txt)
            return True

        except Exception as e:
            logger.error("Error in pipeline: %s", e)
            return False
