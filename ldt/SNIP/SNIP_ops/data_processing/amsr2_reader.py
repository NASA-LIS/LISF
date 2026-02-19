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
SCRIPT: amsr2_reader.py

Script for reading AMSR2 data files.

REVISION HISTORY:
15 Aug 2025: Kehan Yang, Initial Specification.
18 Aug 2025: Eric Kemp, Code cleanup.
"""

# Standard modules
from datetime import datetime, timedelta
import glob
import logging
import os
import time
# Third party modules
import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import cKDTree
import xarray as xr


# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('SnowDepthPredictor')

# Pylint flags catching and raising general exceptions.  This is too
# strict, so we disable the warnings.
# pylint: disable=W0718
# pylint: disable=W0719

class AMSR2DataProcessor:
    """Class for handling AMSR2 data"""
    def __init__(self, config=None):
        self.config = config
        self.target_resolution = self.config.target_resolution
        self.amsr2_files = []

    def process_l1r_data(self, target_datetime, max_retries=3):
        """Main processing pipeline for AMSR2 L1R data with retry logic"""
        for attempt in range(max_retries+1):
            try:
                logger.info("Processing attempt %s / %s for %s", \
                        (attempt + 1), (max_retries+1), target_datetime)

                # 1.1 Check available data
                available_files = self.check_available_data(target_datetime)

                # 1.2 Process each file
                processed_data = []
                for file_path in available_files:
                    tb_data = self.get_amsr_l1r(file_path)
                    # extract descending data only for NOAA/NESDIS AMSR2 L1R dataset
                    if self.config.source == 'NOAA':
                        tb_data = self.get_descending_orbit(tb_data)
                    processed_data.append(tb_data)

                combined_data = None
                combine_attempts = 3
                for i in range(combine_attempts):
                    try:
                        logger.debug("Combining data (attempt %s/%s)",
                                     (i+1), combine_attempts)
                        combined_data = self.combine_data(processed_data)
                        break  # Success
                    except Exception as combine_err:
                        logger.warning("Combine data failed on attempt %s: %s",
                                       (i+1), combine_err)

                        if i == combine_attempts - 1:  # Last attempt
                            raise combine_err  # Re-raise the exception
                        time.sleep(0.5)  # Short delay before retry

                # 1.3 Save to NetCDF
                output_filename = \
                    self.generate_output_filename(target_datetime)
                saved_file = \
                    self.save_to_netcdf(combined_data, output_filename, \
                                        target_datetime)

                # 1.4 Simple verification - check if file exists and has reasonable size
                if os.path.exists(saved_file) and \
                   os.path.getsize(saved_file) > 1000:  # At least 1KB
                    logger.info("Successfully processed: %s",saved_file)
                    return combined_data, saved_file

                raise Exception(f"File not properly saved: {saved_file}")

            except Exception as e:
                logger.error("Attempt %s failed: %s", (attempt+1), str(e))
                if attempt < max_retries:
                    logger.info("Retrying in 2 seconds...")
                    time.sleep(2)
                else:
                    logger.error("Failed after %s attempts: %s", \
                                 (max_retries+1), str(e))
                    raise
        return None, None

    def get_descending_orbit(self, ds):
        """Extract all segments where latitude is decreasing"""
        # Get latitude values
        logger.info('Process to get descending orbit for data from NOAA')
        if 'lat' in ds.keys():
            lat_values = ds['lat']
        else:
            lat_values = ds['latitude']
        # Handle 2D latitude arrays
        if hasattr(lat_values, 'ndim') and lat_values.ndim == 2:
            # Use the first row to find the pattern (all rows should be
            # the same)
            lat_row = lat_values[:, 0]  # first row
            lat_diff = np.diff(lat_row)
            is_decreasing = lat_diff < 0
            decreasing_mask = np.zeros(len(lat_row), dtype=bool)
            decreasing_indices = np.where(is_decreasing)[0]
            # Mark both start and end points
            decreasing_mask[decreasing_indices] = True  # Start points
            decreasing_mask[decreasing_indices + 1] = True  # End points
        else:
            # Original 1D case
            lat_diff = np.diff(lat_values)
            is_decreasing = lat_diff < 0
            decreasing_mask = np.zeros(len(lat_values), dtype=bool)
            decreasing_indices = np.where(is_decreasing)[0]
            decreasing_mask[decreasing_indices] = True
            decreasing_mask[decreasing_indices + 1] = True

        # Get final indices
        final_indices = np.where(decreasing_mask)[0]
        ## Apply to all arrays along the latitude axis
        ds_result = {}
        for var_name, var_array in ds.items():
            if hasattr(var_array, 'ndim') and var_array.ndim == 2 \
               and var_array.shape == lat_values.shape:
                ds_result[var_name] = var_array[final_indices, :]
            else:
                ds_result[var_name] = var_array
        return ds_result

    def check_available_data(self, target_datetime):
        """Check for available AMSR2 files within 6-hour window
        target_datetime: UTC time"""
        # Implementation to search discover for files
        start_time = target_datetime - \
            timedelta(hours=self.config.time_window_hours)
        amsr2_path_root = self.config.project_path / \
            self.config.amsr2_path
        logger.info('Search files in: %s', amsr2_path_root)

        # Current month path
        year_str = target_datetime.strftime('%Y')
        day_of_year = target_datetime.strftime('%j')
        month = target_datetime.strftime('%m')
        hour = target_datetime.strftime('%H')  # 24-hour format (00-23)
        amsr2_path = ''
        if hour == '00':
            # Go back one day using datetime arithmetic
            prev_datetime = target_datetime - timedelta(days=1)
            year_str = prev_datetime.strftime('%Y')
            month = prev_datetime.strftime('%m')
            day_of_year = prev_datetime.strftime('%j')

        if self.config.source == 'NOAA':
            amsr2_path = os.path.join(amsr2_path_root, year_str, day_of_year)
        elif self.config.source == 'JAXA':
            amsr2_path = os.path.join(amsr2_path_root, year_str, month)
        else:
            logger.error("Wrong Source provided (either NOAA or JAXA)")


        all_files = []

        if os.path.exists(amsr2_path):
            all_files.extend(glob.glob(os.path.join(amsr2_path, "*.h5")))
            logger.info("Found %s files", len(all_files))
        else:
            logger.warning("Path does not exist")

        file_list = self._get_file_list(all_files, start_time, \
                                        target_datetime)
        if not file_list:
            error_msg = (f"No AMSR2 descending files found for time window "
                         f"{start_time.strftime('%Y-%m-%d %H:%M')} to "
                         f"{target_datetime.strftime('%Y-%m-%d %H:%M')}")
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        txt = \
             f"Found {len(file_list)} AMSR2 descending files for time window " + \
             f"{start_time.strftime('%Y-%m-%d %H:%M')} to " + \
             f"{target_datetime.strftime('%Y-%m-%d %H:%M')}"
        logger.info(txt)
        for file_path in file_list:
            logger.info("   - %s", file_path)
        self.amsr2_files = file_list
        return file_list

    def _get_file_list(self, all_files, start_time, target_datetime):
        """Internal function to get file list."""
        file_list = []
        for file_path in all_files:
            try:
                filename = os.path.basename(file_path)
                if self.config.source == 'JAXA':  # "JAXA" or "NOAA":
                    orbit = filename.split('_')[2][3]
                    # Only read descending pass
                    # NOAA (sample): 20250120203202_GW1AM2_202501201817_126A_L1DLRTBR_2210210.h5
                    # NOAA (new): GW1AM2_202501201817_126A_L1DLRTBR_2210210.h5
                    # JAXA: GW1AM2_201904122301_172D_L1SGRTBR_2220220.h5
                    if orbit == 'D':
                        time_str = filename.split('_')[1]
                        if len(time_str) == 12:  # YYYYMMDDHHMM UTC time
                            file_datetime = datetime.strptime(time_str, '%Y%m%d%H%M')
                        else:
                            logger.warning('Unexpected time format in filename: %s', filename)
                            continue
                        if start_time <= file_datetime <= target_datetime:
                            file_list.append(file_path)
                elif self.config.source == 'NOAA':  # "JAXA" or "NOAA":
                    time_str = filename.split('_')[1][0:12]
                    if len(time_str) == 12:  # YYYYMMDD UTC time
                        file_datetime = datetime.strptime(time_str, '%Y%m%d%H%M')
                    else:
                        logger.warning('Unexpected time format in filename: %s', filename)
                        continue
                    if start_time <= file_datetime <= target_datetime:
                        file_list.append(file_path)

            except (IndexError, ValueError) as e:
                logger.warning("Could not parse datetime from filename %s: %s", \
                               filename, e)
                continue

        file_list.sort()
        return file_list

    def get_amsr_l1r(self, filename):
        """
        Read AMSR L1R data from HDF5 file.

        Parameters:
        -----------
        filename : str
            Path to the HDF5 file

        Returns:
        --------
        dict : Dictionary containing all AMSR data arrays and metadata
        """

        # Check if file exists
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Cannot find file {filename}")

        try:
            with h5py.File(filename, 'r') as file:
                # Define TB channels with their dataset names
                tb_channels = {
                    'tb_6v': "Brightness Temperature (res06,6.9GHz,V)",
                    'tb_6h': "Brightness Temperature (res06,6.9GHz,H)",
                    'tb_7v': "Brightness Temperature (res06,7.3GHz,V)",
                    'tb_7h': "Brightness Temperature (res06,7.3GHz,H)",
                    'tb_10v': "Brightness Temperature (res10,10.7GHz,V)",
                    'tb_10h': "Brightness Temperature (res10,10.7GHz,H)",
                    'tb_18v': "Brightness Temperature (res10,18.7GHz,V)",
                    'tb_18h': "Brightness Temperature (res10,18.7GHz,H)",
                    'tb_23v': "Brightness Temperature (res10,23.8GHz,V)",
                    'tb_23h': "Brightness Temperature (res10,23.8GHz,H)",
                    'tb_36v': "Brightness Temperature (res10,36.5GHz,V)",
                    'tb_36h': "Brightness Temperature (res10,36.5GHz,H)",
                    'tb_89v': "Brightness Temperature (res10,89.0GHz,V)",
                    'tb_89h': "Brightness Temperature (res10,89.0GHz,H)"
                }

                results = {}
                # Get 89GHz lat/lon data for dimensions
                results = self._get_amsr_l1r_step1(results, file)
                # Read TB data and get main dimensions
                results, tb_data = self._get_amsr_l1r_step2(results,
                                                            file, \
                                                            tb_channels)
                # Read time data
                results = self._get_amsr_l1r_step3(results, file)
                # Read land/water fraction
                results = self._get_amsr_l1r_step4(results, file)
                # Generate interpolated lat/lon grid
                results = self._get_amsr_l1r_step5(results)
                # Add RFI flags
                results = self._get_amsr_l1r_step6(results, file)
                # Generate flags
                results = self._get_amsr_l1r_step7(results, tb_data)
                return results
        except Exception as e:
            logger.error("Error reading HDF5 file: %s", e)
            raise

    def _get_amsr_l1r_step1(self, results, file):
        """Internal function for Step 1 of get_amsr_l1r"""
        # STEP 1: Get 89GHz lat/lon data for dimensions
        try:
            lat89_data = file["Latitude of Observation Point for 89A"][:]
            lon89_data = file["Longitude of Observation Point for 89A"][:]
            n89, m89 = lat89_data.shape
            results['lat89'] = lat89_data.astype(np.float32)
            results['lon89'] = lon89_data.astype(np.float32)
            results['n89'] = n89
            results['m89'] = m89
            logger.info("89GHz lat/lon dimensions: %s x %s", \
                        n89, m89)
        except KeyError as e:
            raise ValueError(f"Failed to read 89GHz coordinate data: {e}") from e
        return results

    def _get_amsr_l1r_step2(self, results, file, tb_channels):
        """Internal method to handle Step 2 of reading AMSR2 data"""
        # STEP 2: Read TB data and get main dimensions
        tb_data = {}
        n, m = None, None
        for tb_name, dataset_name in tb_channels.items():
            try:
                if dataset_name in file:
                    # Read and scale TB data (original is uint16 with
                    # scale factor 0.01)
                    raw_data = file[dataset_name][:]
                    tb_data[tb_name] = (raw_data * 0.01).astype(np.float32)
                    # Get dimensions from first successful read
                    if n is None or m is None:
                        n, m = tb_data[tb_name].shape
                        logger.info("TB data dimensions: %s x %s", n, m)
                    logger.info("Successfully read %s", tb_name)
                else:
                    logger.warning("dataset not found: %s", dataset_name)
                    tb_data[tb_name] = \
                        np.zeros((n, m), \
                                 dtype=np.float32) if n and m else None
            except Exception as e:
                logger.info("[WARN] Failed to read %s: %s", tb_name, e)
                tb_data[tb_name] = \
                    np.zeros((n, m), dtype=np.float32) if n and m else None
        # Add TB data to results
        results.update(tb_data)
        results['n'] = n
        results['m'] = m
        return results, tb_data

    def _get_amsr_l1r_step3(self, results, file):
        """Internal method for Step 3 of reading AMSR2 data"""
        # STEP 3: Read time data
        try:
            if "Scan Time" in file:
                tb_time_seconds = file["Scan Time"][:].astype(np.float64)
                results['tb_time_seconds'] = tb_time_seconds
                logger.info("Successfully read scan time data")
            else:
                logger.warning("Scan Time dataset not found")
                results['tb_time_seconds'] = \
                    np.zeros(results['m'],
                             dtype=np.float64)
        except Exception as e:
            logger.warning("Failed to read scan time: %s", e)
            results['tb_time_seconds'] = \
                np.zeros(results['m'], dtype=np.float64)
        return results

    def _get_amsr_l1r_step4(self, results, file):
        """Internal method for Step 3 of reading AMSR2 data"""
        # STEP 4: Read land/water fraction
        try:
            if "Land_Ocean Flag 6 to 36" in file:
                # Read 3D array and extract layer 4 (36GHz channel)
                land_ocean_3d = file["Land_Ocean Flag 6 to 36"][:]
                if land_ocean_3d.ndim == 3:
                    land_water_frac = \
                        land_ocean_3d[3, :, :].astype(np.int32)  # Layer 4 (0-indexed)
                else:
                    land_water_frac = land_ocean_3d.astype(np.int32)
                results['land_water_frac'] = land_water_frac
                logger.info("Successfully read land/water fraction")
            else:
                logger.warning("Land_Ocean Flag dataset not found")
                results['land_water_frac'] = \
                    np.zeros((results['n'], results['m']),
                             dtype=np.int32)
        except Exception as e:
            logger.warning("Failed to read land/water fraction: %s", e)
            results['land_water_frac'] = \
                np.zeros((results['n'], results['m']), \
                         dtype=np.int32)
        return results

    def _get_amsr_l1r_step5(self, results):
        """Internal method for Step 5 of reading AMSR2 data."""
        # STEP 5: Generate interpolated lat/lon grid
        try:
            # Interpolate 89GHz lat/lon to TB grid dimensions
            lat_interp, lon_interp = self.resample_and_grid(
                results['lat89'], results['lon89'],
                results['n'], results['m'])
            results['lat'] = lat_interp
            results['lon'] = lon_interp
            logger.info("Successfully interpolated lat/lon coordinates")
        except Exception as e:
            logger.warning("Failed to interpolate coordinates: %s", e)
            results['lat'] = \
                np.zeros((results['n'], results['m']),
                         dtype=np.float32)
            results['lon'] = \
                np.zeros((results['n'], results['m']),
                         dtype=np.float32)
        return results

    def _get_amsr_l1r_step6(self, results, file):
        """Internal method for Step 6 of reading AMSR2 data."""
        # STEP 6: add RFI flags
        layer = 'Pixel Data Quality 6 to 36'
        if layer in file:
            qflag = file[layer][:]
            # Stack along first axis (axis=0) to create channels dimension first
            rfi_data = np.stack([
                (qflag >> i & 3) > 2 for i in [0, 2, 4, 6]], axis=0)
            # Create dimensions with channels first
            dims = ['channels'] + \
                [f'phony_dim_{i}' for i in range(qflag.ndim)]
            rfi = xr.DataArray(
                rfi_data,
                dims=dims,
                coords={'channels': ['6v', '6h', '7v', '7h']}
            )
            rfi_v = rfi.sel(channels='7v') | rfi.sel(channels='6v')
            rfi_h = rfi.sel(channels='7h') | rfi.sel(channels='6h')
            # Input arrays are 2035x486, need to resample to 2035x243
            if 'phony_dim_1' in rfi_h.dims:
                rfi_h = \
                    rfi_h.coarsen(phony_dim_1=2).any(). \
                    rename({'phony_dim_1': 'phony_dim_3'})
                rfi_v = rfi_v.coarsen(phony_dim_1=2).any(). \
                    rename({'phony_dim_1': 'phony_dim_3'})
            results['rfi_h'] = np.array(rfi_h)
            results['rfi_v'] = np.array(rfi_v)
        return results

    def _get_amsr_l1r_step7(self, results, tb_data):
        """Internal method for Step 7 of reading AMSR2 data."""
        # STEP 7: Generate flags
        rain, cold_deserts, frozen_ground, glacier, \
            snow, precip = self.apply_flags(
                tb_data, results['land_water_frac'], results['n'],
                results['m'])
        # Bit value for flag
        # 1 - rain - NOAA
        # 2 - cold deserts
        # 3 - frozen ground
        # 4 - glacier
        # 5 - snow cover - Rajat
        # 6 - precip - Rajat, rain
        # 7 - rfi_h
        # 8 - rfi_v

        results['rain'] = rain
        results['cold_deserts'] = cold_deserts
        results['frozen_ground'] = frozen_ground
        results['glacier'] = glacier
        results['snow'] = snow
        results['precip'] = precip
        logger.info("Successfully completed AMSR L1R data reading")
        return results

    def resample_and_grid(self, lat89, lon89, n_out, m_out):
        """
        Bilinear interpolation of 2D lat/lon arrays to new dimensions.

        Parameters:
        -----------
        lat89, lon89 : array_like
            Input coordinate arrays
        n_out, m_out : int
            Output dimensions

        Returns:
        --------
        tuple : (lat_interp, lon_interp) interpolated arrays
        """
        n_in, m_in = lat89.shape

        # Create coordinate arrays for interpolation
        x_in = np.arange(n_in)
        y_in = np.arange(m_in)
        x_out = np.linspace(0, n_in - 1, n_out)
        y_out = np.linspace(0, m_in - 1, m_out)

        # Create interpolation functions - same as bilinear when kx=1 ky=1
        lat_interp_func = RectBivariateSpline(x_in, y_in, lat89, \
                                              kx=1, ky=1)
        lon_interp_func = RectBivariateSpline(x_in, y_in, lon89, \
                                              kx=1, ky=1)

        # Interpolate to new grid
        lat_interp = lat_interp_func(x_out, y_out).astype(np.float32)
        lon_interp = lon_interp_func(x_out, y_out).astype(np.float32)
        return lat_interp, lon_interp

    def apply_flags(self, tb_data, land_water_frac, n, m):
        """
        Generate snow and precipitation flags based on TB values.

        Parameters:
        -----------
        tb_data : dict
            Dictionary containing TB arrays
        land_water_frac : array_like
            Land/water fraction array
        n, m : int
            Array dimensions

        Returns:
        --------
        tuple : (rain, cold_deserts, frozen_ground, glacier) flag arrays
        """

        rain = np.zeros((n, m), dtype=np.int32)
        cold_deserts = np.zeros((n, m), dtype=np.int32)
        frozen_ground = np.zeros((n, m), dtype=np.int32)
        glacier = np.zeros((n, m), dtype=np.int32)
        snow = np.zeros((n, m), dtype=np.int32)
        precip = np.zeros((n, m), dtype=np.int32)

        # Get required TB channels
        tb_18v = tb_data.get('tb_18v', np.zeros((n, m)))
        tb_18h = tb_data.get('tb_18h', np.zeros((n, m)))
        tb_23v = tb_data.get('tb_23v', np.zeros((n, m)))
        tb_36v = tb_data.get('tb_36v', np.zeros((n, m)))
        tb_89v = tb_data.get('tb_89v', np.zeros((n, m)))

        logger.info("Generating rain, cold deserts, frozen ground,"
                       " and glacier flags using NOAA method;"
                       "Generating snow and precip flags using Rajat method ...")

        for i in range(n):
            for j in range(m):
                # Only process data over land (land_water_frac >= 50)
                if land_water_frac[i, j] >= 50:
                    # Check if all required TB values are valid (> 0)
                    if (tb_18v[i, j] > 0 and tb_18h[i, j] > 0 and
                            tb_23v[i, j] > 0 and tb_89v[i, j] > 0 and
                            tb_36v[i, j] > 0):

                        # Apply NOAA's method in flag rain, cold deserts,
                        # frozen ground and glacier (Grody's SCA algorithm)
                        scat = tb_23v[i, j] - tb_89v[i, j]
                        sc37 = tb_18v[i, j] - tb_36v[i, j]
                        pd19 = tb_18v[i, j] - tb_18h[i, j]  # tt18
                        scx = tb_36v[i, j] - tb_89v[i, j]
                        scat = max(scat, sc37)
                        tt = (165 + 0.49 * tb_89v[i, j])
                        if ((tb_23v[i, j] >= 254.0 and scat <= 2.0) or
                                tb_23v[i, j] >= 258.0 or tb_23v[i, j] >= tt):
                            rain[i, j] = 1
                        else:
                            rain[i, j] = 0
                        if pd19 >= 18.0 and sc37 <= 10.0 and scx <= 10.0:
                            cold_deserts[i, j] = 1
                        else:
                            cold_deserts[i, j] = 0
                        if scat <= 6.0 and pd19 >= 8.0:
                            frozen_ground[i, j] = 1
                        else:
                            frozen_ground[i, j] = 0
                        if tb_23v[i, j] <= 210.0 or \
                           (tb_23v[i, j] <= 229.0 and pd19 >= 23.0):
                            glacier[i, j] = 1
                        else:
                            glacier[i, j] = 0
                        # Ehsan's method in classify snow and precipitation
                        # Calculate scattering index and polarization difference
                        sil = (451.88 - 0.44 * tb_18v[i, j] - 1.775 * tb_23v[i, j] +
                               0.00574 * tb_23v[i, j] ** 2 - tb_89v[i, j])
                        tt18 = tb_18v[i, j] - tb_18h[i, j]
                        if sil > 10:
                            if (tb_23v[i, j] <= 264.0 and
                                    tb_23v[i, j] <= (175.0 + 0.49 * tb_89v[i, j])):
                                # Snow branch
                                snow[i, j] = 1
                                # Additional snow checks
                                if (tt18 >= 18 and
                                        (tb_18v[i, j] - tb_36v[i, j]) <= 10 and
                                        (tb_36v[i, j] - tb_89v[i, j]) <= 10):
                                    snow[i, j] = 0
                                if (tt18 >= 8 and
                                        (tb_18v[i, j] - tb_36v[i, j]) <= 2 and
                                        (tb_23v[i, j] - tb_89v[i, j]) <= 6):
                                    snow[i, j] = 1
                            else:
                                # Precipitation branch
                                snow[i, j] = 0
                                precip[i, j] = 1
                                # Additional precip checks
                                if tt18 > 20:
                                    precip[i, j] = 0
                                if tb_89v[i, j] > 253 and tt18 > 7:
                                    precip[i, j] = 0

        logger.info("Rain, cold_deserts, frozen_ground, "
                       "glacier, snow, and precip flags generated")
        return rain, cold_deserts, frozen_ground, glacier, snow, precip

    def combine_data(self, processed_data):
        """
        Combine multiple AMSR2 L1R datasets with potentially different extents.

        Parameters:
        -----------
        processed_data : list
            List of dictionaries, each containing AMSR2 data from get_amsr_l1r()

        Returns:
        --------
        dict : Combined dataset on a common grid
        """
        if not processed_data:
            logger.warning("No processed data to combine")
            return {}
        if len(processed_data) == 1:
            logger.info("Only one dataset found, no combining needed")
            combined_result_flag = self._merge_flag(processed_data[0])
            return combined_result_flag
        logger.info("Combining %s AMSR2 datasets", len(processed_data))

        # Step 1: Determine the common geographic extent
        combined_extent = self._get_combined_extent(processed_data)

        # Step 2: Create target grid
        target_grid = self._create_target_grid(combined_extent)

        # Step 3: Regrid and combine all datasets
        combined_result = self._merge_datasets(processed_data, target_grid)

        # step4: merge flag into one variable
        combined_result_flag = self._merge_flag(combined_result)

        logger.info("Successfully combined AMSR2 datasets")
        return combined_result_flag

    def _merge_flag(self, results):
        """
        Combine multiple data quality flag to one variable.

        Parameters:
        -----------
        results : dict
            Dictionary containing AMSR2 data from get_amsr_l1r()

        Returns:
        --------
        dict : Combined dataset on a common grid
        """
        logger.info("Combine multiple data quality flag to pixel_qual_flag")

        pixel_qual_flag = np.zeros_like(results['rain'], dtype=np.uint8)
        # Set bits for each condition (1 means condition is present)
        pixel_qual_flag += (results['rain'] > 0).astype(np.uint8) * 1  # Bit 0: rain
        pixel_qual_flag += (results['cold_deserts'] > 0).astype(np.uint8) * 2  # Bit 1: cold deserts
        pixel_qual_flag += (results['frozen_ground'] > 0). \
            astype(np.uint8) * 4  # Bit 2: frozen ground
        pixel_qual_flag += (results['glacier'] > 0).astype(np.uint8) * 8  # Bit 3: glacier
        pixel_qual_flag += (results['snow'] > 0).astype(np.uint8) * 16  # Bit 4: snow cover
        pixel_qual_flag += (results['precip'] > 0).astype(np.uint8) * 32  # Bit 5: precip
        pixel_qual_flag += (results['rfi_h']).astype(np.uint8) * 64  # Bit 6: rfi_h
        pixel_qual_flag += (results['rfi_v']).astype(np.uint8) * 128  # Bit 7: rfi_v

        # Add to results
        results['pixel_qual_flag'] = pixel_qual_flag

        # Add metadata
        if hasattr(results['pixel_qual_flag'], 'attrs'):
            results['pixel_qual_flag'].attrs[
                'flag_meanings'] = \
                    'rain cold_deserts frozen_ground glacier snow_cover precip rfi_h rfi_v'
            results['pixel_qual_flag'].attrs['flag_values'] = [1, 2, 4, 8, 16, 32, 64, 128]
            results['pixel_qual_flag'].attrs['flag_descriptions'] = {
                1: 'rain - NOAA',
                2: 'cold deserts',
                4: 'frozen ground',
                8: 'glacier',
                16: 'snow cover - Rajat',
                32: 'precip - Rajat, rain',
                64: 'rfi_h',
                128: 'rfi_v'
            }

        return results

    def _get_combined_extent(self, processed_data):
        """Determine the geographic extent that encompasses all datasets"""
        lat_min, lat_max = float('inf'), float('-inf')
        lon_min, lon_max = float('inf'), float('-inf')
        for data in processed_data:
            if 'lat' in data and 'lon' in data:
                lat = data['lat']
                lon = data['lon']
                # Handle valid data only (mask out zeros/invalid values)
                valid_mask = (lat != 0) & (lon != 0) & np.isfinite(lat) & np.isfinite(lon)
                if np.any(valid_mask):
                    lat_valid = lat[valid_mask]
                    lon_valid = lon[valid_mask]
                    lat_min = min(lat_min, np.min(lat_valid))
                    lat_max = max(lat_max, np.max(lat_valid))
                    lon_min = min(lon_min, np.min(lon_valid))
                    lon_max = max(lon_max, np.max(lon_valid))

        extent = {
            'lat_min': lat_min,
            'lat_max': lat_max,
            'lon_min': lon_min,
            'lon_max': lon_max
        }
        txt = f"Combined extent: lat [{extent['lat_min']:.2f}, {extent['lat_max']:.2f}], " + \
            f"lon [{extent['lon_min']:.2f}, {extent['lon_max']:.2f}]"
        logger.info(txt)
        return extent

    def _create_target_grid(self, extent):
        """Create target latitude/longitude grid"""
        # Calculate grid dimensions based on target resolution
        n_lat = int(np.ceil((extent['lat_max'] - extent['lat_min']) / \
                            self.target_resolution))
        n_lon = int(np.ceil((extent['lon_max'] - extent['lon_min']) / \
                            self.target_resolution))
        # Create coordinate arrays
        lat_1d = np.linspace(extent['lat_min'], extent['lat_max'], n_lat)
        lon_1d = np.linspace(extent['lon_min'], extent['lon_max'], n_lon)
        # Create 2D grids
        target_lon, target_lat = np.meshgrid(lon_1d, lat_1d)
        logger.info("Target grid dimensions: %s x %s", n_lat, n_lon)
        return {
            'lat': target_lat.astype(np.float32),
            'lon': target_lon.astype(np.float32),
            'n': n_lat,
            'm': n_lon
        }

    def _merge_datasets(self, processed_data, target_grid):
        """
        Regrid and merge all datasets onto common grid
        """
        target_shape = (target_grid['n'], target_grid['m'])
        # Initialize combined result with target grid
        combined_result = {
            'lat': target_grid['lat'],
            'lon': target_grid['lon'],
            'n': target_grid['n'],
            'm': target_grid['m']
        }
        # Define TB channels and other variables to combine
        tb_channels = ['tb_6v', 'tb_6h', 'tb_7v', 'tb_7h', 'tb_10v', 'tb_10h',
                       'tb_18v', 'tb_18h', 'tb_23v', 'tb_23h', 'tb_36v', 'tb_36h',
                       'tb_89v', 'tb_89h']
        other_vars = ['rain', 'cold_deserts', 'frozen_ground', 'glacier', 'snow', 'precip',
                      'land_water_frac', 'pixel_qual_flag', 'rfi_h', 'rfi_v']
        # other_vars = ['pixel_qual_flag']
        all_vars = tb_channels + other_vars
        # Initialize combined arrays
        for var in all_vars:
            combined_result[var] = np.full(target_shape, np.nan, dtype=np.float32)
        # Initialize weight array to track data coverage
        weight_total = np.zeros(target_shape, dtype=np.float32)
        # Process each dataset
        for idx, data in enumerate(processed_data):
            logger.info("Processing dataset %s/%s", idx+1, len(processed_data))
            # Get source coordinates
            src_lat = data['lat'].flatten()
            src_lon = data['lon'].flatten()
            # Create mask for valid coordinates
            valid_mask = (src_lat != 0) & (src_lon != 0) \
                & np.isfinite(src_lat) & np.isfinite(src_lon)
            if not np.any(valid_mask):
                logger.warning("No valid coordinates in dataset %s, skipping", idx+1)
                continue
            # Get valid source points
            src_points = np.column_stack([src_lat[valid_mask], src_lon[valid_mask]])
            # Create target points
            target_lat_flat = target_grid['lat'].flatten()
            target_lon_flat = target_grid['lon'].flatten()
            target_points = np.column_stack([target_lat_flat, target_lon_flat])
            # Build KDTree for nearest neighbor interpolation
            tree = cKDTree(src_points)
            distances, indices = tree.query(target_points, k=1)
            # Set maximum distance threshold (in degrees)
            max_distance = self.target_resolution * 2  # Allow up to 2 grid cells distance
            valid_interp = distances < max_distance
            # Map valid source indices back to original array indices
            valid_src_indices = np.where(valid_mask)[0]
            # Process each variable
            for var in all_vars:
                if var in data:
                    src_data = data[var].flatten()
                    # Initialize target array for this variable
                    target_var = np.full(target_shape[0] * \
                                         target_shape[1], np.nan)
                    # Interpolate valid points
                    valid_targets = valid_interp & \
                        np.isfinite(src_data[valid_src_indices[indices]])
                    if np.any(valid_targets):
                        target_var[valid_targets] = \
                            src_data[valid_src_indices[indices[valid_targets]]]
                    # Reshape and combine with existing data
                    target_var_2d = target_var.reshape(target_shape)
                    # Use weighted average where data overlaps
                    current_valid = np.isfinite(target_var_2d)
                    existing_valid = np.isfinite(combined_result[var])
                    # Where no existing data, use new data
                    no_existing = ~existing_valid & current_valid
                    combined_result[var][no_existing] = \
                        target_var_2d[no_existing]
                    # Where both exist, take average (could implement
                    # more sophisticated weighting)
                    both_valid = existing_valid & current_valid
                    if np.any(both_valid):
                        # Simple average - you might want to implement
                        # distance-based weighting
                        combined_result[var][both_valid] = \
                            (combined_result[var][both_valid] + \
                             target_var_2d[both_valid]) / 2
            # Update weight tracking
            valid_points_2d = valid_interp.reshape(target_shape)
            weight_total += valid_points_2d.astype(np.float32)
        # Convert NaN back to appropriate fill values for integer arrays
        for var in other_vars:
            if var in combined_result:
                nan_mask = np.isnan(combined_result[var])
                combined_result[var][nan_mask] = 0
                combined_result[var] = combined_result[var].astype(np.int32)
        # Log coverage statistics
        total_pixels = target_shape[0] * target_shape[1]
        covered_pixels = np.sum(weight_total > 0)
        coverage_percent = (covered_pixels / total_pixels) * 100
        txt = \
            f"Combined grid coverage: {covered_pixels}/{total_pixels}"
        txt += \
            f" pixels ({coverage_percent:.1f}%)"
        logger.info(txt)
        return combined_result

    def save_to_netcdf(self, combined_data, output_path, target_datetime):
        """
        Save combined AMSR2 dataset to NetCDF file with proper metadata
        and CF conventions.
        """
        try:
            logger.info("Saving combined dataset to %s", output_path)
            # Create coordinate arrays
            n_lat, n_lon = combined_data['n'], combined_data['m']
            # Extract 1D coordinate arrays from 2D grids
            lat_1d = combined_data['lat'][:, 0]  # First column
            lon_1d = combined_data['lon'][0, :]  # First row
            # Create time coordinate - use pandas datetime for better
            # compatibility
            time_coord = pd.to_datetime([target_datetime])
            # Define data variables with attributes
            data_vars = {}
            # TB channels
            tb_channels = {
                'tb_6v': {'long_name': 'Brightness Temperature 6.9GHz V-pol', 'units': 'K'},
                'tb_6h': {'long_name': 'Brightness Temperature 6.9GHz H-pol', 'units': 'K'},
                'tb_7v': {'long_name': 'Brightness Temperature 7.3GHz V-pol', 'units': 'K'},
                'tb_7h': {'long_name': 'Brightness Temperature 7.3GHz H-pol', 'units': 'K'},
                'tb_10v': {'long_name': 'Brightness Temperature 10.7GHz V-pol', 'units': 'K'},
                'tb_10h': {'long_name': 'Brightness Temperature 10.7GHz H-pol', 'units': 'K'},
                'tb_18v': {'long_name': 'Brightness Temperature 18.7GHz V-pol', 'units': 'K'},
                'tb_18h': {'long_name': 'Brightness Temperature 18.7GHz H-pol', 'units': 'K'},
                'tb_23v': {'long_name': 'Brightness Temperature 23.8GHz V-pol', 'units': 'K'},
                'tb_23h': {'long_name': 'Brightness Temperature 23.8GHz H-pol', 'units': 'K'},
                'tb_36v': {'long_name': 'Brightness Temperature 36.5GHz V-pol', 'units': 'K'},
                'tb_36h': {'long_name': 'Brightness Temperature 36.5GHz H-pol', 'units': 'K'},
                'tb_89v': {'long_name': 'Brightness Temperature 89.0GHz V-pol', 'units': 'K'},
                'tb_89h': {'long_name': 'Brightness Temperature 89.0GHz H-pol', 'units': 'K'}
            }
            # Add TB variables
            for var_name, attrs in tb_channels.items():
                if var_name in combined_data:
                    data_vars[var_name] = (
                        ['time', 'lat', 'lon'],
                        combined_data[var_name][np.newaxis, :, :],  # Add time dimension
                        {
                            'long_name': attrs['long_name'],
                            'units': attrs['units'],
                            'valid_range': [0., 350.],
                            '_FillValue': np.nan,
                            'grid_mapping': 'crs'
                        }
                    )
            # Add auxiliary variables
            aux_vars = {
                'land_water_frac': {
                    'long_name': 'Land Water Fraction',
                    'units': 'percent',
                    'valid_range': [0, 100],
                    'dtype': np.int32
                },
                'pixel_qual_flag': {
                    'long_name': 'Pixel Quality Flag',
                    'units': 'dimensionless',
                    'dtype': np.uint8
                }
            }
            for var_name, attrs in aux_vars.items():
                if var_name in combined_data:
                    fill_val = 0 if 'flag' in var_name or \
                        var_name == 'land_water_frac' else -9999
                    data_vars[var_name] = (
                        ['time', 'lat', 'lon'],
                        combined_data[var_name][np.newaxis, :, :]. \
                        astype(attrs['dtype']),
                        {
                            'long_name': attrs['long_name'],
                            'units': attrs['units'],
                            '_FillValue': fill_val,
                            'grid_mapping': 'crs',
                            **{k: v for k, v in attrs.items() \
                               if k not in ['long_name', 'units', 'dtype']}
                        }
                    )
            # Add coordinate reference system
            data_vars['crs'] = (
                [],
                0,  # Scalar variable
                {
                    'grid_mapping_name': 'latitude_longitude',
                    'longitude_of_prime_meridian': 0.0,
                    'semi_major_axis': 6378137.0,
                    'inverse_flattening': 298.257223563,
                    'spatial_ref': \
                    'GEOGCS["WGS 84",DATUM["WGS_1984",' + \
                    'SPHEROID["WGS 84",6378137,298.257223563]],' + \
                    'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
                }
            )
            # Create coordinates dictionary - REMOVE units and calendar from time
            coords = {
                'time': (
                    ['time'],
                    time_coord,
                    {
                        'long_name': 'Time',
                        'standard_name': 'time'
                        # Remove both 'units' and 'calendar' -
                        # xarray will handle these automatically
                    }
                ),
                'lat': (
                    ['lat'],
                    lat_1d,
                    {
                        'long_name': 'Latitude',
                        'standard_name': 'latitude',
                        'units': 'degrees_north',
                        'valid_range': [-90., 90.],
                        'axis': 'Y'
                    }
                ),
                'lon': (
                    ['lon'],
                    lon_1d,
                    {
                        'long_name': 'Longitude',
                        'standard_name': 'longitude',
                        'units': 'degrees_east',
                        'valid_range': [-180., 180.],
                        'axis': 'X'
                    }
                )
            }
            # Create xarray Dataset
            ds = xr.Dataset(
                data_vars=data_vars,
                coords=coords,
                attrs={
                    'title': 'AMSR2 L1R Combined Brightness Temperature Data',
                    'summary': \
                    'Combined AMSR2 Level 1R brightness temperature data regridded to common grid',
                    'source': 'AMSR2 Level 1R',
                    'institution': 'NASA GSFC',
                    'creator_name': 'AMSR2 Data Processor',
                    'date_created': datetime.now().strftime('%Y-%m-%dT%H:%M:%SZ'),
                    'target_datetime': target_datetime.strftime('%Y-%m-%dT%H:%M:%SZ'),
                    'spatial_resolution': f'{self.target_resolution} degrees',
                    'geospatial_lat_min': float(np.min(lat_1d)),
                    'geospatial_lat_max': float(np.max(lat_1d)),
                    'geospatial_lon_min': float(np.min(lon_1d)),
                    'geospatial_lon_max': float(np.max(lon_1d)),
                    'geospatial_lat_units': 'degrees_north',
                    'geospatial_lon_units': 'degrees_east',
                    'time_coverage_start': target_datetime.strftime('%Y-%m-%dT%H:%M:%SZ'),
                    'time_coverage_end': target_datetime.strftime('%Y-%m-%dT%H:%M:%SZ'),
                    'Conventions': 'CF-1.8',
                    'processing_info': f'Combined from {len(self.amsr2_files)} AMSR2 L1R files',
                    'input_files': '; '.join([os.path.basename(f) for f in self.amsr2_files])
                }
            )
            # Set encoding for compression and chunking
            encoding = {}
            for var in ds.data_vars:
                if var not in ['crs', 'pixel_qual_flag']:
                    encoding[var] = {
                        'zlib': True,
                        'complevel': 6,
                        'shuffle': True,
                        'chunksizes': (1, min(256, n_lat), \
                                       min(256, n_lon))
                    }
            encoding['pixel_qual_flag'] = {
                'dtype': 'uint8',
                'zlib': True,
                'complevel': 6,
                'shuffle': True,
                'chunksizes': (1, min(256, n_lat), min(256, n_lon)),
                '_FillValue': None
            }
            # Add time encoding
            encoding['time'] = {
                'units': f'seconds since {target_datetime.strftime("%Y-%m-%d %H:%M:%S")}',
                'calendar': 'gregorian'
            }
            # Create output directory if it doesn't exist
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            # Save to NetCDF
            ds.to_netcdf(
                output_path,
                encoding=encoding,
                format='NETCDF4'
            )
            # Log file info
            file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
            logger.info("Successfully saved NetCDF file: %s", \
                        output_path)
            txt = f"File size: {file_size_mb:.2f} MB"
            logger.info(txt)
            logger.info("Grid dimensions: %s x %s", n_lat, n_lon)
            txt = f"Variables saved: {list(data_vars.keys())}"
            logger.info(txt)
            return output_path

        except Exception as e:
            logger.error("Error saving NetCDF file: %s", e)
            raise

    def generate_output_filename(self, target_datetime, output_dir=None):
        """
        Generate standardized output filename for NetCDF files.

        Parameters:
        -----------
        target_datetime : datetime
            Target datetime for the dataset
        output_dir : str, optional
            Output directory path

        Returns:
        --------
        str : Full path to output file
        """

        if output_dir is None:
            output_dir = getattr(self.config, 'amsr2_merge_path', \
                                 './amsr2_merge_path')

        # Add project_path to output_dir
        full_output_dir = self.config.project_path / output_dir

        # Create filename with timestamp
        timestamp_str = target_datetime.strftime('%Y%m%d%H%M')
        filename = f'AMSR2_L1R_combined_{timestamp_str}.nc'
        return full_output_dir / filename
