#!/usr/bin/python3
"""
Create GEOS-S2S-V3 Downward Radiation Components

Workflow
--------
1. Prepare files and directories (input and output)
2. Read surface variables
3. Compute emissivity from vegtype, leaf area index, and fractional snow cover
4. Compute downward longwave radiation from emissivity, surface temperature, 
   and fractional land cover
5. Read albedo variables
6. Compute downward shortwave radiation from albedo, net shortwave radiation, 
   and fractional land cover

Adapted from: LWSWpostcal.py
Author: Ryan Zamora; LWSWpostcal.py written by the GMAO Team
"""
#
# Modules
#
import sys
import os
import calendar
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
from netCDF4 import Dataset # pylint: disable=no-name-in-module
from scipy.ndimage import distance_transform_edt

# Custom Modules
from bcsd_geosv3_functions import (generate_geosv3_3hr_surface_filename,
                                   generate_geosv3_3hr_albedo_filename,
                                   generate_geosv3_3hr_radiation_filename)

def fill_radiation_gaps_nearest_neighbor(geos_mask, ldt_mask, lw_data, sw_data):
    """
    Fill gaps in radiation data where LDT has land but GEOS mask has water/NaN.
    Uses nearest neighbor filling.
    Parameters:
    -----------
    lw_data : numpy array (time, lat, lon)
        Longwave radiation data on GEOS mask grid
    sw_data : numpy array (time, lat, lon)
        Shortwave radiation data on GEOS mask grid
    Returns:
    --------
    lw_filled : numpy array (time, lat, lon)
        Filled longwave data
    sw_filled : numpy array (time, lat, lon)
        Filled shortwave data
    """
    # Identify cells that need filling
    ldt_land = ldt_mask == 1
    geos_not_land = (geos_mask != 1) | np.isnan(geos_mask)
    cells_need_filling = ldt_land & geos_not_land

    num_to_fill = np.sum(cells_need_filling)
    print(f"Cells needing gap-fill: {num_to_fill:,} ({100*num_to_fill/ldt_mask.size:.2f}%)")

    if num_to_fill == 0:
        print("No gaps to fill!")
        return lw_data, sw_data

    # Create mask of valid GEOS data
    geos_valid = (geos_mask == 1) & ~np.isnan(geos_mask)

    # Function to fill a single 2D slice using nearest neighbor
    def fill_slice_nearest(data_slice, valid_mask, fill_mask):
        """Fill using scipy's distance transform"""
        filled = data_slice.copy()

        if np.sum(fill_mask) > 0 and np.sum(valid_mask) > 0:
            indices = distance_transform_edt(
                ~valid_mask,
                return_distances=False,
                return_indices=True
            )

            filled[fill_mask] = data_slice[tuple(indices[:, fill_mask])]

        return filled

    # Process each time step
    print("Filling longwave radiation...")
    lw_filled = np.zeros_like(lw_data)
    for t in range(lw_data.shape[0]):
        if t % 50 == 0:
            print(f"  Processing time step {t}/{lw_data.shape[0]}")
        lw_filled[t, :, :] = fill_slice_nearest(
            lw_data[t, :, :],
            geos_valid,
            cells_need_filling
        )

    print("Filling shortwave radiation...")
    sw_filled = np.zeros_like(sw_data)
    for t in range(sw_data.shape[0]):
        if t % 50 == 0:
            print(f"  Processing time step {t}/{sw_data.shape[0]}")
        sw_filled[t, :, :] = fill_slice_nearest(
            sw_data[t, :, :],
            geos_valid,
            cells_need_filling
        )

    # Mask to LDT land only
    lw_filled[:, ~ldt_land] = np.nan
    sw_filled[:, ~ldt_land] = np.nan

    print("✓ Gap filling complete")

    return lw_filled, sw_filled

# Name for downward shortwave and longwave radiation (Names follow GEOSv3 convention)
longwave_down_name = 'LWS'
shortwave_down_name = 'SWGDWN'

# Coordinates
nlat, nlon = 361, 720
latitudes = np.linspace(-90, 90, nlat)
longitudes = np.linspace(-180, 179.5, nlon)

# Constants
SIGMA = 5.67e-8
NMONTHS = 3
EMSVEG = np.array([0.99560, 0.99000, 0.99560, 0.99320, 0.99280, 0.99180])
EMSBARESOIL = 0.94120
EMSSNO = 0.99999

#
# Functions
#
def _usage():
    """Print command line usage"""
    txt = (f"[INFO] Usage: {(sys.argv[0])} dir_geos dir_rad fcst_init_year"
           "fcst_init_month ensemble_number")
    print(txt)
    print("[INFO] where")
    print('[INFO] dir_geos: Directory where the GEOS files are located (string)')
    print('[INFO] dir_rad: Directory where the radiation files will be written (string)')
    print('[INFO] fcst_init_year: Initial forecast year (integer)')
    print('[INFO] fcst_init_month: Initial forecast month (integer)')
    print('[INFO] ensemble_number: Ensemble number (integer)')

def compute_land_emissivity(vegtype, lai, asnow):
    """Compute emissivity from vegetation type, leaf area index, and fractional snow cover"""
    emis_veg = EMSVEG[vegtype.astype(int) - 1]
    emis_land = emis_veg + (EMSBARESOIL - emis_veg) * np.exp(-lai)
    emissivity = emis_land * (1 - asnow) + EMSSNO * asnow
    return emissivity

def compute_lw_upward(emissivity, surface_temp):
    """Compute upward longwave radiation from emissivity, surface temperature"""
    lw_upward = emissivity * SIGMA * surface_temp**4
    return lw_upward

def compute_lw_downward(lw_net, lw_upward):
    """Compute downward longwave radiation from upward and net longwave radiation"""
    lw_downward = lw_upward + lw_net
    return lw_downward

def compute_sw_downward(albedo, sw_net, frland):
    """Compute downward shortwave radiation from albedo, net shortwave radiation, and land cover"""
    threshold = 0.999
    mask = albedo > threshold
    if np.any(mask):
        sys.exit(f"Aborting: {np.sum(mask)} grid points have albedo > {threshold}")

    land_mask = frland > 0.1

    sw_net = np.where(mask, 0.0, sw_net)
    sw_downward = sw_net / (1 - albedo)
    sw_downward = np.where(~land_mask, np.nan, sw_downward)
    sw_downward = np.where(land_mask & np.isnan(sw_downward), 0.0, sw_downward)
    return sw_downward

def read_sfc(filename):
    """Reads the surface variables needed to calculate radiation components"""
    if not os.path.isfile(filename):
        return (None,) * 7
    try:
        with Dataset(filename, 'r') as nc:
            vegtype = np.where(nc.variables['VEGTYPE'][:, :, :] > 6, 0,
                               nc.variables['VEGTYPE'][:, :, :])
            return (
                vegtype,
                nc.variables['LAI'][:, :, :],
                nc.variables['ASNOW'][:, :, :],
                nc.variables['TS'][:, :, :],
                nc.variables['FRLAND'][:, :, :],
                nc.variables['LWLAND'][:, :, :],
                nc.variables['SWLAND'][:, :, :]
            )
    except:
        return (None,) * 7

def read_albedo(filename, is_feb_and_not_leap_year = False):
    """Reads the albedo variable needed to calculated radiation components"""
    if not os.path.isfile(filename):
        return None
    try:
        with Dataset(filename, 'r') as nc:
            # Feb albedo file has 29 days; Remove last day (8 time steps) if Feb and not a leap year
            if is_feb_and_not_leap_year:
                return nc.variables['ALBEDO'][:-8, 0, :, :]
            else:
                return nc.variables['ALBEDO'][:, 0, :, :]
    except:
        return None

def write_nc4_file(output_filename, lw_var_name, lw_data, lw_long_name, sw_var_name,
                   sw_data, sw_long_name, filename_sfc):
    """Writes out the components for Downward Shortwave Radiation & Downward Longwave Radiation"""
    try:
        with Dataset(filename_sfc, 'r') as sample_nc:
            time = np.arange(lw_data.shape[0])*180.
            latitudes = np.linspace(-90.0, 90.0, 361)
            longitudes = np.linspace(-180.0, 179.5, 720)
            lat_attrs = {a: sample_nc.variables['lat'].getncattr(a)
                         for a in sample_nc.variables['lat'].ncattrs()}
            lon_attrs = {a: sample_nc.variables['lon'].getncattr(a)
                         for a in sample_nc.variables['lon'].ncattrs()}
            time_attrs = {a: sample_nc.variables['time'].getncattr(a)
                          for a in sample_nc.variables['time'].ncattrs()}
        out_dir = os.path.dirname(output_filename)
        print(f"Writing NetCDF to directory: {out_dir}")
        print(f"Filename: {os.path.basename(output_filename)}")

        with Dataset(output_filename, 'w', format='NETCDF4') as nc:
            nc.createDimension('time', len(time))
            nc.createDimension('lat', len(latitudes))
            nc.createDimension('lon', len(longitudes))

            lat_var = nc.createVariable('lat', np.float64, ('lat',))
            lon_var = nc.createVariable('lon', np.float64, ('lon',))
            time_var = nc.createVariable('time', np.float64, ('time',))
            lw_var = nc.createVariable(lw_var_name, np.float32,
                                       ('time', 'lat', 'lon'), fill_value=1.0e20)
            sw_var = nc.createVariable(sw_var_name, np.float32,
                                       ('time', 'lat', 'lon'), fill_value=1.0e20)

            lat_var[:] = latitudes
            lon_var[:] = longitudes
            time_var[:] = time
            lw_var[:, :, :] = lw_data
            sw_var[:, :, :] = sw_data

            lw_var.units = 'W/m^2'
            lw_var.long_name = lw_long_name
            sw_var.units = 'W/m^2'
            sw_var.long_name = sw_long_name

            for attr, val in lat_attrs.items(): lat_var.setncattr(attr, val)
            for attr, val in lon_attrs.items(): lon_var.setncattr(attr, val)
            for attr, val in time_attrs.items(): time_var.setncattr(attr, val)

            nc.title = 'NetCDF file using for radiations'
            nc.Conventions = 'CF-1.6'
    except Exception as e:
        print(f"Failed to write {output_filename}: {e}")

def main():
    """Main driver"""

    # Check for correct number of input arguments
    if len(sys.argv) != 5:
        print('[ERR] Invalid number of command line arguments!')
        _usage()
        sys.exit(1)

    # Read input arguments
    dir_geos = sys.argv[1]
    dir_rad = sys.argv[2]
    sup_dir = sys.argv[3]
    fcst_init_year = int(sys.argv[4])
    ldt_mask_file = f'{sup_dir}/bcsd_fcst/LDT_mask/GEOSv3_5km_landmask.nc4'
    geos_mask_file = f'{sup_dir}/bcsd_fcst/GEOS_mask/GEOSv3_5km_landmask.nc4'
    print("Loading masks...")
    ldt_ds = xr.open_dataset(ldt_mask_file)
    geos_ds = xr.open_dataset(geos_mask_file)

    ldt_mask = ldt_ds['LANDMASK'].values
    geos_mask = geos_ds['LANDMASK'].values

    for fcst_init_month in range(1,13):
        for  ensemble_number in range (1,11):
            # Generate 3-Hourly Radiation Files per Lead Month
            for lead_month in range(NMONTHS):
                # Generate datetime of forecast
                fcst_datetime = datetime(fcst_init_year, fcst_init_month, 1) +\
                    relativedelta(months=lead_month)

                # Formatted forecast month
                fcst_year = fcst_datetime.year
                fcst_month = fcst_datetime.month

                # check if 3-Hourly Radiation File  exists
                filename_rad = generate_geosv3_3hr_radiation_filename(
                    fcst_init_year, fcst_init_month, ensemble_number, lead_month, dir_rad)
                if os.path.exists(filename_rad):
                    print(f'EXISTS: {filename_rad}')
                    continue

                # Read 3-Hourly Surface File
                filename_sfc = generate_geosv3_3hr_surface_filename(
                    fcst_init_year, fcst_init_month, ensemble_number, lead_month, dir_geos)
                if not os.path.exists(filename_sfc):
                    print(f'MISSING: {filename_sfc}')
                    continue
                vegtype, lai, asnow, ts, frland, lw_net, sw_net = read_sfc(filename_sfc)
                print(filename_sfc)
                if vegtype is None:
                    continue

                # Read 3-Hourly Albedo File; Check if Albedo file needs to be trimmed
                is_feb_and_not_leap_year = (fcst_month == 2) & (not calendar.isleap(fcst_year))

                filename_albedo = generate_geosv3_3hr_albedo_filename(fcst_month, dir_geos)
                albedo = read_albedo(filename_albedo, is_feb_and_not_leap_year)

                # Compute Radiation Components
                emissivity = compute_land_emissivity(vegtype, lai, asnow)
                lw_upward = compute_lw_upward(emissivity, ts)
                lw_geos = compute_lw_downward(lw_net, lw_upward)

                sw_geos = compute_sw_downward(albedo, sw_net, frland)
                lw_downward, sw_downward = \
                    fill_radiation_gaps_nearest_neighbor(geos_mask, ldt_mask, lw_geos,
                                                         sw_geos)

                # write 3-Hourly Radiation File
                write_nc4_file(filename_rad, longwave_down_name, lw_downward,
                               'surface_downward_longwave_land', shortwave_down_name, sw_downward,
                               'surface_downward_shortwave_land', filename_sfc)

#
# Main
#
if __name__ == "__main__":
    main()
