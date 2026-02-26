#!/usr/bin/python3
"""
Create GEOS-S2S-V3 Downward Radiation Components

Workflow
--------
1. Prepare files and directories (input and output)
2. Read surface variables
3. Compute emissivity from vegtype, leaf area index, and fractional snow cover
4. Compute downward longwave radiation from emissivity, surface temperature, and fractional land cover
5. Read albedo variables
6. Compute downward shortwave radiation from albedo, net shortwave radiation, and fractional land cover

Adapted from: LWSWpostcal.py
Author: Ryan Zamora; LWSWpostcal.py written by the GMAO Team
"""
#
# Modules
#
import sys
import os
import calendar
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from datetime import datetime
from dateutil.relativedelta import relativedelta

#
# Custom Modules
#
from bcsd_geosv3_functions import (generate_geosv3_3hr_surface_filename, generate_geosv3_3hr_albedo_filename, generate_geosv3_3hr_radiation_filename)

# Name for downward shortwave and longwave radiation (Names follow GEOSv3 convention)
longwave_down_name = 'LWS'
shortwave_down_name = 'SWGDWN'

# Coordinates
nlat, nlon = 361, 720
latitudes = np.linspace(-90, 90, nlat)
longitudes = np.linspace(-180, 179.5, nlon)

# Constants
SIGMA = 5.67e-8
nmonths = 4
EMSVEG = np.array([0.99560, 0.99000, 0.99560, 0.99320, 0.99280, 0.99180])
EMSBARESOIL = 0.94120
EMSSNO = 0.99999

#
# Functions
#
def _usage():
    """Print command line usage"""
    txt = f"[INFO] Usage: {(sys.argv[0])} dir_geos fcst_init_year fcst_init_month ensemble_number"
    print(txt)
    print("[INFO] where")
    print('[INFO] dir_geos: Directory where the GEOS files are located (string)')
    print('[INFO] fcst_init_year: Initial forecast year (integer)')
    print('[INFO] fcst_init_month: Initial forecast month (integer)')
    print('[INFO] ensemble_number: Ensemble number (integer)')

def compute_land_emissivity(vegtype, lai, asnow):
    """Compute emissivity from vegetation type, leaf area index, and fractional snow cover"""
    emis_veg = EMSVEG[vegtype.astype(int) - 1]
    emis_land = emis_veg + (EMSBARESOIL - emis_veg) * np.exp(-lai)
    emissivity = emis_land * (1 - asnow) + EMSSNO * asnow
    return emissivity

def compute_lw_upward(emissivity, surface_temp, frland):
    """Compute upward longwave radiation from emissivity, surface temperature, and fractional land cover"""
    lw_upward = emissivity * SIGMA * surface_temp**4 * frland
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
            VEGTYPE = np.where(nc.variables['VEGTYPE'][:, :, :] > 6, 0, nc.variables['VEGTYPE'][:, :, :])
            return (
                VEGTYPE,
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

def write_nc4_file(output_filename, lw_var_name, lw_data, lw_long_name, sw_var_name, sw_data, sw_long_name, filename_sfc):
    """Writes out the components for Downward Shortwave Radiation & Downward Longwave Radiation"""
    try:
        with Dataset(filename_sfc, 'r') as sample_nc:
            time = np.arange(lw_data.shape[0])*180.
            latitudes = np.linspace(-90.0, 90.0, 361)
            longitudes = np.linspace(-180.0, 179.5, 720)
            lat_attrs = {a: sample_nc.variables['lat'].getncattr(a) for a in sample_nc.variables['lat'].ncattrs()}
            lon_attrs = {a: sample_nc.variables['lon'].getncattr(a) for a in sample_nc.variables['lon'].ncattrs()}
            time_attrs = {a: sample_nc.variables['time'].getncattr(a) for a in sample_nc.variables['time'].ncattrs()}
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
            lw_var = nc.createVariable(lw_var_name, np.float32, ('time', 'lat', 'lon'), fill_value=1.0e20)
            sw_var = nc.createVariable(sw_var_name, np.float32, ('time', 'lat', 'lon'), fill_value=1.0e20)
            
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
        pass


def main():
    """Main driver"""

    # Check for correct number of input arguments  
    if len(sys.argv) != 5:
        print('[ERR] Invalid number of command line arguments!')
        _usage()
        sys.exit(1)
        
    # Read input arguments
    dir_geos = sys.argv[1]
    fcst_init_year = int(sys.argv[2])
    fcst_init_month = int(sys.argv[3])
    ensemble_number = int(sys.argv[4])
    
    # Generate 3-Hourly Radiation Files per Lead Month
    for lead_month in range(nmonths):
    
        # Generate datetime of forecast
        fcst_datetime = datetime(fcst_init_year, fcst_init_month, 1) + relativedelta(months=lead_month)
        
        # Formatted forecast month
        fcst_year = fcst_datetime.year
        fcst_month = fcst_datetime.month
        
        # Read 3-Hourly Surface File
        filename_sfc = generate_geosv3_3hr_surface_filename(
            fcst_init_year, fcst_init_month, ensemble_number, lead_month, dir_geos)
        vegtype, lai, asnow, ts, frland, lw_net, sw_net = read_sfc(filename_sfc)
        if vegtype is None:
            continue
    
        # Read 3-Hourly Albedo File; Check if Albedo file needs to be trimmed
        is_feb_and_not_leap_year = (fcst_month == 2) & (not calendar.isleap(fcst_year))
                                                         
        filename_albedo = generate_geosv3_3hr_albedo_filename(fcst_month, dir_geos)
        albedo = read_albedo(filename_albedo, is_feb_and_not_leap_year)
    
        # Compute Radiation Components
        emissivity = compute_land_emissivity(vegtype, lai, asnow)
        lw_upward = compute_lw_upward(emissivity, ts, frland)
        lw_downward = compute_lw_downward(lw_net, lw_upward)

        sw_downward = compute_sw_downward(albedo, sw_net, frland)
    
        # Write 3-Hourly Radiation File
        filename_rad = generate_geosv3_3hr_radiation_filename(
            fcst_init_year, fcst_init_month, ensemble_number, lead_month, dir_geos)
        write_nc4_file(filename_rad, longwave_down_name, lw_downward, 'surface_downward_longwave_land', shortwave_down_name, sw_downward, 'surface_downward_shortwave_land', filename_sfc)

#
# Main
#
if __name__ == "__main__":
    main()
