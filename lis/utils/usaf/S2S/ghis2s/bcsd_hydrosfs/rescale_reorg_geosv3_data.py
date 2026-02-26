#!/usr/bin/python3
"""Spatially disaggregates the downloaded GEOSv3 forecast data

Workflow
--------
1. Get grib2 filename and determine if file needs to be patched/replaced
2. Rename GEOSv3 variable names to HydroSFS convention
3. Apply Regridder to file
4. [OPTIONAL] Apply lapse-rate correction
5. [OPTIONAL] Apply slope-aspect correction
6. [OPTIONAL] Apply landmask
7. [OPTIONAL] Clip unfeasable values
9. Write 6-hourly variables to NetCDF
"""

#
# Modules
#
import sys
import os
from datetime import datetime
import math
import yaml
import numpy as np
import pandas as pd
import xarray as xr
from dateutil.relativedelta import relativedelta

#
# Custom Modules
#
from bcsd_filename_functions import generate_raw_6hourly_forecast_filename
from bcsd_geosv3_functions import (
    generate_geosv3_3hr_surface_filename,
    generate_geosv3_3hr_radiation_filename,
    create_geosv3_elevation_difference_file
)
from bcsd_helper_functions import (
    apply_regridder_weights_file,
    create_lat_lon_xr_dataset,
    get_domain_info
)
from bcsd_helper_functions import VarLimits as lim
from dict_variables import (
    get_geosv3_name,
    get_hydrosfs_name,
    get_all_geosv3_names,
    get_all_hydrosfs_names,
    get_all_key_variable_names,
    get_interp_method
)
#
# Functions
#
def _usage():
    """Print command line usage."""
    txt =  f'[INFO] Usage: {sys.argv[0]}'
    txt += 'config_filename fcst_init_year fcst_init_month ensemble_number lead_month'
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""

    with open(sys.argv[1], 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 6:
        print('[ERR] Invalid number of command line arguments!')
        _usage()
        sys.exit(1)

    args = {
        'config_filename'    : str(sys.argv[1]),
        'fcst_init_year'     : int(sys.argv[2]), #2015
        'fcst_init_month'    : int(sys.argv[3]), #5
        'ensemble_number'    : int(sys.argv[4]), #1
        'lead_month'         : int(sys.argv[5]), #0
        'dir_main'           : str(config['SETUP']['DIR_MAIN']),  #'/discover/nobackup/projects/..'
        'dir_geos'           : str(config['BCSD']['fcst_download_dir']),
        'dir_supplementary'  : str(config['SETUP']['supplementarydir']),
        'file_name_landmask' : str(config['SETUP']['ldtinputfile']),
        'fcst_name'          : str(config['BCSD']['fcst_data_type']),
        'dataset_type'       : str(config['SETUP']['DATATYPE']), #'hindcast'
        'nmonths'            : int(config['EXP']['lead_months']), #9
    }
    args['config'] = config

    return args

def read_and_resample_geos_data(
    fcst_init_year, fcst_init_month, ensemble_number, lead_month, dir_geos):

    # GEOSv3 names
    lwdown_geosv3_name = get_geosv3_name('LWdown')
    swdown_geosv3_name = get_geosv3_name('SWdown')

    # Variables to drop (not needed)
    vars_to_drop = ['VEGTYPE', 'ASNOW', 'FRLAND', 'LAI', 'LWLAND', 'SWLAND', 'TS']
    
    # Read 3-Hourly Surface File
    filename_sfc = generate_geosv3_3hr_surface_filename(
        fcst_init_year, fcst_init_month, ensemble_number, lead_month, dir_geos)
    # vegtype, lai, asnow, ts, frland, lw_net, sw_net = read_sfc(filename_sfc)
    ds_sfc = xr.open_dataset(filename_sfc, drop_variables=vars_to_drop)

    # Write 3-Hourly Radiation File
    filename_rad = generate_geosv3_3hr_radiation_filename(
        fcst_init_year, fcst_init_month, ensemble_number, lead_month, dir_geos)
    ds_rad = xr.open_dataset(filename_rad)
    
    ds_sfc[lwdown_geosv3_name] = (('time','lat','lon'), ds_rad[lwdown_geosv3_name].values)
    ds_sfc[swdown_geosv3_name] = (('time','lat','lon'), ds_rad[swdown_geosv3_name].values)
    ds_sfc = ds_sfc.resample(time='6h').mean(dim='time')

    return ds_sfc

def apply_mask_to_data_array(ds_in : xr.Dataset, dir_supplementary : str, file_name_landmask : str):
    file_directory_landmask = f"{dir_supplementary}/lis_darun/"
    file_path_landmask = file_directory_landmask + file_name_landmask

    # Open land mask file
    ds_landmask = xr.open_dataset(file_path_landmask)[['LANDMASK']]
    
    # Rename land mask dimensions
    ds_landmask = ds_landmask.rename_dims({'north_south' : 'lat', 'east_west' : 'lon'})
    
    # Apply land mask
    ds_out = ds_in.where(ds_landmask.LANDMASK > 0)
    
    return ds_out

def apply_lapse_rate_correction_to_data(ds_in : xr.Dataset, elevdiff):

    ds_out = ds_in.copy(deep=True)
    
    # Constants
    rdry = 287.
    lapse = -0.0065
    bb = 2016
    LIS_CONST_G = 9.80616
    LIS_CONST_TKFRZ = 273.16

    # HydroSFS Names
    tmp_name = get_hydrosfs_name('T')
    hum_name = get_hydrosfs_name('Q')
    lwd_name = get_hydrosfs_name('LWdown')
    prs_name = get_hydrosfs_name('PS')

    for itime in range(ds_in.time.size):

        # Convert to numpy values
        force_tmp = ds_in[tmp_name].values[itime,:,:]
        force_hum = ds_in[hum_name].values[itime,:,:]
        force_lwd = ds_in[lwd_name].values[itime,:,:]
        force_prs = ds_in[prs_name].values[itime,:,:]

        # Temperature
        tcforce = force_tmp + (lapse*elevdiff)

        # Pressure
        tbar = (force_tmp + tcforce)/2
        pcforce = force_prs / (np.exp((LIS_CONST_G * elevdiff) / (rdry * tbar)))

        # Humidity
        force_hum = np.where(force_hum == 0, 1e-08, force_hum)
        ee = (force_hum * force_prs) / 0.622               
        esat = 611.2 * np.exp(
            (17.67 * (force_tmp - LIS_CONST_TKFRZ)) / ((force_tmp - LIS_CONST_TKFRZ) + 243.5)
        )
        qsat = (0.622 * esat) / (force_prs - (0.378 * esat))
        rh = (force_hum / qsat) * 100.0
        fesat = 611.2 * np.exp(
            (17.67 * (tcforce - LIS_CONST_TKFRZ))/((tcforce - LIS_CONST_TKFRZ) + 243.5)
        )
        fqsat = (0.622 * fesat) / (pcforce - (0.378 * fesat))
        hcforce = (rh * fqsat) / 100.0

        # Longwave Radiation
        fe = (hcforce * pcforce) / 0.622
        mee = ee / 100.0
        mfe = fe / 100.0
        
        # Correct for negative vapor pressure at very low temperatures at high latitudes
        mee = np.where(mee < 0, 1e-08, mee)
        mfe = np.where(mfe < 0, 1e-08, mfe)
        emiss  = 1.08 * (1 - np.exp(-mee**(force_tmp / bb)))
        femiss = 1.08 * (1 - np.exp(-mfe**(tcforce / bb)))
        ratio = (femiss * (tcforce**4))/(emiss * (force_tmp**4))
        lcforce = force_lwd * ratio
        
        # Output
        ds_out[tmp_name].values[itime,:,:] = tcforce
        ds_out[hum_name].values[itime,:,:] = hcforce
        ds_out[lwd_name].values[itime,:,:] = lcforce
        ds_out[prs_name].values[itime,:,:] = pcforce

    return ds_out

def get_timezone(lon_value):
    """Returns the timezone integer (1-24) """
    if lon_value > 180:
        raise ValueError('Longitude value may not be > 180')
    if lon_value < -180:
        raise ValueError('Longitude value may not be < -180')
    
    zone = int((lon_value + 202.51) / 15)
    zone = zone - 24 if zone > 24 else zone
    
    return zone

def correct_swddirect(aslope, aspect, saz, solzen, swddirect):
    deg2rad = math.pi / 180.0
    delaz = abs(aspect - saz) * deg2rad
    sdircorr = np.sin(solzen) * np.sin(aslope) * np.cos(delaz)
    sdircorr = np.where(solzen <= 85 * deg2rad, sdircorr / np.cos(solzen), sdircorr / math.cos(85))
    swddirect = swddirect * (np.cos(aslope) + sdircorr)
    swddirect = np.where(swddirect < 0, 0, swddirect)
    return swddirect


def apply_slope_aspect_correction_to_data(ds_in : xr.Dataset):
    
    ds_out = ds_in.copy(deep=True)
    
    # HydroSFS Names
    swd_name = get_hydrosfs_name('SWdown')
    
    # Constants
    deg2rad = math.pi / 180.0

    # Read in Sl
    file_directory_merit = f"{dir_supplementary}/lis_darun"
    file_path_merit = f"{file_directory_merit}/{file_name_landmask}"
    ds_merit = xr.open_dataset(file_path_merit)

    lon_d  = ds_merit.lon.values
    lat_d  = ds_merit.lat.values
    lat_r  = lat_d*deg2rad
    slope  = ds_merit['SLOPE'].values / deg2rad
    aspect = ds_merit['ASPECT'].values / deg2rad

    for index_time in range(ds_in.time.size):
        # Setup variables
        swd = ds_in[swd_name].isel(time = index_time).values

        # We include the full calculation here, even though minute & second are 0
        timestamp = pd.Timestamp(ds_in.time.isel(time = index_time).values)
        gmt = timestamp.hour + timestamp.minute / 60 + timestamp.second / 3600
        hr = timestamp.hour
        yr = timestamp.year
        mn = timestamp.minute
        doy = timestamp.day_of_year

        # Calculate time zone (1-24)
        zone = np.vectorize(get_timezone)(lon_d)

        # Calculate local hour (0-23) 0 = midnight, 23 = 11:00 p.m.
        change = zone - 13
        lhour = gmt + change
        lhour = np.where(lhour < 0, lhour + 24, lhour)
        lhour = np.where(lhour > 23, lhour - 24, lhour)

        # Generate cosz and decl
        gamma = 2 * math.pi * (doy - 1) / 365.
        decl = (
            0.006918
            - 0.399912 * math.cos(gamma)
            + 0.070257 * math.sin(gamma)
            - 0.006758 * math.cos(2 * gamma)
            + 0.000907 * math.sin(2 * gamma)
            - 0.002697 * math.cos(3 * gamma)
            + 0.00148 * math.sin(3 * gamma)
        )
        
        # Equation of Time
        et = 229.18 * (
            0.000075
            + 0.001868 * math.cos(gamma)
            - 0.032077 * math.sin(gamma)
            - 0.014615 * math.cos(2 * gamma)
            - 0.04089 * math.sin(2 * gamma)
        )

        #
        ls = ((zone - 1) * 15) - 180.
        lcorr = 4.*(ls - lon_d) * (-1)

        #
        latime = lhour + lcorr / 60. + et / 60.
        latime = np.where(latime < 0, latime + 24, latime)
        latime = np.where(latime > 24, latime - 24, latime)

        #
        omegad = (latime - 12) * (-15)
        omega = omegad * deg2rad

        #
        cosz = math.sin(decl) * np.sin(lat_r) + math.cos(decl) * np.cos(lat_r) * np.cos(omega)
        cosz = np.where(cosz < 0, 0, cosz)
        cosz = np.where(cosz > 1, 1, cosz)

        #
        sunang = np.where(cosz < 0.01764, 0.01764, cosz)
        cloud  = (1160.0 * sunang - swd) / (963.0 * sunang)
        cloud  = np.where(cloud < 0, 0, cloud)
        cloud  = np.where(cloud > 1, 1, cloud)

        #
        difrat = np.where(abs(sunang - 0.0223) > 0.0001, 0.0604 / (sunang - 0.0223) + 0.0683, 1)
        difrat = np.where(difrat < 0, 0, difrat)
        difrat = np.where(difrat > 1, 1, difrat)
        difrat = difrat + (1.0 - difrat) * cloud

        #
        vnrat = (580.0 - cloud * 464.0) / ((580.0 - cloud * 499.0) + (580.0 - cloud * 464.0))
        swddirect  = swd * ((1.0 - difrat) * vnrat + (1.0 - difrat)*(1.0 - vnrat))
        swddiffuse = swd * (difrat * vnrat + difrat * (1.0 - vnrat))

        #
        thour = 0 if (hr - 24) <= 0 else hr
        mody = yr - math.floor(yr * .25)*4
        if abs(mody) > 0:
            fyear = (2.0*math.pi / 365.0) * (doy - 1.0 + (thour - 12.0) / 24.0)
        else:
            fyear = (2.0*math.pi / 366.0) * (doy - 1.0 + (thour - 12.0) / 24.0)

        # Calculate the equation of time in minutes
        eqtime = 229.18 * (
            7.5e-5
            + 1.868e-3 * math.cos(fyear)
            - 3.2077e-2 * math.sin(fyear)
            - 1.4615e-2 * math.cos(2 * fyear)
            - 4.0849e-2 * math.sin(2 * fyear)
        )

        # Calculate the true solar time 
        time_offset = eqtime + 4.0 * lon_d + 60.0 * abs(lhour - hr)
        tst = thour * 60.0 + mn + time_offset

        # Solar hour angle
        ha = (tst * 0.25 - 180.0) * deg2rad

        cosphi = np.sin(lat_r) * math.sin(decl) + np.cos(lat_r) * math.cos(decl) * np.cos(ha)
        phi = np.acos(cosphi)

        #
        costheta = np.where(np.cos(lat_r) * np.sin(phi) == 0,
                            np.where(np.sin(lat_r)*cosphi - math.sin(decl) > 0, 1, -1),
                            np.sin(lat_r) * cosphi - math.sin(decl)) / (np.cos(lat_r) * np.sin(phi))
        
        # Avoid floating point errors
        costheta = np.where(abs(costheta) > 1,
                            np.where(costheta > 0, 1, -1),
                            costheta)

        #
        saz = np.where(lat_r >= 0,
                       np.where(ha < 0,
                                180.0 - np.acos(costheta) / deg2rad,
                                180.0 + np.acos(costheta) / deg2rad),
                       np.where(ha < 0,
                                np.acos(costheta) / deg2rad,
                                360 - np.acos(costheta) / deg2rad))


        #
        aslope = np.vectorize(min)(1.57, slope * deg2rad)
        aslope = np.vectorize(max)(0, aslope)
        solzen = np.acos(cosz)

        swddirect = np.where((swd > 0) & (aslope > 0) & (aslope < 90 * deg2rad),
                       correct_swddirect(aslope, aspect, saz, solzen, swddirect),
                       swddirect)

        #
        ds_out[swd_name].values[index_time, :, :] = swddirect + swddiffuse

    return ds_out

def apply_value_clipping(ds_in : xr.Dataset):
    
    ds_out = ds_in.copy(deep=True)
    
    limits = lim()
    
    for variable_name in get_all_hydrosfs_names():
        # clip limits
        is_precip = True if variable_name == get_hydrosfs_name('Pr') else False
        ds_out[variable_name].values[:] = limits.clip_array(
            np.array(ds_out[variable_name].values[:]),
            var_name = variable_name,
            precip = is_precip)
    return ds_out
        
    
#
# Main
#

sys.stdout.write(f'Start. {datetime.now().time()} \n')
args = _read_cmd_args()

# Input Args
fcst_init_year  = args['fcst_init_year']
fcst_init_month = args['fcst_init_month']
ensemble_number = args['ensemble_number']
lead_month      = args['lead_month']

# Config Args
dir_main     = args['dir_main']
dir_geos     = args['dir_geos']
dir_supplementary = args['dir_supplementary']
file_name_landmask = args['file_name_landmask']
fcst_name    = args['fcst_name']
dataset_type = args['dataset_type']
nmonths      = args['nmonths']


# Options to set (by default: True)
option_landmask = True
option_lapse_rate_correction = True
option_slope_aspect_correction = True
option_clip_values = True


variable_dict = {
    get_hydrosfs_name(key_variable):{'method': get_interp_method(key_variable)} 
    for key_variable in get_all_key_variable_names()}
    
LATS, LONS = get_domain_info(args['config_filename'], coord = True)

# Create ds_out grid
lon_start = LONS[0]
lon_end = LONS[-1]
lat_start = LATS[0]
lat_end = LATS[-1]
cell_step = np.unique(np.round(np.diff(LONS),2))

ds_out = create_lat_lon_xr_dataset(lon_start, lon_end, lat_start, lat_end, cell_step)

netcdf_attribute_dict = {
    'zlib': True, 'complevel': 6, 'shuffle': True, 'missing_value': -9999., '_FillValue': -9999.}

new_netcdf_attribute_dict = {var_name: netcdf_attribute_dict for var_name in get_all_hydrosfs_names()}

# Generate datetime information
dt_fcst_begin = datetime(fcst_init_year, fcst_init_month, 1)
dt_fcst_end   = dt_fcst_begin + relativedelta(months = nmonths, days = -1)
dt_fcst_range = pd.date_range(start = dt_fcst_begin, end = dt_fcst_end)
sys.stdout.write(f'DT FCST RANGE: {dt_fcst_range} \n')

dt_begin = datetime(fcst_init_year, fcst_init_month, 1) + relativedelta(months = lead_month)
dt_end   = dt_begin + relativedelta(months = 1, days = -1)
dt_range = pd.date_range(start=dt_begin, end = dt_end)

# Generate elevation difference dataset
elevdiff = create_geosv3_elevation_difference_file(dir_supplementary)

# 1. Get filenames, read in, and resample to 6-hourly
ds_6hourly = read_and_resample_geos_data(fcst_init_year, fcst_init_month, ensemble_number, lead_month, dir_geos)

sys.stdout.write(f'Read raw files. {datetime.now().time()} \n')
# Iterate by day, and work on 6-hourly data   
for day_time in dt_range:

    # Select data for current day
    lead_day = dt_fcst_range.get_loc(day_time)
    dt_current_day = pd.date_range(day_time, freq='6h', periods = 4)    
    ds_current = ds_6hourly.sel(time = dt_current_day)
    sys.stdout.write(f'Select day: {lead_day} \n')

    # 2. Rename GEOSv3 variables to use HydroSFS naming convention
    rename_dict = dict(zip(get_all_geosv3_names(), get_all_hydrosfs_names()))
    ds_current = ds_current.rename(rename_dict)
    
    # 3. Apply Regridder to data
    ds_current = apply_regridder_weights_file(ds_current, ds_out, variable_dict, dir_supplementary)
    sys.stdout.write(f'Applied Regridder. {datetime.now().time()} \n')
    
    # 4. [OPTIONAL] Apply landmask
    if option_landmask:
        ds_current = apply_mask_to_data_array(ds_current, dir_supplementary, file_name_landmask)
        sys.stdout.write(f'Applied land mask. {datetime.now().time()} \n')
    
    # 5. [OPTIONAL] Apply lapse-rate correction
    if option_lapse_rate_correction:
        ds_current = apply_lapse_rate_correction_to_data(ds_current, elevdiff)
        sys.stdout.write(f'Applied Lapse Rate. {datetime.now().time()} \n')

    # 6. [OPTIONAL] Apply slope-aspect correction
    if option_slope_aspect_correction:
        ds_current = apply_slope_aspect_correction_to_data(ds_current)
        sys.stdout.write(f'Applied Slope Aspect. {datetime.now().time()} \n')
    
    # 7. [OPTIONAL] Apply landmask
    if option_landmask:
        ds_current = apply_mask_to_data_array(ds_current, dir_supplementary, file_name_landmask)
        sys.stdout.write(f'Applied land mask. {datetime.now().time()} \n')
    
    # 8. [OPTIONAL] Clip unfeasible values
    if option_clip_values:
        ds_current = apply_value_clipping(ds_current)

    # 9. Write 6-hourly variables to daily NetCDF files
    # Create file directory
    file_directory_6hourly_out = generate_raw_6hourly_forecast_filename(
        dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number, lead_day, 
        'file_directory')
    if not os.path.exists(file_directory_6hourly_out):
        os.makedirs(file_directory_6hourly_out)

    # Create file path
    file_path_6hourly_out = generate_raw_6hourly_forecast_filename(
        dir_main, fcst_name, dataset_type, fcst_init_year, fcst_init_month, ensemble_number, lead_day)

    ds_current.to_netcdf(file_path_6hourly_out, format='NETCDF4', mode = 'w',
                               encoding = new_netcdf_attribute_dict)
    sys.stdout.write(f'Wrote File. {datetime.now().time()} \n')

# if __name__ == "__main__":
#     _driver()
