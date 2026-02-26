#!/usr/bin/python3
""" Spatially disaggregates the downloaded CFSv2 forecast data

Workflow
--------
1. Get grib2 filename and determine if file needs to be patched/replaced
2. Apply Regridder to file
3. [OPTIONAL] Apply lapse-rate correction
4. [OPTIONAL] Apply slope-aspect correction
5. [OPTIONAL] Apply landmask
6. [OPTIONAL] Clip unfeasable values
7. Select only time steps for the number of months wanted
8. Generate monthly mean of current lead_month data
9. Write variables to NetCDF
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
import xarray as xr
import pandas as pd
import cfgrib
from dateutil.relativedelta import relativedelta
from bcsd_filename_functions import generate_raw_6hourly_forecast_filename
from bcsd_cfsv2_functions import generate_downloaded_cfsv2_filename, read_grib2_file, create_cfsv2_elevation_difference_file
from bcsd_helper_functions import apply_regridder_weights_file, create_lat_lon_xr_dataset
from bcsd_helper_functions import VarLimits as lim

#
# Methods
#

def _usage():
    """Print command line usage."""
    txt =  f'[INFO] Usage: {sys.argv[0]}'
    txt += 'config_filename fcst_init_year fcst_init_month ensemble_number'
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
        'config_filename'   : str(sys.argv[1]),
        'fcst_init_year'    : int(sys.argv[2]), #2015
        'fcst_init_month'   : int(sys.argv[3]), #5
        'ensemble_number'   : int(sys.argv[4]), #1
        'lead_month'        : int(sys.argv[5]), #0
        'dir_main'          : str(config['SETUP']['DIR_MAIN']),  #'/discover/nobackup/projects/..'
        'fcst_name'         : str(config['BCSD']['fcst_data_type']),
        'dataset_type'      : str(config['SETUP']['DATATYPE']), #'hindcast'
        'nmonths'           : int(config['EXP']['lead_months']), #9
        'dir_supplementary' : str(config['SETUP']['supplementarydir'],
        'dir_patch_files'   : str(config['SETUP']['supplementarydir'] + '/bcsd_fcst/patch_files/'),
    }
    args['config'] = config

    return args

def get_original_or_patch_file_path(
    fcst_init_year : int, fcst_init_month : int, ensemble_number : int, variable_name : str,
    dir_patch_files : str,
    dir_downloads : str = '/discover/nobackup/projects/lis/MET_FORCING/CFSv2') -> str:

    # read patch_files_list.txt
    patch_list = f"{dir_patch_files}/patch_files_list.txt"
    df_ = pd.read_csv(patch_list, sep=',', engine='python',
                     header=0, names=['Time','Bad','Replace'])

    # check replace if necessary
    file_name = generate_downloaded_cfsv2_filename(
        fcst_init_year, fcst_init_month, ensemble_number, variable_name, dir_downloads, 'file_name')
    
    
    patch_list_file_path = args['dir_patch_files'] + '/patch_files_list.txt'
    df_patch_list = pd.read_csv(
        patch_list_file_path, sep=', ', engine='python', header=0, names=['Time','Bad','Replace'])
    df_patch_list = df_patch_list.sort_values(by=['Time'])
    df_sub = df_patch_list[df_patch_list.Bad == file_name]


    # If there is no 'bad' file, return the normal file_path
    if df_sub.empty:
        file_path = generate_downloaded_cfsv2_filename(
            fcst_init_year, fcst_init_month, ensemble_number, variable_name)
        print(f"[INFO] Using file: {file_name}")
    else:
        # If there is a single instance of a 'bad' file in the list
        if len(df_sub.index) == 1:
            patch_file_name = df_sub.Replace.values[0]
            
            file_path = args['dir_patch_files'] + patch_file_name
            print(f"[INFO] File not available: {file_name}")
            print(f"[INFO] Using alternate file: {patch_file_name}")
        else:
            # Error if there are multiple instances of a 'bad' file in the list
            raise RuntimeError('There are multiple instances of bad values in the patch files list.')
    
    return file_path

def apply_mask_to_data_array(ds_in : xr.Dataset):
    file_directory_landmask = '/discover/nobackup/projects/servir-s2s/razamora/supplementary_files/lis_darun/'
    file_name_landmask = 'lis_input.global_5km_hydroscs_mask_topo_merit1km_wclimppt.nc'
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

    for itime in range(ds_in.time.size):

        # Convert to numpy values
        force_tmp = ds_in['T2M'].values[itime,:,:]
        force_hum = ds_in['Q2M'].values[itime,:,:]
        force_lwd = ds_in['LWS'].values[itime,:,:]
        force_prs = ds_in['PS'].values[itime,:,:]

        # Temperature
        tcforce = force_tmp + (lapse*elevdiff)

        # Pressure
        tbar = (force_tmp + tcforce)/2
        pcforce = force_prs/(np.exp((LIS_CONST_G*elevdiff)/(rdry*tbar)))

        # Humidity
        force_hum = np.where(force_hum == 0, 1e-08, force_hum)
        ee = (force_hum*force_prs)/0.622               
        esat = 611.2*np.exp((17.67*(force_tmp - LIS_CONST_TKFRZ))/((force_tmp - LIS_CONST_TKFRZ) + 243.5))
        qsat = (0.622*esat)/(force_prs - (0.378*esat))
        rh = (force_hum/qsat)*100.0
        fesat = 611.2*np.exp((17.67*(tcforce - LIS_CONST_TKFRZ))/((tcforce - LIS_CONST_TKFRZ) + 243.5))
        fqsat = (0.622*fesat)/(pcforce-(0.378*fesat))
        hcforce = (rh*fqsat)/100.0

        # Longwave Radiation
        fe = (hcforce*pcforce)/0.622
        mee = ee / 100.0
        mfe = fe / 100.0
        #----------------------------------------------------------------------
        # correct for negative vapor pressure at very low temperatures at
        # high latitudes
        #----------------------------------------------------------------------
        mee = np.where(mee < 0, 1e-08, mee)
        mfe = np.where(mfe < 0, 1e-08, mfe)
        emiss  =1.08*(1 - np.exp(-mee**(force_tmp/bb)))
        femiss =1.08*(1 - np.exp(-mfe**(tcforce/bb)))
        ratio=(femiss*(tcforce**4))/(emiss * (force_tmp**4))
        lcforce = force_lwd * ratio
        
        # Output
        ds_out['T2M'].values[itime,:,:] = tcforce
        ds_out['Q2M'].values[itime,:,:] = hcforce
        ds_out['LWS'].values[itime,:,:] = lcforce
        ds_out['PS'].values[itime,:,:]  = pcforce

    return ds_out

def get_timezone(lon_value):
    if lon_value > 180:
        raise ValueError('Longitude value may not be > 180')
    if lon_value < -180:
        raise ValueError('Longitude value may not be < -180')
    
    zone = int((lon_value + 202.51)/15)
    zone = zone - 24 if zone > 24 else zone
    
    return zone

def correct_swddirect(aslope, aspect, saz, solzen, swddirect):
    deg2rad = math.pi/180.0
    delaz = abs(aspect - saz) * deg2rad
    sdircorr = np.vectorize(math.sin)(solzen)*np.vectorize(math.sin)(aslope)*np.vectorize(math.cos)(delaz)
    sdircorr = np.where(solzen <= 85*deg2rad, 
                        sdircorr/np.vectorize(math.cos)(solzen),
                        sdircorr/math.cos(85))
    swddirect = swddirect*(np.vectorize(math.cos)(aslope) + sdircorr)
    swddirect = np.where(swddirect < 0, 0, swddirect)
    return swddirect


def apply_slope_aspect_correction_to_data(ds_in : xr.Dataset):

	ds_out = ds_in.copy(deep=True)

	# Constants
	deg2rad = math.pi/180.0

	# Read in Sl
	file_path_merit = '/discover/nobackup/projects/servir-s2s/razamora/supplementary_files/lis_darun/lis_input.global_5km_hydroscs_mask_topo_merit1km_wclimppt.nc'
	ds_merit = xr.open_dataset(file_path_merit)

	lon_d  = ds_merit.lon.values
	lat_d  = ds_merit.lat.values
	lat_r  = lat_d*deg2rad
	slope  = ds_merit['SLOPE'].values/deg2rad
	aspect = ds_merit['ASPECT'].values/deg2rad

	for index_time in range(ds_in.time.size):
		# Setup variables
		swd = ds_in['SLRSF'].isel(time = index_time).values

		# swd = np.where(swd > 0, swd, np.nan)

		# We include the full calculation here, even though minute & second are 0
		timestamp = pd.Timestamp(ds_in.time.isel(time = index_time).values)
		gmt = timestamp.hour + timestamp.minute/60 + timestamp.second/3600
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
		gamma = 2*math.pi*(doy - 1)/365.
		decl = (0.006918 - 0.399912*math.cos(gamma) + 0.070257*math.sin(gamma) - 0.006758*math.cos(2*gamma) + 0.000907*math.sin(2*gamma) - 0.002697*math.cos(3*gamma)+0.00148*math.sin(3*gamma))
		et = (0.000075 + 0.001868*math.cos(gamma) - 0.032077*math.sin(gamma) - 0.014615*math.cos(2*gamma) - 0.04089*math.sin(2*gamma))*229.18

		ls = ((zone - 1)*15) - 180.
		lcorr = 4.*(ls - lon_d)*(-1)

		latime = lhour + lcorr/60. + et/60.
		latime = np.where(latime < 0, latime + 24, latime)
		latime = np.where(latime > 24, latime - 24, latime)

		omegad = (latime - 12)*(-15)
		omega = omegad*deg2rad

		#
		cosz = math.sin(decl)*np.vectorize(math.sin)(lat_r) + math.cos(decl)*np.vectorize(math.cos)(lat_r)*np.vectorize(math.cos)(omega)
		cosz = np.where(cosz < 0, 0, cosz)
		cosz = np.where(cosz > 1, 1, cosz)

		#
		sunang = np.where(cosz < 0.01764, 0.01764, cosz)
		cloud  = (1160.0*sunang - swd)/(963.0*sunang)
		cloud  = np.where(cloud < 0, 0, cloud)
		cloud  = np.where(cloud > 1, 1, cloud)

		#
		difrat = np.where(abs(sunang - 0.0223) > 0.0001, 0.0604 / (sunang - 0.0223) + 0.0683, 1)
		difrat = np.where(difrat < 0, 0, difrat)
		difrat = np.where(difrat > 1, 1, difrat)
		difrat = difrat + (1.0 - difrat)*cloud

		#
		vnrat = (580.0 - cloud*464.0)/((580.0 - cloud*499.0) + (580.0 - cloud*464.0))
		swddirect  = swd * ((1.0 - difrat)*vnrat + (1.0 - difrat)*(1.0 - vnrat))
		swddiffuse = swd * (difrat*vnrat + difrat*(1.0 - vnrat))

		#
		thour = 0 if (hr - 24) <= 0 else hr
		mody = yr - math.floor(yr*.25)*4
		if abs(mody) > 0:
		    fyear = (2.0*math.pi/365.0)*(doy - 1.0 + (thour - 12.0)/24.0)
		else:
		    fyear = (2.0*math.pi/366.0)*(doy - 1.0 + (thour-12.0)/24.0)

		#
		#calculate the equation of time in minutes
		eqtime = 229.18*(7.5e-5 + 1.868e-3*math.cos(fyear) - 3.2077e-2*math.sin(fyear) - 1.4615e-2*math.cos(2*fyear) - 4.0849e-2*math.sin(2*fyear))

		#calculate the true solar time 
		#CHECK if this the timeoffset = lhour - LIS_rc%hr is correct
		time_offset = eqtime + 4.0*lon_d + 60.0*abs(lhour-hr)

		tst = thour*60.0 + mn + time_offset

		#solar hour angle
		ha = (tst*0.25 - 180.0)*deg2rad

		cosphi = np.vectorize(math.sin)(lat_r)*math.sin(decl) + np.vectorize(math.cos)(lat_r)*math.cos(decl)*np.vectorize(math.cos)(ha)
		phi = np.vectorize(math.acos)(cosphi)

		#
		costheta = np.where(np.vectorize(math.cos)(lat_r)*np.vectorize(math.sin)(phi) == 0,
		                    np.where(np.vectorize(math.sin)(lat_r)*cosphi - math.sin(decl) > 0, 1, -1),
		                    np.vectorize(math.sin)(lat_r)*cosphi - math.sin(decl))/(np.vectorize(math.cos)(lat_r)*np.vectorize(math.sin)(phi))
		#avoid floating point errors
		costheta = np.where(abs(costheta) > 1,
		                    np.where(costheta > 0, 1, -1),
		                    costheta)

		#
		saz = np.where(lat_r >= 0,
		               np.where(ha < 0,
		                        180.0 - np.vectorize(math.acos)(costheta)/deg2rad,
		                        180.0 + np.vectorize(math.acos)(costheta)/deg2rad),
		               np.where(ha < 0,
		                        np.vectorize(math.acos)(costheta)/deg2rad,
		                        360 - np.vectorize(math.acos)(costheta)/deg2rad))


		#
		aslope = np.vectorize(min)(1.57,slope*deg2rad)
		aslope = np.vectorize(max)(0, aslope)
		solzen = np.vectorize(math.acos)(cosz)

		swddirect = np.where((swd > 0) & (aslope > 0) & (aslope < 90*deg2rad),
		               correct_swddirect(aslope, aspect, saz, solzen, swddirect),
		               swddirect)

		#
		ds_out['SLRSF'].values[index_time,:,:] = swddirect + swddiffuse

	return ds_out

def apply_value_clipping(ds_in : xr.Dataset):
    
    ds_out = ds_in.copy(deep=True)
    
    limits = lim()
    
    for variable_name in variable_dict:
        # clip limits
        is_precip = True if variable_name == 'PRECTOT' else False
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
fcst_name    = args['fcst_name']
dataset_type = args['dataset_type']
nmonths      = args['nmonths']

variable_dict = {
    'PRECTOT': {'grib_name':'prate',   'var_name':'prate',  'method':'conservative'},
    'SLRSF':   {'grib_name':'dswsfc',  'var_name':'dswsfc', 'method':'conservative'},
    'PS':      {'grib_name':'pressfc', 'var_name':'pressfc','method':'conservative'},
    'T2M':     {'grib_name':'tmp2m',   'var_name':'tmp2m',  'method':'bilinear'},
    'LWS':     {'grib_name':'dlwsfc',  'var_name':'dlwsfc', 'method':'conservative'},
    'Q2M':     {'grib_name':'q2m',     'var_name':'q2m',    'method':'conservative'},
    'U10M':    {'grib_name':'wnd10m',  'var_name':'u10',    'method':'bilinear'},
    'V10M':    {'grib_name':'wnd10m',  'var_name':'v10',    'method':'bilinear'}
}

# Create ds_out grid
lon_start = -179.975
lon_end = 179.975
lat_start = -59.975
lat_end = 89.975
cell_step = 0.05
ds_out = create_lat_lon_xr_dataset(lon_start, lon_end, lat_start, lat_end, cell_step)

netcdf_attribute_dict = {
    'zlib': True, 'complevel': 6, 'shuffle': True, 'missing_value': -9999., '_FillValue': -9999.}

new_netcdf_attribute_dict = {
    'PRECTOT': netcdf_attribute_dict,
    'SLRSF': netcdf_attribute_dict,
    'PS': netcdf_attribute_dict,
    'T2M': netcdf_attribute_dict,
    'LWS': netcdf_attribute_dict,
    'Q2M': netcdf_attribute_dict,
    'U10M': netcdf_attribute_dict,
    'V10M': netcdf_attribute_dict
}

# Generate datetime information
dt_fcst_begin = datetime(fcst_init_year, fcst_init_month, 1)
dt_fcst_end   = dt_fcst_begin + relativedelta(months = nmonths, days = -1)
dt_fcst_range = pd.date_range(start = dt_fcst_begin, end = dt_fcst_end)

dt_begin = datetime(fcst_init_year, fcst_init_month, 1) + relativedelta(months = lead_month)
dt_end   = dt_begin + relativedelta(months = 1, days = -1)
dt_range = pd.date_range(start=dt_begin, end = dt_end)


elevdiff = create_cfsv2_elevation_difference_file(dir_supplementary)

ds_6hourly = xr.Dataset()

# 1. Get grib2s filename, determine if file needs to be patched/replaced, and read in
for variable_name in variable_dict:

    filename = get_original_or_patch_file_path(fcst_init_year, fcst_init_month, ensemble_number, variable_dict[variable_name]['grib_name'],args['dir_patch_files'])
    varname = variable_dict[variable_name]['var_name'] if variable_name in ['U10M', 'V10M'] else None
    
    fcst_input = read_grib2_file(filename, varname).sel(time = slice(dt_begin.strftime('%Y-%m-%d'), dt_end.strftime('%Y-%m-%d')))
    ds_6hourly = ds_6hourly.merge(fcst_input)

sys.stdout.write(f'Read raw files. {datetime.now().time()} \n')
# Iterate by day, and work on 6-hourly data   
for day_time in dt_range:
    lead_day = int(np.where(day_time == dt_fcst_range)[0])
    sys.stdout.write(f'{lead_day} \n')
    dt_current_day = pd.date_range(day_time, freq='6h', periods = 4)
    
    ds_current = ds_6hourly.sel(time = dt_current_day)

    # 2. Apply Regridder to data
    ds_current = apply_regridder_weights_file(ds_current, ds_out, variable_dict, dir_supplementary)
    sys.stdout.write(f'Applied Regridder. {datetime.now().time()} \n')
    
    # 3. [OPTIONAL] Apply landmask
    ds_current = apply_mask_to_data_array(ds_current)
    sys.stdout.write(f'Applied land mask. {datetime.now().time()} \n')
    
    # 4. [OPTIONAL] Apply lapse-rate correction
    ds_current = apply_lapse_rate_correction_to_data(ds_current, elevdiff)
    sys.stdout.write(f'Applied Lapse Rate. {datetime.now().time()} \n')

    # 5. [OPTIONAL] Apply slope-aspect correction
    ds_current = apply_slope_aspect_correction_to_data(ds_current)
    sys.stdout.write(f'Applied Slope Aspect. {datetime.now().time()} \n')
    
    # 6. [OPTIONAL] Apply landmask
    ds_current = apply_mask_to_data_array(ds_current)
    sys.stdout.write(f'Applied land mask again. {datetime.now().time()} \n')
    
    # 7. [OPTIONAL] Clip unfeasable values
    ds_current = apply_value_clipping(ds_current)

    # 8. Write 6-hourly variables to daily NetCDF files
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
