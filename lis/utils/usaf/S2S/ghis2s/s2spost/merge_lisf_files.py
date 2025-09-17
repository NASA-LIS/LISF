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
#------------------------------------------------------------------------------
#
# SCRIPT: merge_lisf_files.py
#
# PURPOSE:  Merges daily netCDF output from LIS-NoahMP and LIS-HYMAP2 into
# single, CF-compliant netCDF4 file for distribution.  Rewritten to use
# NetCDF4 Python library instead of NCO software, to reduce runtime.
#
# REQUIREMENTS as of 02 Jun 2023:
# * Python 3.9 or higher.
# * UNIDATA NetCDF4 Python library
# * Numpy Python library
#
# REFERENCES:
# https://cfconventions.org for specifications of NetCDF Climate and Forecast
#   (CF) Metadata Conventions.
# https://unidata.ucar.edu/software/udunits for documentation on UDUNITS2
#   library, which CF is generally consistent with for unit specifications.
#
# REVISION HISTORY:
# 09 Sep 2022: Eric Kemp/SSAI, first version.
# 14 Nov 2022: K. Arsenault/SAIC, removed fields for FOC.
# 02 Jun 2023: K. Arsenault, updated the s2spost filenaming conventions
#------------------------------------------------------------------------------
"""

# Standard modules
import xarray as xr
import configparser
import copy
import datetime
import os
import sys

# Third-party libraries
# NOTE: pylint cannot see the Dataset class in netCDF4 since the latter is not
# written in Python.  We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module, too-many-branches, too-many-statements, too-many-locals
from netCDF4 import Dataset as nc4_dataset
import numpy as np
import yaml

# pylint: enable=no-name-in-module

# Private methods.
def create_final_filename(output_dir, fcst_date, curdt, model_forcing, domain):
    """Create final filename, following 557 convention."""
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing}"
    name += "_GR.C0P25DEG"
    if domain == 'AFRICOM':
        name += "_AR.AFRICA"
    if domain == 'GLOBAL':
        name += "_AR.GLOBAL"

    name += "_PA.ALL"
    name += f"_DD.{fcst_date.year:04d}{fcst_date.month:02d}{fcst_date.day:02d}"
    name += f"_DT.{curdt.hour:02d}00"
    name += f"_FD.{curdt.year:04d}{curdt.month:02d}{curdt.day:02d}"
    name += f"_DT.{curdt.hour:02d}00"
    name += "_DF.NC"
    if len(os.path.basename(name)) > 128:
        print("[ERR] Output file name is too long!")
        print(f"[ERR] {os.path.basename(name)} exceeds 128 characters!")
        sys.exit(1)
    return name

def merge_files_xarray(ldtfile, noahmp_file, hymap2_file, merge_file, fcst_date, logger, subtask):
    """Copy LDT, NoahMP and HYMAP2 fields into same file using xarray."""
    
    # Define dimension mapping
    dimension_dict = {
        "east_west": "lon",
        "north_south": "lat", 
        "SoilMoist_profiles": "soil_layer",
        "SoilTemp_profiles": "soil_layer",
        "RelSMC_profiles": "soil_layer",
    }
    
    # Define the attribute dictionaries
    _cell_methods = {
        "Evap_tavg": "time: mean area: point where land",
        "Qs_acc": "time: sum area: point where land",
        "Qsb_acc": "time: sum area: point where land",
        "AvgSurfT_tavg": "time: mean area: point where land",
        "SWE_tavg": "time: mean area: point where land",
        "SnowDepth_tavg": "time: mean area: point where land",
        "SoilMoist_tavg": "time: mean area: point where land",
        "SoilTemp_tavg": "time: mean area: point where land",
        "TWS_tavg": "time: mean area: point where land",
        "GWS_tavg": "time: mean area: point where land",
        "Snowcover_tavg": "time: mean area: point where land",
        "Wind_f_tavg": "time: mean",
        "Tair_f_tavg": "time: mean",
        "Tair_f_min": "time: minimum",
        "Tair_f_max": "time: maximum",
        "Qair_f_tavg": "time: mean",
        "Psurf_f_tavg": "time: mean",
        "SWdown_f_tavg": "time: mean",
        "LWdown_f_tavg": "time: mean",
        "RelSMC_tavg": "time: mean area: point where land",
        "TotalPrecip_acc": "time: sum",
        "Streamflow_tavg": "time: mean area: point where land",
        "RiverStor_tavg": "time: mean area: point where land",
        "FloodStor_tavg": "time: mean area: point where land",
        "FloodedFrac_tavg": "time: mean area: point where land",
        "FloodedArea_tavg": "time: mean area: point where land",
        "SWS_tavg": "time: mean area: point where land",
        "Elevation_inst": "area: point where land",
        "Greenness_inst": "area: point where land"
    }
    
    _new_standard_names = {
        "Evap_tavg": "water_evapotranspiration_flux",
        "SnowDepth_tavg": "surface_snow_thickness",
        "Elevation_inst": "height_above_mean_sea_level"
    }
    
    _remove_standard_names_list = []
    
    _new_units = {
        "SoilMoist_tavg": "1",
        "Qair_f_tavg": "1"
    }
    
    # Open datasets with xarray
    ds_noahmp = xr.open_dataset(noahmp_file)
    ds_hymap2 = xr.open_dataset(hymap2_file)
    ds_ldt = xr.open_dataset(ldtfile)
        
    def safe_rename_dims(dataset, rename_dict):
        existing_dims = {k: v for k, v in rename_dict.items() if k in dataset.sizes}
        if existing_dims:
            return dataset.rename(existing_dims)
        return dataset
    
    # Rename dimensions in all datasets (only if they exist)
    ds_noahmp = safe_rename_dims(ds_noahmp, dimension_dict)
    ds_hymap2 = safe_rename_dims(ds_hymap2, dimension_dict)
    ds_ldt = safe_rename_dims(ds_ldt, dimension_dict)
        
    # Start with NoahMP as the base dataset
    merged_ds = ds_noahmp.copy(deep=True)
    
    # Handle coordinate variables for CF convention
    if 'lat' in merged_ds.variables and len(merged_ds['lat'].dims) == 2:
        merged_ds = merged_ds.assign_coords({
            'lat': merged_ds['lat'].isel(lon=0),
            'lon': merged_ds['lon'].isel(lat=0)
        })
    
    # Find soil variables by checking what actually exists
    potential_soil_vars = ["SoilMoist_tavg", "SoilTemp_tavg", "RelSMC_tavg"]
    actual_soil_vars = [var for var in potential_soil_vars if var in merged_ds.data_vars]
    
    #print("Found soil variables:", actual_soil_vars)
    
    # Transpose soil variables to match CF convention
    for var in actual_soil_vars:
        #print(f"Processing soil variable {var} with dimensions: {merged_ds[var].dims}")
        
        dims = list(merged_ds[var].dims)
        if 'soil_layer' in dims:
            if 'ensemble' in dims and 'time' not in dims:
                merged_ds[var] = merged_ds[var].expand_dims('time', axis=1)
                target_dims = ['ensemble', 'time', 'soil_layer', 'lat', 'lon']
                current_dims = list(merged_ds[var].dims)
                final_dims = [dim for dim in target_dims if dim in current_dims]
                merged_ds[var] = merged_ds[var].transpose(*final_dims)
            elif 'time' not in dims:
                merged_ds[var] = merged_ds[var].expand_dims('time')
    
    # Handle other variables - add time dimension where needed
    for var in merged_ds.data_vars:
        if var not in actual_soil_vars and var not in ['lat', 'lon', 'time', 'ensemble']:
            current_dims = list(merged_ds[var].dims)
            #print(f"Processing variable {var} with dimensions: {current_dims}")
            
            if var in ["Landcover_inst", "Soiltype_inst", "Elevation_inst"]:
                # Keep 2D for these specific variables - remove ensemble if present
                if 'ensemble' in current_dims and len(current_dims) > 2:
                    merged_ds[var] = merged_ds[var].isel(ensemble=0, drop=True)
            elif len(current_dims) == 3 and 'time' not in current_dims:
                # Add time dimension for 3D variables
                merged_ds[var] = merged_ds[var].expand_dims('time', axis=1)
    
    # Apply landmask to TotalPrecip_acc and ensure correct dimension order
    if 'TotalPrecip_acc' in merged_ds and 'LANDMASK' in ds_ldt:
        landmask = ds_ldt['LANDMASK']
        precip = merged_ds['TotalPrecip_acc']
        masked_precip = xr.where(landmask == 1, precip, -9999.)
        
        if 'time' not in masked_precip.dims:
            masked_precip = masked_precip.expand_dims('time', axis=1)
        
        # Ensure correct dimension order (ensemble, time, lat, lon)
        target_dims = ['ensemble', 'time', 'lat', 'lon']
        current_dims = list(masked_precip.dims)
        final_dims = [dim for dim in target_dims if dim in current_dims]
        merged_ds['TotalPrecip_acc'] = masked_precip.transpose(*final_dims)
    
    # Add HYMAP2 variables
    src2_excludes = ["lat", "lon", "time", "ensemble", "RunoffStor_tavg", "BaseflowStor_tavg"]
    hymap2_vars = [var for var in ds_hymap2.data_vars if var not in src2_excludes]
    
    for var_name in hymap2_vars:
        var_data = ds_hymap2[var_name]
        current_dims = list(var_data.dims)
        
        # Add time dimension if needed
        if len(current_dims) == 3 and 'time' not in current_dims:
            var_data = var_data.expand_dims('time', axis=1)
        
        merged_ds[var_name] = var_data
    
    # Add LANDMASK from LDT
    if 'LANDMASK' in ds_ldt:
        merged_ds['LANDMASK'] = ds_ldt['LANDMASK']
    
    # Create soil layer coordinate if it doesn't exist
    if 'soil_layer' not in merged_ds.coords and 'soil_layer' in merged_ds.sizes:
        merged_ds = merged_ds.assign_coords({
            'soil_layer': xr.DataArray(np.array([1, 2, 3, 4], dtype=np.int32), dims=['soil_layer'])
        })
    
    # Add soil layer thickness
    if 'soil_layer' in merged_ds.sizes:
        merged_ds['soil_layer_thickness'] = xr.DataArray(
            np.array([0.1, 0.3, 0.6, 1.0], dtype=np.float32), 
            dims=['soil_layer'],
            attrs={
                "long_name": "soil layer thicknesses",
                "units": "m"
            }
        )
    
    # Add forecast reference time
    merged_ds['atime'] = xr.DataArray(
        0., 
        attrs={
            "standard_name": "forecast_reference_time",
            "units": "hours since " + fcst_date.strftime("%Y-%m-%d") + " 00:00"
        }
    )
        
    # Create time bounds
    if 'time' in merged_ds.sizes:
        time_size = merged_ds.sizes['time']
        time_bnds_data = np.full((time_size, 2), [[-1440., 0.]], dtype=np.float32)
        
        time_bnds = xr.DataArray(
            time_bnds_data,
            dims=['time', 'nv'],
            coords={'time': merged_ds.coords['time']},
            name='time_bnds'
        )
        merged_ds['time_bnds'] = time_bnds
    
    # Set global attributes
    attrs = copy.deepcopy(ds_noahmp.attrs)
    attrs["Conventions"] = "CF-1.8"
    if "conventions" in attrs:
        del attrs["conventions"]
    if "missing_value" in attrs:
        del attrs["missing_value"]
    if "SOIL_LAYER_THICKNESSES" in attrs:
        del attrs["SOIL_LAYER_THICKNESSES"]
    attrs["source"] = "Noah-MP.4.0.1+template open water+HYMAP2"
    merged_ds.attrs = attrs
        
    # Set coordinate attributes - let xarray handle all encoding automatically
    for coord_name in merged_ds.coords:
        coord_attrs = copy.deepcopy(merged_ds[coord_name].attrs)
        
        # Remove any encoding-related attributes - let xarray handle everything
        encoding_related = ['scale_factor', 'add_offset', '_FillValue', 'missing_value', 'units', 'calendar']
        for attr in encoding_related:
            if attr in coord_attrs:
                del coord_attrs[attr]
        
        if coord_name == "time":
            coord_attrs.update({
                "axis": "T",
                "bounds": "time_bnds",
                "begin_date": (fcst_date + datetime.timedelta(days=1)).strftime("%Y%m%d"),
                "begin_time": "000000",
                "long_name": "time",
            })
                
        elif coord_name == "lat":
            coord_attrs.update({
                "axis": "Y",
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_north",
                "vmin": 0.0,
                "vmax": 0.0
            })
            
        elif coord_name == "lon":
            coord_attrs.update({
                "axis": "X", 
                "standard_name": "longitude",
                "long_name": "longitude",
                "units": "degrees_east",
                "vmin": 0.0,
                "vmax": 0.0
            })
            
        elif coord_name == "ensemble":
            coord_attrs.update({
                "axis": "E",
                "units": "1",
                "long_name": "Ensemble numbers"
            })
                
        elif coord_name == "soil_layer":
            coord_attrs.update({
                "long_name": "soil layer level",
                "axis": "Z",
                "positive": "down"
            })
        
        merged_ds[coord_name].attrs = coord_attrs
    
    # Set data variable attributes 
    for var_name in merged_ds.data_vars:
        var_attrs = copy.deepcopy(merged_ds[var_name].attrs)
        
        # Remove any existing encoding-related attributes from data variables
        encoding_related = ['scale_factor', 'add_offset', '_FillValue']
        for attr in encoding_related:
            if attr in var_attrs:
                del var_attrs[attr]
        
        # Set missing_value for data variables (but let encoding handle _FillValue)
        if var_name not in ['atime', ]:
            var_attrs['missing_value'] = -9999.0
        
        # Specific attribute modifications for each variable type
        if var_name == "SoilMoist_tavg":
            var_attrs["long_name"] = "volumetric soil moisture content"
            var_attrs["units"] = "1"
            
        elif var_name == "Soiltype_inst":
            var_attrs.update({
                "flag_meanings": "sand loamy_sand sandy_loam silt_loam silt loam sandy_clay_loam silty_clay_loam clay_loam sandy_clay silty_clay clay organic_material water bedrock other+land-ice",
                "valid_range": [1., 16.],
                "flag_values": [np.float32(i) for i in range(1, 17)]
            })
            if "units" in var_attrs:
                del var_attrs["units"]
                
        elif var_name == "Landcover_inst":
            var_attrs.update({
                "flag_meanings": "evergreen_needleleaf_forest evergreen_broadleaf_forest deciduous_needleleaf_forest deciduous_broadleaf_forest mixed_forests closed_shrublands open_shrublands woody_savannas savannas grasslands permanent_wetlands croplands urban_and_built-up cropland+natural_vegetation_mosaic snow_and_ice barren_or_sparsely_vegetated water wooded_tundra mixed_tundra barren_tundra water",
                "valid_range": [1., 21.],
                "flag_values": [np.float32(i) for i in range(1, 22)]
            })
            if "units" in var_attrs:
                del var_attrs["units"]
                
        elif var_name == "LANDMASK":
            var_attrs.update({
                "flag_values": [np.float32(i) for i in range(0, 2)],
                "flag_meanings": "water land",
                "long_name": "land mask from LDT",
                "units": ""
            })
            if "standard_name" in var_attrs:
                del var_attrs["standard_name"]
                
        elif var_name == "TotalPrecip_acc":
            var_attrs.update({
                "units": "kg m-2",
                "standard_name": "precipitation_amount", 
                "long_name": "total precipitation amount",
                "vmin": -1.e+15,
                "vmax": 1.e+15
            })
        
        elif var_name == "atime":
            var_attrs.update({
                "standard_name": "forecast_reference_time",
                "units": "hours since " + fcst_date.strftime("%Y-%m-%d") + " 00:00"
            })
        
        # Apply common attribute modifications
        if var_name in _cell_methods:
            var_attrs["cell_methods"] = _cell_methods[var_name]
        if var_name in _new_standard_names:
            var_attrs["standard_name"] = _new_standard_names[var_name]
        if var_name in _remove_standard_names_list and "standard_name" in var_attrs:
            del var_attrs["standard_name"]
        if var_name in _new_units:
            var_attrs["units"] = _new_units[var_name]
            
        merged_ds[var_name].attrs = var_attrs
    
    # Remove time_increment attribute if it exists
    if 'time' in merged_ds and 'time_increment' in merged_ds['time'].attrs:
        del merged_ds['time'].attrs['time_increment']
    
    # Define very minimal encoding - only valid netCDF4 backend parameters
    true_coords = ['time', 'ensemble', 'soil_layer']
    special_data_vars = ['atime', 'time_bnds', 'soil_layer_thickness']
    encoding = {}
    
    # Handle true coordinates
    for var in true_coords:
        if var in merged_ds:
            if var == 'time':
                encoding[var] = {
                    '_FillValue': None,
                    'dtype': 'float32',
                }
            else:
                encoding[var] = {
                    '_FillValue': None
                }

    # Handle special data variables (no compression but need _FillValue=None)
    for var in special_data_vars:
        if var in merged_ds:
            encoding[var] = {'_FillValue': None}
            #if var == 'soil_layer':
            #    encoding[var]['dtype'] = 'int32'
            #elif var in ['soil_layer_thickness', 'time_bnds']:
            #    encoding[var]['dtype'] = 'float32'
            # Remove any existing _FillValue and missing_value from attributes
            if '_FillValue' in merged_ds[var].attrs:
                del merged_ds[var].attrs['_FillValue']
                if 'missing_value' in merged_ds[var].attrs:
                    del merged_ds[var].attrs['missing_value']

    # Handle normal data variables
    for var in merged_ds.data_vars:
        if var not in special_data_vars:
            encoding[var] = {
                'dtype': 'float32',
                'zlib': True,
                'complevel': 6,
                'shuffle': True,
                '_FillValue': -9999.0
            }
            # Add offset and scale factor attributes
            merged_ds[var].attrs['add_offset'] = 0.0
            merged_ds[var].attrs['scale_factor'] = 1.0

    try:
        merged_ds.to_netcdf(merge_file, format='NETCDF4', encoding=encoding)
    except Exception as e:
        logger.error(f"Error saving file: {e}", subtask=subtask)
        raise
    finally:
        ds_noahmp.close()
        ds_hymap2.close() 
        ds_ldt.close()

    return

