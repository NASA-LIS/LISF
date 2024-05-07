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
# SCRIPT:  convert_forecast_data_to_netcdf.py
#
# PURPOSE: Convert CFSv2 GRIB2 files to netCDF.  Based on
# convert_forecast_data_to_netcdf.scr by Ryan Zamora.
#
# REQUIREMENTS as of 04 Nov 2021:
# * Python 3.9
#
# REVISION HISTORY:
# 04 Nov 2021: Eric Kemp/SSAI, first version.
#
#------------------------------------------------------------------------------
"""


# Standard modules
import sys
from datetime import datetime
import numpy as np
import cfgrib
import xarray as xr
import pandas as pd

# Internal functions
def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} indir file_pfx, fcst_timestring file_sfx"
    txt += " outdir temp_name varname patchdir"
    print(txt)

def _make_settings_dict(patched_timestring,
                        index_subdaily_start,
                        index_subdaily_end,
                        index_patch_start,
                        index_patch_end):
    settings = {
        'patched_timestring' : patched_timestring,
        'index_subdaily_start' : index_subdaily_start,
        'index_subdaily_end' : index_subdaily_end,
        'index_patch_start' : index_patch_start,
        'index_patch_end' : index_patch_end,
    }
    return settings

def wgrib2_to_netcdf(grib2_file):
    new_name = {
        'prate':'PRECTOT','sp':'PS','t2m':'T2M',
        'dlwrf':'LWS','dswrf':'SLRSF','sh2':'Q2M',
        'u10':'U10M','v10':'V10M','q':'Q2M',}
    ds_ = cfgrib.open_dataset(grib2_file, indexpath ="")
    for varname, da_ in ds_.data_vars.items():
        ds_ = ds_.rename({varname : new_name.get(varname)})
    return ds_

def magnitude(_a, _b):
    ''' computes wind magnitude u^2 + v^2'''
    func = lambda x, y: np.sqrt(x**2 + y**2)
    return xr.apply_ufunc(func, _a, _b)

def _check_replace_missing (args):

    ds_ = None
    # read patch_files_list.txt
    patch_list = args['patchdir'] + '/patch_files_list.txt'
    df_ = pd.read_csv(patch_list, sep=',', engine='python',
                     header=0, names=['Time','Bad','Replace'])
    df_sort = df_.sort_values(by=['Time'])
    df_sort['Time'] = pd.to_datetime(df_sort['Time'],format='%Y%m%d%H')

    # check replace if necessary
    replace_file = None
    check_file = args['file_pfx'] + '.' + args['fcst_timestring'] + '.' + args['file_sfx']
    fcst_timestring = args['fcst_timestring']
    this_time = datetime.strptime(fcst_timestring, '%Y%m%d%H')
    df_sub = df_sort[(df_sort.Time == this_time)]
    if not df_sub.empty:
        for index, row in df_sub.iterrows():
            if check_file in row['Bad']:
                replace_file = row['Replace']

    # read replace file if need be
    if replace_file is not None:
        patch_file = args['patchdir'] + '/' + replace_file.strip()
        ds_ = wgrib2_to_netcdf (patch_file)
        print(f"[INFO] File not available: {args['subdaily_file']}")
        print(f"[INFO] Using alternate data: {patch_file}")

    return ds_

def read_wgrib (argv1, argv2, argv3, argv4, argv5, argv6, argv7, argv8):
    """Main driver."""
    args = {
        "indir" : argv1,
        "file_pfx" : argv2,
        "fcst_timestring" : argv3,
        "file_sfx" : argv4,
        "outdir" : argv5,
        "temp_name" : argv6,
        "varname" : argv7,
        "patchdir" : argv8,
    }
    args['comparison_name'] = \
        f"{args['file_pfx']}." + \
        f"{args['fcst_timestring']}." + \
        f"{args['file_sfx']}"
    args['subdaily_file'] = \
        f"{args['indir']}/" + \
        f"{args['comparison_name']}"

    ds_ = None
    ds_ = _check_replace_missing (args)

    if ds_ is None:
        # If we reach this point, we assume the file is fine.
        print("[INFO] " + args['subdaily_file'])
        print("[INFO] File is normal.")
        ds_ = wgrib2_to_netcdf(args['subdaily_file'])

    if args["varname"] == "wnd10m":
        u10 = magnitude(ds_.U10M, ds_.V10M)
        ds_['WIND10M'] = u10
        ds_['WIND10M'].attrs['units'] = 'm/s'
        ds_['WIND10M'].attrs['short_name']='wnd10m'
        ds_['WIND10M'].attrs['long_name']='Wind Speed'
        ds_['WIND10M'].attrs['level']='10 m above ground'

    return ds_
