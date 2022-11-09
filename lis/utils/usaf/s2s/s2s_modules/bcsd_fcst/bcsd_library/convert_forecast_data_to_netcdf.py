#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT:  convert_forecast_data_to_netcdf.py
#
# PURPOSE: Convert NNME GRIB2 files to netCDF.  Based on
# convert_forecast_data_to_netcdf.scr by Ryan Zamora.
#
# REQUIREMENTS as of 04 Nov 2021:
# * Python 3.9
# * Climate Data Operators (CDO) software
# * wgrib2 software
#
# REVISION HISTORY:
# 04 Nov 2021: Eric Kemp/SSAI, first version.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import numpy as np
import cfgrib
import xarray as xr
import subprocess
import sys

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

def _wgrib2_to_netcdf(grib2_file):
    ds = cfgrib.open_dataset(grib2_file, indexpath ="")
    return ds
def _patch_corrupted_file(settings, args):
    """Patch corrupted GRIB2 file."""

    patched_timestring = settings['patched_timestring']
    index_subdaily_start = settings['index_subdaily_start']
    index_subdaily_end = settings['index_subdaily_end']
    index_patch_start = settings['index_patch_start']
    index_patch_end = settings['index_patch_end']
    patch_file = f"{args['patchdir']}/{args['file_pfx']}." + \
                 f"{patched_timestring}.{args['file_sfx']}"

    print(f"[INFO] Patching corrupted data: {args['subdaily_file']}")
    print(f"[INFO] Using data from: {patch_file}")

    # Subset corrupted file
    cmd = f"cdo -seltimestep,{index_subdaily_start}/{index_subdaily_end}" + \
        f" {args['subdaily_file']}" + \
        f" {args['outdir']}/junk1_{args['file_pfx']}_patch_A.grb2" + \
        " >/dev/null"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = "wgrib2 -v0" + \
        f" {args['outdir']}/junk1_{args['file_pfx']}_patch_A.grb2" + \
        f" -netcdf {args['outdir']}/junk1_{args['file_pfx']}_patch_A.nc" + \
        " >/dev/null"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)
    print("[INFO] Subset corrupted file.")

    # Subset patch file
    cmd = f"cdo -seltimestep,{index_patch_start}/{index_patch_end}" + \
        f" {args['patch_file']}" + \
        f" {args['outdir']}/junk1_{args['file_pfx']}_patch_B.grb2" + \
        " >/dev/null"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = "wgrib2 -v0" + \
        f" {args['outdir']}/junk1_{args['file_pfx']}_patch_B.grb2" + \
        f" -netcdf {args['outdir']}/junk1_{args['file_pfx']}_patch_B.nc" + \
        " >/dev/null"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)
    print("[INFO] Subset corrupted file.")

    # Merge files
    cmd = "cdo mergetime" + \
        f" {args['outdir']}/junk1_{args['file_pfx']}_patch_A.nc" + \
        f" {args['outdir']}/junk1_{args['file_pfx']}_patch_B.nc" + \
        f" {args['outdir']}/junk1_{args['file_pfx']}_{args['temp_name']}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)
    print("[INFO] Merged times")

def _replace_missing_file(patched_timestring, args):
    """Replace missing GRIB2 file."""
    patch_file = f"{args['patchdir']}/{args['file_pfx']}." + \
        f"{patched_timestring}.{args['file_sfx']}"
    cmd = f"wgrib2 -v0 {patch_file} -netcdf" + \
        f" {args['outdir']}/junk1_{args['file_pfx']}_{args['temp_name']}" + \
        " >/dev/null"
    ds = _wgrib2_to_netcdf (patch_file)
    print(f"[INFO] File not available: {args['subdaily_file']}")
    print(f"[INFO] Using alternate data: {patch_file}")
    return ds

def _handle_cases_with_corrupted_wind_files(settings, args):
    """Handle cases with corrupted wind GRIB2 files"""
    patched_timestring = settings['patched_timestring']
    if args['varname'] == "wnd10m":
        _replace_missing_file(patched_timestring, args)
    else:
        _patch_corrupted_file(settings, args)

def _patch_one_file(args):
    """Check and patch files where only one file needs to be replaced."""
    comparison_name = args['comparison_name']
    if comparison_name == "dlwsfc.01.2018122700.daily.grb2":
        settings = \
            _make_settings_dict(2018122600, 1, 1015, 1020, 1116)
        _patch_corrupted_file(settings, args)
        sys.exit(0)        
    if comparison_name == "q2m.2011022512.time.grb2":
        settings = \
            _make_settings_dict(2011022506, 1, 403, 405, 1239)
        _patch_corrupted_file(settings, args)
        sys.exit(0)
    if comparison_name == "q2m.2011022518.time.grb2":
        settings = \
            _make_settings_dict(2011022506, 1, 829, 832, 1239)
        _patch_corrupted_file(settings, args)
        sys.exit(0)
    if comparison_name == "dlwsfc.2011032706.time.grb2":
        settings = \
            _make_settings_dict(2011032700, 1, 121, 123, 1244)
        _patch_corrupted_file(settings, args)
        sys.exit(0)

def _patch_all_files(args):
    """Check and patch files where all variables need to be replaced"""
    fcst_timestring = args['fcst_timestring']
    if fcst_timestring == "2016032718":
        settings = \
            _make_settings_dict(2016032518, 1, 180, 189, 1125)
        _handle_cases_with_corrupted_wind_files(settings, args)
        sys.exit(0)
    if fcst_timestring == "2016042618":
        settings = \
            _make_settings_dict(2016042518, 1, 32, 37, 1125)
        _handle_cases_with_corrupted_wind_files(settings, args)
        sys.exit(0)
    if fcst_timestring == "2016063018":
        settings = \
            _make_settings_dict(2016070218, 1, 180, 173, 1209)
        _handle_cases_with_corrupted_wind_files(settings, args)
        sys.exit(0)
    if fcst_timestring == "2016073018":
        settings = \
            _make_settings_dict(2016072918, 1, 180, 185, 1225)
        _handle_cases_with_corrupted_wind_files(settings, args)
        sys.exit(0)
    if fcst_timestring == "2016082918":
        settings = \
            _make_settings_dict(2016090318, 1, 180, 161, 1201)
        _handle_cases_with_corrupted_wind_files(settings, args)
        sys.exit(0)

def _patch_missing_files(args):
    """Check and patch files that are completely missing."""
    fcst_timestring = args['fcst_timestring']
    ds = None
    if fcst_timestring == "2011082406":
        ds = _replace_missing_file(2011082306, args)
    if fcst_timestring == "2018062518":
        ds = _replace_missing_file(2018062418, args)
    if fcst_timestring == "2018063018":
        ds = _replace_missing_file(2018062818, args)
    if fcst_timestring == "2018073018":
        ds = _replace_missing_file(2018072618, args)
    if fcst_timestring == "2018082918":
        ds = _replace_missing_file(2018082818, args)
    if fcst_timestring == "2022012100":
        ds = _replace_missing_file(2022012000, args)
    return ds

def magnitude(a, b):
    func = lambda x, y: np.sqrt(x**2 + y**2)
    return xr.apply_ufunc(func, a, b)

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
    
    ds = None
    #_patch_one_file(args) # Will return if not applicable
    #_patch_all_files(args) # Will return if not applicable
    ds = _patch_missing_files(args) # Will return if not applicable

    if ds is None:
        # If we reach this point, we assume the file is fine.
        print("[INFO] File is normal.")
        ds = _wgrib2_to_netcdf(args['subdaily_file'])
        #cmd = f"wgrib2 -v0 {args['subdaily_file']} -netcdf"
        #cmd += f" {args['outdir']}/junk1_{args['file_pfx']}_{args['temp_name']}"
        #cmd += " >/dev/null"
        #print(cmd)
        #subprocess.run(cmd, shell=True, check=True)
    if args["varname"] == "wnd10m":
        U10 = magnitude(ds.u10, ds.v10)
        ds['WIND10M'] = U10
        ds['WIND10M'].attrs['units'] = 'm/s'
        ds['WIND10M'].attrs['short_name']='wnd10m'
        ds['WIND10M'].attrs['long_name']='Wind Speed'
        ds['WIND10M'].attrs['level']='10 m above ground'
        
    return ds

#if __name__ == "__main__":
#    ds = _driver()
#    return rs
