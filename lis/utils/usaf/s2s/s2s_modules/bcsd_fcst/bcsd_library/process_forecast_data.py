#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT:  process_forecast_data.py
#
# PURPOSE: Convert NNME GRIB2 files to netCDF.  Based on
# process_forecast_data.scr by Ryan Zamora.
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
import os
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
import xesmf as xe
import yaml
# pylint: disable=import-error
from convert_forecast_data_to_netcdf import read_wgrib
from bcsd_stats_functions import get_domain_info
# pylint: enable=import-error

# Internal functions
def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} syr eyr fcst_init_monthday outdir"
    txt += " forcedir grid_description patchdir ic1 ic2 ic3"
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""

    with open(sys.argv[5], 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 9:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    args = {
        "syr" : sys.argv[1],
        "ens_num" : sys.argv[2],
        "fcst_init_monthday" : sys.argv[3],
        "outdir" : sys.argv[4],
        "forcedir" : config['BCSD']["fcst_download_dir"],
        "patchdir" : config['SETUP']['supplementarydir'] + '/bcsd_fcst/patch_files/',
        "ic1" : sys.argv[6],
        "ic2" : sys.argv[7],
        "ic3" : sys.argv[8],
        "configfile": sys.argv[5],
    }
    ic1 = args['ic1']
    ic2 = args['ic2']
    ic3 = args['ic3']
    args['all_ensmembers'] = ["00", "06", "12", "18", \
                              "00", "06", "12", "18", \
                              "00", "06", "12", "18"]
    args['all_monthdays'] = [ic1, ic1, ic1, ic1, \
                             ic2, ic2, ic2, ic2, \
                             ic3, ic3, ic3, ic3]
    args['config'] = config
    return args

def _set_input_file_info(input_fcst_year, input_fcst_month, input_fcst_var):
    """Set input file information"""
    cutoff_refor_yyyymm = 201103
    cutoff_oper_yyyymm = 209912
    current_yyyymm = (100 * input_fcst_year) + input_fcst_month

    # Up to apr1 2011 - Refor_HPS (dlwsfc, dswsfc, q2m, wnd10m), Refor_FL
    # (prate, pressfc, tmp2m)
    if current_yyyymm <= cutoff_refor_yyyymm:
        if input_fcst_var in ["dlwsfc", "dswsfc", "q2m", "wnd10m"]:
            subdir = "Refor_HPS"
        elif input_fcst_var in ["prate", "pressfc", "tmp2m"]:
            subdir = "Refor_FL"
        file_pfx = input_fcst_var
        file_sfx = "time.grb2"
    #From may01 2011 to jan01 2021 - Oper_TS (all fields)
    elif current_yyyymm <= cutoff_oper_yyyymm:
        subdir = "Oper_TS"
        file_pfx = f"{input_fcst_var}.01"
        file_sfx = "daily.grb2"
    # From feb01 2021 - OperRT_TS (all fields)
    elif current_yyyymm > cutoff_oper_yyyymm:
        subdir = "OperRT_TS"
        file_pfx = f"{input_fcst_var}.01"
        file_sfx = "daily.grb2"
    return subdir, file_pfx, file_sfx

def _migrate_to_monthly_files(cfsv2, outdirs, fcst_init, args):
    outdir_6hourly = outdirs["outdir_6hourly"]
    outdir_monthly = outdirs["outdir_monthly"]
    final_name_pfx = f"{fcst_init['monthday']}.cfsv2."
    reftime = \
         f"{fcst_init['year']}-{fcst_init['month']}-{fcst_init['day']}"
    reftime += f",{fcst_init['hour']}:00:00,1hour"

    if args["run_regrid"] is True:
        # resample to the S2S grid
        cfsv2["slice"] = cfsv2["T2M"].isel(step=0)

        # build regridder
        lats, lons = get_domain_info(args["configfile"], coord=True)
        ds_out = xr.Dataset(
            {
                "lat": (["lat"], lats),
                "lon": (["lon"], lons),
        }
        )
        args["regridder"] = xe.Regridder(cfsv2, ds_out, "bilinear", periodic=True)
        args["run_regrid"] = False

    # apply to the entire data set
    ds_out = args["regridder"](cfsv2)
    ds_out2 = ds_out.drop_vars('slice', errors="ignore")

    mmm = fcst_init['monthday'].split("0")[0].capitalize()
    dt1 = datetime.strptime('{} 1 {}'.format(mmm,fcst_init["year"]), '%b %d %Y')
    for month in range(1,10):
        file_6h = outdir_6hourly + '/' + final_name_pfx + '{:04d}{:02d}.nc'.format (dt1.year,dt1.month)
        file_mon = outdir_monthly + '/' + final_name_pfx + '{:04d}{:02d}.nc'.format (dt1.year,dt1.month)
        dt2 = dt1 + relativedelta(months=1)
        dt1s = np.datetime64(dt1.strftime('%Y-%m-%d'))
        dt2s = np.datetime64(dt2.strftime('%Y-%m-%d'))

        this_6h1 = ds_out2.sel(step = (ds_out2['valid_time']  >= dt1s) &
                                (ds_out2['valid_time']  < dt2s), drop=True)
        this_6h2 = this_6h1.rename_vars({"time": "time_step"})
        this_6h = this_6h2.rename_dims({"step": "time"})
        this_6h.to_netcdf(
            file_6h,format="NETCDF4",
            encoding = {
                "PRECTOT": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "PS": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "T2M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "LWS": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "SLRSF": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "Q2M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "U10M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "V10M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "WIND10M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.}})
        this_mon = this_6h.mean (dim='time')
        this_mon.to_netcdf(
            file_mon,format="NETCDF4",
            encoding = {
                "PRECTOT": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "PS": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "T2M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "LWS": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "SLRSF": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "Q2M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "U10M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "V10M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
                "WIND10M": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.}})
        dt1 = dt2

def _print_reftime(fcst_init, ens_num):
    """Print reftime to standard out"""
    reftime = \
        f"{fcst_init['year']}-{fcst_init['month']}-{fcst_init['day']}"
    reftime += f",{fcst_init['hour']}:00:00,1hour"
    txt = f"[INFO] ENS{ens_num}: " + \
        f"{fcst_init['year']}-{fcst_init['monthday']}:" + \
        f"{fcst_init['hour']}"
    print(txt)

def _driver():
    """Main driver."""
    args = _read_cmd_args()
    args['run_regrid'] = True
    args['regridder'] = 0.
    fcst_init = {}
    fcst_init["monthday"] = args['fcst_init_monthday']
    outdirs = {}

    year = int(args['syr'])
    ens_num = int(args['ens_num'])
    print(f"[INFO] {fcst_init['monthday']} {year}")

    if fcst_init['monthday'] == "jan01":
        fcst_init["year"] = year - 1
    else:
        fcst_init["year"] = year

    mmm = fcst_init['monthday'].split("0")[0].capitalize()
    dt1 = datetime.strptime('{} 1 {}'.format(mmm,fcst_init["year"]), '%b %d %Y')
    dt2 = dt1 + relativedelta(months=9)
    dt1 = np.datetime64(dt1.strftime('%Y-%m-%d'))
    dt2 = np.datetime64(dt2.strftime('%Y-%m-%d'))

    monthday = args['all_monthdays'][ens_num - 1]
    temp_name = f"cfsv2.{fcst_init['year']}{monthday}.nc"
    fcst_init['date'] = f"{fcst_init['year']}{monthday}"
    fcst_init['month'] = int(monthday[0:2])
    fcst_init['day'] = int(monthday[2:4])
    fcst_init['hour'] = args['all_ensmembers'][ens_num - 1]
    fcst_init['timestring'] = f"{fcst_init['date']}{fcst_init['hour']}"
    wanted_months = []
    for i in range(int(fcst_init['month']), 13):
        wanted_months.append(i)
    for i in range(1, int(fcst_init['month'])):
        wanted_months.append(i)
    wanted_months = wanted_months[1:10]

    _print_reftime(fcst_init, ens_num)

    outdirs['outdir_6hourly'] = \
        f"{args['outdir']}/6-Hourly/" + \
        f"{fcst_init['monthday']}/{year}/ens{ens_num}"
    if not os.path.exists(outdirs['outdir_6hourly']):
        os.makedirs(outdirs['outdir_6hourly'])
    outdirs['outdir_monthly'] = \
        f"{args['outdir']}/Monthly/" + \
        f"{fcst_init['monthday']}/{year}/ens{ens_num}"
    if not os.path.exists(outdirs['outdir_monthly']):
        os.makedirs(outdirs['outdir_monthly'])

    cfsv2 = []
    for varname in ["prate", "pressfc", "tmp2m", "dlwsfc", "dswsfc",
                    "q2m", "wnd10m"]:
        print(f"[INFO] {varname}")
        subdir, file_pfx, file_sfx = \
            _set_input_file_info(fcst_init['year'],
                                 fcst_init['month'],
                                 varname)
        indir = f"{args['forcedir']}/{subdir}/"
        if subdir == "Oper_TS" and not os.path.exists(indir):
            indir = f"{args['forcedir']}/"

        indir += f"{fcst_init['year']}/{fcst_init['date']}"

        # Convert GRIB file to netCDF and handle missing/corrupted data
        cfsv2.append(read_wgrib (indir, file_pfx, fcst_init['timestring'], file_sfx, outdirs['outdir_6hourly'], temp_name, varname, args['patchdir']))

    cfsv2 = xr.merge (cfsv2, compat='override')
    _migrate_to_monthly_files(cfsv2.sel (step = (cfsv2['valid_time']  >= dt1) &
                                         (cfsv2['valid_time']  < dt2)),
                              outdirs, fcst_init, args)

    print("[INFO] Done processing CFSv2 forecast files")

if __name__ == "__main__":
    _driver()
