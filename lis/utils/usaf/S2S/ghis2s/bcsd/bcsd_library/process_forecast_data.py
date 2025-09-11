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
# SCRIPT:  process_forecast_data.py
#
# PURPOSE: Convert CFSv2 GRIB2 files to netCDF.  Based on
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
import gc

#from concurrent.futures import ProcessPoolExecutor
# pylint: disable=import-error
from ghis2s.bcsd.bcsd_library.convert_forecast_data_to_netcdf import read_wgrib
from ghis2s.shared.utils import get_domain_info
from ghis2s.bcsd.bcsd_library.bcsd_function import apply_regridding_with_mask
from ghis2s.bcsd.bcsd_library.bcsd_function import VarLimits as lim
from ghis2s.shared.logging_utils import TaskLogger
# pylint: enable=import-error
limits = lim()
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
        "config": config,
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

def write_monthly_files(this_6h1, file_6h, file_mon, logger, subtask):
    this_6h2 = this_6h1.rename_vars({"time": "time_step"})
    this_6h = this_6h2.rename_dims({"step": "time"})

    encoding = {
        var: {"zlib": True, "complevel": 6, "shuffle": True, "missing_value": -9999.}
        for var in ["PRECTOT", "PS", "T2M", "LWS", "SLRSF", "Q2M", "U10M", "V10M", "WIND10M"]
    }

    this_6h.to_netcdf(file_6h, format="NETCDF4", encoding=encoding)
    this_mon = this_6h.mean(dim="time")
    this_mon.to_netcdf(file_mon, format="NETCDF4", encoding=encoding)
    this_6h2.close()
    this_6h.close()
    this_mon.close()
    del this_6h2, this_6h, this_mon
    return

def _migrate_to_monthly_files(cfsv2_in, outdirs, fcst_init, args, rank, logger, subtask):
    regrid_method = {
        '25km': {'PRECTOT':'conservative', 'SLRSF':'bilinear', 'LWS':'bilinear','PS':'bilinear',
                 'Q2M':'bilinear', 'T2M':'bilinear', 'U10M':'bilinear', 'V10M':'bilinear', 'WIND10M':'bilinear'},
        '10km': {'PRECTOT':'conservative', 'SLRSF':'bilinear', 'LWS':'bilinear','PS':'bilinear',
                 'Q2M':'bilinear', 'T2M':'bilinear', 'U10M':'bilinear', 'V10M':'bilinear', 'WIND10M':'bilinear'},
        '5km': {'PRECTOT':'conservative', 'SLRSF':'conservative', 'LWS':'conservative','PS':'conservative',
                 'Q2M':'conservative', 'T2M':'bilinear', 'U10M':'bilinear', 'V10M':'bilinear', 'WIND10M':'bilinear'},}
    
    outdir_6hourly = outdirs["outdir_6hourly"]
    outdir_monthly = outdirs["outdir_monthly"]
    final_name_pfx = f"{fcst_init['monthday']}.cfsv2."
    reftime = \
         f"{fcst_init['year']}-{fcst_init['month']}-{fcst_init['day']}"
    reftime += f",{fcst_init['hour']}:00:00,1hour"

    # resample to the S2S grid
    # build regridder
    weightdir = args["config"]['SETUP']['supplementarydir'] + '/bcsd_fcst/'
    target_land_mask = xr.open_dataset(args["config"]['SETUP']['supplementarydir'] + '/lis_darun/' + \
                             args["config"]['SETUP']['ldtinputfile'])
    target_land_mask = target_land_mask.rename({'north_south': 'lat', 'east_west': 'lon'})
    
    lats, lons = get_domain_info(args["configfile"], coord=True)
    resol = round((lats[1] - lats[0])*100)
    resol = f'{resol}km'
    # read CFSv2 land mask
    cfsv2_land_mask = xr.open_dataset(weightdir + f'CFSv2_{resol}_landmask.nc4')
    
    '''
    resol='25km': NY=720, NX=1440; 
    resol='10km': NY=1500, NX=3600; 
    resol='5km': NY=3000, NX=7200; 
    cfsv2 dimensions: step: ~1151, latitude: 190, longitude: 384
    '''
    
    method = regrid_method.get(resol)
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], lats),
            "lon": (["lon"], lons),
        }
    )
    cfsv2_masked = cfsv2_in.copy()
    cfsv2_masked['mask'] = cfsv2_land_mask['LANDMASK']

    weight_file = weightdir + f'CFSv2_{resol}_bilinear_land.nc'
    bil_regridder = xe.Regridder(cfsv2_masked, ds_out, "bilinear", periodic=True, 
                                 reuse_weights=True,
                                 extrap_method='nearest_s2d', 
                                 filename=weight_file)

    weight_file = weightdir + f'CFSv2_{resol}_conservative.nc'
    con_regridder = xe.Regridder(cfsv2_in, ds_out, "conservative", periodic=True, 
                                 reuse_weights=True, 
                                 filename=weight_file)
    
    bilinear_vars = [var for var in cfsv2_in.data_vars if var in method and method[var] == 'bilinear']
    conservative_vars = [var for var in cfsv2_in.data_vars if var in method and method[var] == 'conservative']

    if bilinear_vars:
        cfsv2_bilinear = cfsv2_masked[bilinear_vars]
        result_bilinear = apply_regridding_with_mask(cfsv2_bilinear, bil_regridder,
                                                     cfsv2_land_mask, target_land_mask, 'bilinear')
        for var in bilinear_vars:
            ds_out[var] = result_bilinear[var]

    if conservative_vars:
        cfsv2_conservative = cfsv2_in[conservative_vars]
        result_conservative = apply_regridding_with_mask(cfsv2_conservative, con_regridder,
                                                         cfsv2_land_mask, None, 'conservative')
        for var in conservative_vars:
            ds_out[var] = result_conservative[var]
            
    mmm = fcst_init['monthday'].split("0")[0].capitalize()
    dt1 = datetime.strptime('{} 1 {}'.format(mmm,fcst_init["year"]), '%b %d %Y')
    # clip limits
    ds_out["PRECTOT"].values[:] = limits.clip_array(np.array(ds_out["PRECTOT"].values[:]),
                                                     var_name = "PRECTOT", precip=True)
    ds_out["PS"].values[:] = limits.clip_array(np.array(ds_out["PS"].values[:]),
                                                var_name = "PS")
    ds_out["T2M"].values[:] = limits.clip_array(np.array(ds_out["T2M"].values[:]),
                                                 var_name = "T2M")
    ds_out["LWS"].values[:] = limits.clip_array(np.array(ds_out["LWS"].values[:]),
                                                 var_name = "LWS")
    ds_out["SLRSF"].values[:] = limits.clip_array(np.array(ds_out["SLRSF"].values[:]),
                                                   var_name = "SLRSF")
    ds_out["Q2M"].values[:] = limits.clip_array(np.array(ds_out["Q2M"].values[:]),
                                                 var_name = "Q2M")
    ds_out["WIND10M"].values[:] = limits.clip_array(np.array(ds_out["WIND10M"].values[:]),
                                                     var_name = "WIND")

    dt1 = dt1 + relativedelta(months=rank)
    file_6h = outdir_6hourly + '/' + \
        final_name_pfx + '{:04d}{:02d}.nc'.format (dt1.year,dt1.month)
    file_mon = outdir_monthly + '/' + \
        final_name_pfx + '{:04d}{:02d}.nc'.format (dt1.year,dt1.month)

    write_monthly_files(ds_out, file_6h, file_mon)
    ds_out.close()
    cfsv2_masked.close()
    del ds_out, cfsv2_masked
    logger.info(f"Done processing CFSv2 forecast files rank: {rank}", subtask = subtask)
    return

def _print_reftime(fcst_init, ens_num):
    """Print reftime to standard out"""
    reftime = \
        f"{fcst_init['year']}-{fcst_init['month']}-{fcst_init['day']}"
    reftime += f",{fcst_init['hour']}:00:00,1hour"
    txt = f"[INFO] CFSv2 ENS-MEM #{ens_num}: " + \
        f"{fcst_init['year']}-{fcst_init['monthday']}:" + \
        f"{fcst_init['hour']} cycle"
    return txt

def _driver(rank):
    """Main driver."""
    args = _read_cmd_args()
    fcst_init = {}
    fcst_init["monthday"] = args['fcst_init_monthday']
    outdirs = {}
    year = int(args['syr'])
    ens_num = int(args['ens_num'])

    fcst_init["year"] = year
    if fcst_init['monthday'] == "jan01":
        fcst_init["year_cfsv2"] = year - 1
    else:
        fcst_init["year_cfsv2"] = year

    mmm = fcst_init['monthday'].split("0")[0].capitalize()
    dt0 = datetime.strptime('{} 1 {}'.format(mmm,fcst_init["year"]), '%b %d %Y')
    dt1 = dt0 + relativedelta(months=rank)
    dt2 = dt1 + relativedelta(months=1)
    task_name = os.environ.get('SCRIPT_NAME')
    subtask = f'{dt1.year:04d}{dt1.month:02d}'
    logger = TaskLogger(task_name,
                        os.getcwd(),
                        f'bcsd/process_forecast_data.py processing CFSv2 forcings for {subtask}')

    dt1 = np.datetime64(dt1.strftime('%Y-%m-%d'))
    dt2 = np.datetime64(dt2.strftime('%Y-%m-%d'))

    monthday = args['all_monthdays'][ens_num - 1]
    temp_name = f"cfsv2.{fcst_init['year_cfsv2']}{monthday}.nc"
    fcst_init['date'] = f"{fcst_init['year_cfsv2']}{monthday}"
    fcst_init['month'] = int(monthday[0:2])
    fcst_init['day'] = int(monthday[2:4])
    fcst_init['hour'] = args['all_ensmembers'][ens_num - 1]
    fcst_init['timestring'] = f"{fcst_init['date']}{fcst_init['hour']}"

    logger.info(f"Forecast Init Date: {fcst_init['monthday']} {year}", subtask = subtask)
    # Print Ensemble member and reference date+cycle time:
    txt_reftime = _print_reftime(fcst_init, ens_num)
    logger.info(f"Reference time: {txt_reftime}", subtask = subtask)

    # 6-hourly output dir:
    outdirs['outdir_6hourly'] = \
        f"{args['outdir']}/6-Hourly/" + \
        f"{fcst_init['monthday']}/{year}/ens{ens_num}"
    if not os.path.exists(outdirs['outdir_6hourly']):
        os.makedirs(outdirs['outdir_6hourly'], exist_ok=True)

    # Monthly output dir:
    outdirs['outdir_monthly'] = \
        f"{args['outdir']}/Monthly/" + \
        f"{fcst_init['monthday']}/{year}/ens{ens_num}"
    if not os.path.exists(outdirs['outdir_monthly']):
        os.makedirs(outdirs['outdir_monthly'], exist_ok=True)

    # Loop over CFSv2 variables:
    cfsv2 = []
    for varname in ["prate", "pressfc", "tmp2m", "dlwsfc", "dswsfc",
                    "q2m", "wnd10m"]:
        logger.info(f"CFSv2 variable: {varname}", subtask = subtask)
        subdir, file_pfx, file_sfx = \
            _set_input_file_info(fcst_init['year_cfsv2'],
                                 fcst_init['month'],
                                 varname)
        indir = f"{args['forcedir']}/{subdir}/"
        if subdir == "Oper_TS" and not os.path.exists(indir):
            indir = f"{args['forcedir']}/"

        indir += f"{fcst_init['year_cfsv2']}/{fcst_init['date']}"

        # Checking CFSv2 member-date subdir presence:
        if not os.path.isdir(indir):
            logger.error(f"CFSV2 directory, {indir}, does NOT exist.", subtask= subtask)
            sys.exit(1)  # Exit with an error code

        # Convert GRIB file to netCDF and handle missing/corrupted data
        cfsv2.append(read_wgrib (indir, file_pfx, fcst_init['timestring'], \
                                 file_sfx, outdirs['outdir_6hourly'], temp_name, varname, args['patchdir'],
                                 [logger,subtask]))

    cfsv2 = xr.merge (cfsv2, compat='override')
    _migrate_to_monthly_files(cfsv2.sel (step = (cfsv2['valid_time']  >= dt1) &
                                         (cfsv2['valid_time']  < dt2)),
                                        outdirs, fcst_init, args, rank, logger, subtask)
    cfsv2.close()
    del cfsv2
    gc.collect()
    return 

if __name__ == "__main__":
    try:
        rank = int(os.environ.get('SLURM_PROCID', '0'))
        size = int(os.environ.get('SLURM_NTASKS', '1'))
    except (ValueError, TypeError):
        rank = 0
        size = 1
    if size==1:
        args = _read_cmd_args()
        with open(sys.argv[5], 'r', encoding="utf-8") as file:
            config = yaml.safe_load(file)
            for rank in range(config["EXP"]["lead_months"]):
                _driver(rank)
    else:
        _driver(rank)
        
