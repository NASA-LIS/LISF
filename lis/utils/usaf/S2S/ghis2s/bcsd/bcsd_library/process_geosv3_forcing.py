import os
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
import xesmf as xe
import yaml
import gc

# pylint: disable=import-error
from ghis2s.shared.utils import get_domain_info, load_ncdata
from ghis2s.bcsd.bcsd_library.bcsd_function import apply_regridding_with_mask
from ghis2s.bcsd.bcsd_library.bcsd_function import VarLimits as lim
from ghis2s.shared.logging_utils import TaskLogger
# pylint: enable=import-error
limits = lim()
# Internal functions

def magnitude(_a, _b):
    ''' computes wind magnitude u^2 + v^2'''
    func = lambda x, y: np.sqrt(x**2 + y**2)
    return xr.apply_ufunc(func, _a, _b)

def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} syr eyr fcst_init_monthday outdir"
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""

    with open(sys.argv[5], 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if len(sys.argv) != 6:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    args = {
        "syr" : sys.argv[1],
        "ens_num" : sys.argv[2],
        "fcst_init_monthday" : sys.argv[3],
        "outdir" : sys.argv[4],
        "forcedir" : config['BCSD']["fcst_download_dir"] + 'sfc_tavg_3hr_glo_L720x361_sfc/',
        "configfile": sys.argv[5],
        "config": config,
    }

    return args

def write_monthly_files(this_6h, file_6h, file_mon):
    encoding = {
        var: {"zlib": True, "complevel": 6, "shuffle": True, "missing_value": -9999.}
        for var in ["PRECTOT", "PS", "T2M", "LWGAB", "SWGDN", "QV2M", "WIND10M"]
    }

    this_6h.to_netcdf(file_6h, format="NETCDF4", encoding=encoding)
    this_mon = this_6h.mean(dim="time")
    this_mon.to_netcdf(file_mon, format="NETCDF4", encoding=encoding)
    this_6h.close()
    this_mon.close()
    del this_6h, this_mon
    return

def _migrate_to_monthly_files(geosv3_in, outdirs, fcst_init, args, rank, logger, subtask):
    regrid_method = {
        '25km': {'PRECTOT':'conservative', 'SWGDN':'bilinear', 'LWGAB':'bilinear','PS':'bilinear',
                 'QV2M':'bilinear', 'T2M':'bilinear', 'WIND10M':'bilinear'},
        '10km': {'PRECTOT':'conservative', 'SWGDN':'bilinear', 'LWGAB':'bilinear','PS':'bilinear',
                 'QV2M':'bilinear', 'T2M':'bilinear', 'WIND10M':'bilinear'},
        '5km': {'PRECTOT':'conservative', 'SWGDN':'conservative', 'LWGAB':'conservative','PS':'conservative',
                 'QV2M':'conservative', 'T2M':'bilinear', 'WIND10M':'bilinear'},}
    
    outdir_6hourly = outdirs["outdir_6hourly"]
    outdir_monthly = outdirs["outdir_monthly"]
    final_name_pfx = f"{fcst_init['monthday']}.geosv3."
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
    # read GEOSv3 land mask
    geosv3_land_mask = xr.open_dataset(weightdir + f'GEOSv3_{resol}_landmask.nc4')
    
    '''
    resol='25km': NY=720, NX=1440; 
    resol='10km': NY=1500, NX=3600; 
    resol='5km': NY=3000, NX=7200; 
    geosv3 dimensions: step: 248, latitude: 361, longitude: 720
    '''
    
    method = regrid_method.get(resol)
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], lats),
            "lon": (["lon"], lons),
        }
    )
    geosv3_masked = geosv3_in.copy()
    geosv3_masked['mask'] = geosv3_land_mask['LANDMASK']

    weight_file = weightdir + f'GEOSv3_{resol}_bilinear_land.nc'
    bil_regridder = xe.Regridder(geosv3_masked, ds_out, "bilinear", periodic=True, 
                                 reuse_weights=True,
                                 extrap_method='nearest_s2d', 
                                 filename=weight_file)

    weight_file = weightdir + f'GEOSv3_{resol}_conservative.nc'
    con_regridder = xe.Regridder(geosv3_in, ds_out, "conservative", periodic=True, 
                                 reuse_weights=True, 
                                 filename=weight_file)
    
    bilinear_vars = [var for var in geosv3_in.data_vars if var in method and method[var] == 'bilinear']
    conservative_vars = [var for var in geosv3_in.data_vars if var in method and method[var] == 'conservative']

    if bilinear_vars:
        geosv3_bilinear = geosv3_masked[bilinear_vars]
        result_bilinear = apply_regridding_with_mask(geosv3_bilinear, bil_regridder,
                                                     geosv3_land_mask, target_land_mask, 'bilinear')
        for var in bilinear_vars:
            ds_out[var] = result_bilinear[var]

    if conservative_vars:
        geosv3_conservative = geosv3_in[conservative_vars]
        result_conservative = apply_regridding_with_mask(geosv3_conservative, con_regridder,
                                                         geosv3_land_mask, None, 'conservative')
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
    ds_out["LWGAB"].values[:] = limits.clip_array(np.array(ds_out["LWGAB"].values[:]),
                                                 var_name = "LWGAB")
    ds_out["SWGDN"].values[:] = limits.clip_array(np.array(ds_out["SWGDN"].values[:]),
                                                   var_name = "SWGDN")
    ds_out["QV2M"].values[:] = limits.clip_array(np.array(ds_out["QV2M"].values[:]),
                                                 var_name = "QV2M")
    ds_out["WIND10M"].values[:] = limits.clip_array(np.array(ds_out["WIND10M"].values[:]),
                                                     var_name = "WIND")

    dt1 = dt1 + relativedelta(months=rank)
    file_6h = outdir_6hourly + '/' + \
        final_name_pfx + '{:04d}{:02d}.nc'.format (dt1.year,dt1.month)
    file_mon = outdir_monthly + '/' + \
        final_name_pfx + '{:04d}{:02d}.nc'.format (dt1.year,dt1.month)

    write_monthly_files(ds_out, file_6h, file_mon)
    logger.info(f"Writing: 3h GEOSv3 file: {file_6h}", subtask = subtask)
    logger.info(f"Writing: monthly GEOSv3 file: {file_mon}", subtask = subtask)
    ds_out.close()
    geosv3_masked.close()
    del ds_out, geosv3_masked
    logger.info(f"Done processing GEOSv3 forecast files for {rank}", subtask = subtask)
    return

def _driver(rank):
    """Main driver. rank = forecast month"""
    new_name = {'PRECTOTCORR': 'PRECTOT', 'T2M': 'T2M', 'PS': 'PS',
                'Q2M': 'QV2M', 'LWS': 'LWGAB', 'SWGDWN': 'SWGDN'}
    args = _read_cmd_args()
    fcst_init = {}
    fcst_init["monthday"] = args['fcst_init_monthday']
    outdirs = {}
    year = int(args['syr'])
    ens_num = int(args['ens_num'])

    fcst_init["year"] = year

    mmm = fcst_init['monthday'].split("0")[0].capitalize()
    dt0 = datetime.strptime('{} 1 {}'.format(mmm,fcst_init["year"]), '%b %d %Y')
    dt1 = dt0 + relativedelta(months=rank)
    dt2 = dt1 + relativedelta(months=1)
    task_name = os.environ.get('SCRIPT_NAME')
    subtask = f'{dt1.year:04d}{dt1.month:02d}'
    logger = TaskLogger(task_name,
                        os.getcwd(),
                        f'bcsd/process_forecast_data.py processing GEOSv3ens{ens_num:02d} for {subtask}')

    yyyymm_ic = f'{dt0.year:04d}{dt0.month:02d}'
    yyyymm_fc = f'{dt1.year:04d}{dt1.month:02d}'
    ensNN = f'/ens{ens_num:02d}/'
    fcst_init["month"] = int(dt0.month)
    fcst_init["day"] = int(dt0.day)
    fcst_init["hour"] = 0
    
    geos_file = args["forcedir"] + yyyymm_ic + ensNN + f'geos_s2s_v3.{yyyymm_fc}.nc'
    rad_file = args["forcedir"] + yyyymm_ic + ensNN + f'rad_geos_s2s_v3.{yyyymm_fc}.nc'
    logger.info(f"Reading GEOS file: {geos_file}", subtask = subtask)
    logger.info(f"Reading Radiation file: {rad_file}", subtask = subtask)
    geos_full = load_ncdata(geos_file, [logger, subtask])[['PRECTOTCORR', 'T2M', 'PS', 'Q2M', 'U10M', 'V10M','time','lat','lon']].load()
    geos_full =  geos_full.rename({'lat': 'latitude', 'lon': 'longitude'})
    rad_full = load_ncdata(rad_file, [logger, subtask])[['LWS', 'SWGDWN']].load()
    rad_full =  rad_full.rename({'lat': 'latitude', 'lon': 'longitude'})
    
    # Create new dataset with renamed variables
    geos_ds = xr.Dataset()
    
    for var in ['PRECTOTCORR', 'T2M', 'PS', 'Q2M']:
        geos_ds[new_name[var]] = geos_full[var]

    geos_ds['time'] = geos_full['time']
    geos_ds['latitude'] = geos_full['latitude']
    geos_ds['longitude'] = geos_full['longitude']
    
    # Calculate wind magnitude
    u10 = magnitude(geos_full['U10M'], geos_full['V10M'])
    geos_ds['WIND10M'] = u10
    geos_ds['WIND10M'].attrs = {
        'units': 'm/s',
        'short_name': 'wnd10m',
        'long_name': 'Wind Speed',
        'level': '10 m above ground'
    }
    
    # Add radiation variables
    for var in ['LWS', 'SWGDWN']:
        geos_ds[new_name[var]] = rad_full[var]
    
    # Close datasets to free memory
    geos_full.close()
    rad_full.close()
    del geos_full, rad_full
    gc.collect()

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

    _migrate_to_monthly_files(geos_ds, outdirs, fcst_init, args, rank, logger, subtask)
    geos_ds.close()
    del geos_ds
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
