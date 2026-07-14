''' GEOSv3 regridder '''
import os
import sys
import gc
from datetime import datetime
os.environ['OMP_NUM_THREADS'] = '8'
os.environ['OPENBLAS_NUM_THREADS'] = '8'
os.environ['MKL_NUM_THREADS'] = '8'
os.environ['VECLIB_MAXIMUM_THREADS'] = '8'
os.environ['NUMEXPR_NUM_THREADS'] = '8'
# pylint: disable=wrong-import-position
import dask
dask.config.set(scheduler='threads', num_workers=8)
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
import xesmf as xe
import yaml
from ghis2s.shared.utils import get_domain_info, load_ncdata, get_chunk_sizes, write_ncfile, pack_dataset_to_int16, add_packing_attributes
from ghis2s.bcsd.bcsd_library.bcsd_functions import apply_regridding_with_mask, add_fcorr_vars, apply_fcorr
from ghis2s.bcsd.bcsd_library.bcsd_functions import VarLimits as lim
from ghis2s.shared.logging_utils import TaskLogger

limits = lim()
# Internal functions
def read_geosv3_elevation(static_file, logger):
    ''' reads GEOS5 surface geopotential (PHIS) and returns the height at the surface '''
    lis_const_g = 9.80616
    phyis = load_ncdata(static_file, logger, var_name='PHIS').isel(time = 0)
    return phyis/lis_const_g

def magnitude(_a, _b):
    ''' computes wind magnitude u^2 + v^2'''
    return np.sqrt(_a**2 + _b**2)

def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} syr eyr fcst_init_monthday outdir"
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""

    with open(sys.argv[5], 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    if len(sys.argv) > 7:
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

def write_monthly_optimized(this_6h, file_6h, file_mon, logger):
    ''' writes regridded raw Monthly and 3-hourly files '''
    var_list = ["PRECTOT", "PS", "T2M", "LWGAB", "SWGDN", "QV2M", "WIND10M"]
    compress_encoding = {
        var: {'dtype': 'int16', "zlib": True, "complevel": 6, "shuffle": True, "_FillValue": -32767,
              "missing_value": -32767}
        for var in var_list
    }

    # ==========================================
    # PART 1: 3-HOURLY FILE
    # ==========================================
    packed_ds, packing_params = pack_dataset_to_int16(this_6h, var_list, logger=logger)

    for var in var_list:
        if var in packed_ds:
            packed_ds[var].attrs.pop('_FillValue', None)

    write_ncfile(packed_ds, file_6h, compress_encoding, logger)
    add_packing_attributes(file_6h, packing_params, logger)

    packed_ds.close()
    this_6h.close()
    del packed_ds, this_6h
    gc.collect()

    # ==========================================
    # PART 2: MONTHLY FILE
    # ==========================================
    ds_disk = load_ncdata(file_6h, logger, decode_cf=True,
                          decode_timedelta=True, chunks={'time': -1})
    this_mon = ds_disk.mean(dim="time", keep_attrs=True)
    this_mon = this_mon.fillna(-9999.0)

    # Pack
    packed_ds, packing_params = pack_dataset_to_int16(this_mon, var_list, logger=logger)

    # Strip from packed_ds
    for var in var_list:
        if var in packed_ds:
            packed_ds[var].attrs.pop('missing_value', None)
            packed_ds[var].attrs.pop('_FillValue', None)
            packed_ds[var].encoding.pop('missing_value', None)
            packed_ds[var].encoding.pop('_FillValue', None)

            packed_ds[var].attrs.pop('scale_factor', None)
            packed_ds[var].attrs.pop('add_offset', None)
            packed_ds[var].encoding.pop('scale_factor', None)
            packed_ds[var].encoding.pop('add_offset', None)

    # Write
    write_ncfile(packed_ds, file_mon, compress_encoding, logger)
    add_packing_attributes(file_mon, packing_params, logger)

    ds_disk.close()
    this_mon.close()
    packed_ds.close()
    del ds_disk, this_mon, packed_ds
    gc.collect()

def write_monthly_files(this_6h, file_6h, file_mon, logger):
    ''' writes regridded raw Monthly and 3-hourly files '''
    encoding = {
        var: {'dtype': 'float32', "zlib": True, "complevel": 6, "shuffle": True, "missing_value": -9999.}
        for var in ["PRECTOT", "PS", "T2M", "LWGAB", "SWGDN", "QV2M", "WIND10M"]
    }

    write_ncfile(this_6h, file_6h, encoding, logger)
    this_mon = this_6h.mean(dim="time")
    write_ncfile(this_mon, file_mon, encoding, logger)
    this_6h.close()
    this_mon.close()
    del this_6h, this_mon

def _migrate_to_monthly_files(geosv3_in, outdirs, fcst_init, args, rank, logger, subtask, static_cache):
    ''' regrids raw GEOSv3 '''
    regrid_method = {
        '25km': {'PRECTOT':'conservative', 'SWGDN':'bilinear', 'LWGAB':'bilinear','PS':'bilinear',
                 'QV2M':'bilinear', 'T2M':'bilinear', 'WIND10M':'bilinear'},
        '10km': {'PRECTOT':'conservative', 'SWGDN':'bilinear', 'LWGAB':'bilinear','PS':'bilinear',
                 'QV2M':'bilinear', 'T2M':'bilinear', 'WIND10M':'bilinear'},
        '5km': {'PRECTOT':'conservative', 'SWGDN':'conservative', 'LWGAB':'conservative',
                'PS':'conservative','QV2M':'conservative', 'T2M':'bilinear',
                'WIND10M':'bilinear'}}

    outdir_6hourly = outdirs["outdir_6hourly"]
    outdir_monthly = outdirs["outdir_monthly"]
    final_name_pfx = f"{fcst_init['monthday']}.geosv3."
    reftime = \
         f"{fcst_init['year']}-{fcst_init['month']}-{fcst_init['day']}"
    reftime += f",{fcst_init['hour']}:00:00,1hour"

    # resample to the S2S grid
    # build regridder
    weightdir = args["config"]['SETUP']['supplementarydir'] + '/bcsd_fcst/LDT_mask/'
    target_land_mask = xr.open_dataset(args["config"]['SETUP']['supplementarydir'] +
                                       '/lis_darun/' + args["config"]['SETUP']['ldtinputfile'])
    target_land_mask = target_land_mask.rename({'north_south': 'lat', 'east_west': 'lon'})
    lats, lons = get_domain_info(args["configfile"], coord=True)
    resol = round((lats[1] - lats[0])*100)
    resol = f'{resol}km'
    lat_chunk, lon_chunk = get_chunk_sizes(None, [len(lats), len(lons)])

    # check for optional keywords
    lapsrate_corr = args["config"].get('BCSD', {}).get('force_corr', {}).get('lapserate') is True
    aspect_corr = args["config"].get('BCSD', {}).get('force_corr', {}).get('aspect') is True

    n_outvar = 0
    if lapsrate_corr: n_outvar += 4
    if aspect_corr: n_outvar += 1

    if static_cache and "bil_regridder" in static_cache:
        logger.info("Using cached regridders and land masks...", subtask=subtask)
        target_land_mask = static_cache["target_land_mask"]
        geosv3_land_mask = static_cache["geosv3_land_mask"]
        bil_regridder = static_cache["bil_regridder"]
        con_regridder = static_cache["con_regridder"]
        ds_out_template = static_cache["ds_out_template"]
        if lapsrate_corr:
            geos_elev = static_cache["geos_elev"]
        geosv3_masked = geosv3_in.copy()
        geosv3_masked['mask'] = geosv3_land_mask['LANDMASK']
    else:
        logger.info("Building regridders and land masks for the first time...", subtask=subtask)
        weightdir = args["config"]['SETUP']['supplementarydir'] + '/bcsd_fcst/LDT_mask/'
        target_land_mask = xr.open_dataset(args["config"]['SETUP']['supplementarydir'] +
                                           '/lis_darun/' + args["config"]['SETUP']['ldtinputfile'])
        target_land_mask = target_land_mask.rename({'north_south': 'lat', 'east_west': 'lon'})
        geosv3_land_mask = xr.open_dataset(weightdir + f'GEOSv3_{resol}_landmask.nc4')

        if lapsrate_corr:
            static_file = args["config"]['SETUP']['supplementarydir'] + '/bcsd_fcst/geosv3_staticfields_20150401.nc'
            geos_elev = read_geosv3_elevation(static_file, [logger, subtask])
            geos_elev = xr.Dataset({'ELEVATION': geos_elev})

        ds_out_template = xr.Dataset({"lat": (["lat"], lats), "lon": (["lon"], lons)})

        geosv3_masked = geosv3_in.copy()
        geosv3_masked['mask'] = geosv3_land_mask['LANDMASK']

        weight_file_bil = weightdir + f'GEOSv3_{resol}_bilinear_land.nc'
        bil_regridder = xe.Regridder(geosv3_masked, ds_out_template, "bilinear", periodic=True,
                                     reuse_weights=True, extrap_method='nearest_s2d', filename=weight_file_bil)

        weight_file_con = weightdir + f'GEOSv3_{resol}_conservative_land.nc'
        con_regridder = xe.Regridder(geosv3_in, ds_out_template, "conservative", periodic=True,
                                     reuse_weights=True, filename=weight_file_con)

        # Save them to the cache for the next iteration!
        if static_cache is not None:
            static_cache["target_land_mask"] = target_land_mask
            static_cache["geosv3_land_mask"] = geosv3_land_mask
            static_cache["bil_regridder"] = bil_regridder
            static_cache["con_regridder"] = con_regridder
            static_cache["ds_out_template"] = ds_out_template
            if lapsrate_corr:
                static_cache["geos_elev"] = geos_elev

    # --- Use saved template
    ds_out = ds_out_template.copy()
    '''
    resol='25km': NY=720, NX=1440; 
    resol='10km': NY=1500, NX=3600; 
    resol='5km': NY=3000, NX=7200; 
    geosv3 dimensions: step: 248, latitude: 361, longitude: 720
    '''

    method = regrid_method.get(resol)
    bilinear_vars = [var for var in geosv3_masked.data_vars
                     if var in method and method[var] == 'bilinear']
    conservative_vars = [var for var in geosv3_in.data_vars
                         if var in method and method[var] == 'conservative']

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

    if lapsrate_corr or aspect_corr:
        # Apply forcing corrections
        # Regrid elevation and add ELEVATION Difference for lapserate correction
        elevation_regridded = bil_regridder(geos_elev['ELEVATION'].where(
            geosv3_land_mask['LANDMASK'] > 0))
        elevation_regridded = elevation_regridded.where(target_land_mask['LANDMASK'] > 0)
        ds_out["ELEV_DIFF"] = target_land_mask["ELEVATION"] - elevation_regridded

        # Add local hour, doy, lat_2d for ASPECT correction
        ds_out = add_fcorr_vars(ds_out)

        # Apply forcing correction
        result = xr.apply_ufunc(
            apply_fcorr,
            target_land_mask["LANDMASK"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            target_land_mask["SLOPE"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            target_land_mask["ASPECT"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["ELEV_DIFF"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["LAT_2D"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["T2M"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["QV2M"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["LWGAB"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["PS"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["SWGDN"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["LHOUR"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            ds_out["DOY"].chunk({'lat': lat_chunk, 'lon': lon_chunk}),
            input_core_dims=[[],[],[],[],[], ['time'],['time'],['time'],['time'],
                             ['time'],['time'],['time']],
            kwargs={'force_corr': args["config"]["BCSD"]["force_corr"]},
            output_core_dims=[['variable', 'time']],
            vectorize=True,
            dask='parallelized',
            output_dtypes=[np.float32],
            dask_gufunc_kwargs={
                'output_sizes': {'time': ds_out.sizes['time'], 'variable': n_outvar},
                'allow_rechunk': True  
            }
        ).transpose('variable', 'time', 'lat', 'lon')

        # update ds_out
        var_no = 0
        if lapsrate_corr:
            ds_out["T2M"] = result.isel(variable = var_no + 0)
            ds_out["QV2M"] = result.isel(variable = var_no + 1)
            ds_out["LWGAB"] = result.isel(variable = var_no + 2)
            ds_out["PS"] = result.isel(variable = var_no + 3)
            var_no = 4
        if aspect_corr:
            ds_out["SWGDN"] = result.isel(variable = var_no)

        # remove temporary variables
        ds_out = ds_out.drop_vars(["ELEV_DIFF", "LHOUR", "DOY", "LAT_2D"])

    mmm = fcst_init['monthday'].split("0")[0].capitalize()
    dt1 = datetime.strptime(f'{mmm} 1 {fcst_init["year"]}', '%b %d %Y')

    # Clip array
    ds_out = ds_out.fillna(-9999.)
    var_configs = {
        'PRECTOT': {'precip': True},
        'PS': {},
        'T2M': {},
        'LWGAB': {},
        'SWGDN': {},
        'QV2M': {},
        'WIND10M': {'var_name': 'WIND'}
    }
    ds_out = limits.clip_forcing_variables(ds_out, var_configs)
    #ds_out = ds_out.chunk({'time':-1, 'lat': -1, 'lon': -1})
    dt1 = dt1 + relativedelta(months=rank)
    file_6h = outdir_6hourly + '/' + \
        final_name_pfx + f'{dt1.year:04d}{dt1.month:02d}.nc'
    file_mon = outdir_monthly + '/' + \
        final_name_pfx + f'{dt1.year:04d}{dt1.month:02d}.nc'

    write_monthly_optimized(ds_out, file_6h, file_mon, [logger, subtask])
    logger.info(f"Writing: 3h GEOSv3 file: {file_6h}", subtask = subtask)
    logger.info(f"Writing: monthly GEOSv3 file: {file_mon}", subtask = subtask)
    ds_out.close()
    geosv3_masked.close()
    del ds_out, geosv3_masked

    # Force garbage collection and purge Xarray's internal file cache
    try:
        xr.backends.file_manager.FILE_CACHE.clear()
    except AttributeError:
        pass
    gc.collect()
    logger.info(f"Done processing GEOSv3 forecast files for {rank}", subtask = subtask)

def driver(rank, logger_task=None, static_cache=None):
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
    dt0 = datetime.strptime(f'{mmm} 1 {fcst_init["year"]}', '%b %d %Y')
    dt1 = dt0 + relativedelta(months=rank)
    subtask = f'{dt1.year:04d}{dt1.month:02d}'

    if logger_task is None:
        task_name = os.environ.get('SCRIPT_NAME')
        logger = TaskLogger(task_name,
                            os.getcwd(),
                            f'bcsd/bcsd_library/process_geosv3_forcing.py processing '
                            f'GEOSv3 ens{ens_num:01d} {subtask}')
    else:
        logger = logger_task

    yyyymm_ic = f'{dt0.year:04d}{dt0.month:02d}'
    yyyymm_fc = f'{dt1.year:04d}{dt1.month:02d}'
    ens_geos = f'ens{ens_num:02d}/'
    ens_num = f'ens{ens_num:01d}/'
    fcst_init["month"] = int(dt0.month)
    fcst_init["day"] = int(dt0.day)
    fcst_init["hour"] = 0

    geos_file = args["forcedir"] + yyyymm_ic + '/' + ens_geos + f'geos_s2s_v3.{yyyymm_fc}.nc'
    rad_file = args["forcedir"] + yyyymm_ic + '/' + ens_geos + f'rad_geos_s2s_v3.{yyyymm_fc}.nc'
    logger.info(f"Reading GEOS file: {geos_file}", subtask = subtask)
    logger.info(f"Reading Radiation file: {rad_file}", subtask = subtask)
    geos_full = load_ncdata(geos_file, [logger, subtask], chunks={'time': 31})[
        ['PRECTOTCORR', 'T2M', 'PS', 'Q2M', 'U10M', 'V10M','time', 'lat','lon']]
    rad_full = load_ncdata(rad_file, [logger, subtask], chunks={'time': 31})[['LWS', 'SWGDWN']]
    # Ensure the same coordinates
    rad_full = rad_full.assign_coords({
        'lat': geos_full['lat'],
        'lon': geos_full['lon'],
        'time': geos_full['time']
    })

    # Create new dataset with renamed variables
    geos_ds = xr.Dataset()

    for var in ['PRECTOTCORR', 'T2M', 'PS', 'Q2M']:
        geos_ds[new_name[var]] = geos_full[var]

    geos_ds['time'] = geos_full['time']
    geos_ds['lat'] = geos_full['lat']
    geos_ds['lon'] = geos_full['lon']

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
        f"{fcst_init['monthday']}/{year}/{ens_num}"
    if not os.path.exists(outdirs['outdir_6hourly']):
        os.makedirs(outdirs['outdir_6hourly'], exist_ok=True)

    # Monthly output dir:
    outdirs['outdir_monthly'] = \
        f"{args['outdir']}/Monthly/" + \
        f"{fcst_init['monthday']}/{year}/{ens_num}"
    if not os.path.exists(outdirs['outdir_monthly']):
        os.makedirs(outdirs['outdir_monthly'], exist_ok=True)

    _migrate_to_monthly_files(geos_ds, outdirs, fcst_init, args, rank, logger, subtask,
                              static_cache)
    geos_ds.close()
    del geos_ds
    gc.collect()

if __name__ == "__main__":
    try:
        RANK = int(os.environ.get('SLURM_PROCID', '0'))
        SIZE = int(os.environ.get('SLURM_NTASKS', '1'))
    except (ValueError, TypeError):
        RANK = 0
        SIZE = 1

    # --- empty dictionary to hold regridders ---
    static_cache_ = {}

    if SIZE==1:
        with open(sys.argv[5], 'r', encoding="utf-8") as _file:
            _config = yaml.safe_load(_file)
        TASK_NAME = os.environ.get('SCRIPT_NAME')
        LOGGER = TaskLogger(TASK_NAME,
                            os.getcwd(),
                            f'bcsd/bcsd_library/process_geosv3_forcing.py '
                            f'processing GEOSv3 ens{sys.argv[2]}')
        loop = [0, _config["EXP"]["lead_months"]]
        if len(sys.argv) == 7:
            start_rank = int(sys.argv[6])
            loop = [start_rank, start_rank + 1]
        for _rank in range(loop[0], loop[1]):
            driver(_rank, logger_task=LOGGER, static_cache=static_cache_)
    else:
        driver(RANK, static_cache=static_cache_)
