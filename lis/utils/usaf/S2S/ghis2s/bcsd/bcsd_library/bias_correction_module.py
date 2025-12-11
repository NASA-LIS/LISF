#!/usr/bin/env python
"""
# Consolidated Shrad Shukla's bias_correction_modulefast.py and bias_correction_nmme_modulefast.py
#   dask paralalellized the process
# This module bias corrects a forecasts following
# probability mapping approach as described in Wood et al. 2002
"""

import calendar
import os
import sys
import gc
from datetime import datetime
from dateutil.relativedelta import relativedelta
import xarray as xr
import numpy as np
import yaml
from ghis2s.shared.utils import get_domain_info, load_ncdata, write_ncfile, get_chunk_sizes
from ghis2s.shared.logging_utils import TaskLogger
from ghis2s.bcsd.bcsd_library.bcsd_functions import calc_stats, lookup
from ghis2s.bcsd.bcsd_library.bcsd_functions import VarLimits as lim

limits = lim()
CMDARGS = str(sys.argv)
OBS_VAR = str(sys.argv[1])
FCST_VAR = str(sys.argv[2])
BC_VAR = str(sys.argv[3])
UNIT = str(sys.argv[4])
INIT_FCST_MON = int(sys.argv[5])
INIT_FCST_YEAR = int(sys.argv[6])
MODEL_NAME = str(sys.argv[7])
ENS_NUM = int(sys.argv[8])
CLIM_SYR = int(sys.argv[9])
CLIM_EYR = int(sys.argv[10])
CONFIG_FILE = str(sys.argv[11])
OBS_CLIM_INDIR = str(sys.argv[12])
FCST_DIR = str(sys.argv[13])
LEAD_NUM = int(sys.argv[14])
FCST_INIT_DATE = f'{calendar.month_abbr[INIT_FCST_MON].lower()}01'
NUM_YRS = (CLIM_EYR-CLIM_SYR)+1
TINY = ((1/(NUM_YRS))/ENS_NUM)/2
OUTDIR = f'{FCST_DIR}/bcsd/Monthly/{FCST_INIT_DATE}/'
SUBTASK = f'{FCST_VAR}_{LEAD_NUM+1:02d}'
os.makedirs(OUTDIR, exist_ok=True)

with open(CONFIG_FILE, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)
ldt_xr = xr.open_dataset(config['SETUP']['supplementarydir'] + '/lis_darun/' + \
        config['SETUP']['ldtinputfile'])
land_mask_ldt = ldt_xr.rename({'north_south': 'latitude', 'east_west': 'longitude'})['LANDMASK']
LATS, LONS = get_domain_info(CONFIG_FILE, coord=True)
land_mask_ldt = land_mask_ldt.assign_coords(latitude=LATS, longitude=LONS)

LAT_CHUNK, LON_CHUNK = get_chunk_sizes(land_mask_ldt)

METFORCE = False
if MODEL_NAME == "METF":
    METFORCE = True

OBS_CLIM_FILE_TEMPLATE = '{}/{}_obs_clim.nc'
OUTFILE_TEMPLATE = '{}/{}.{}.{}_{:04d}{:02d}.nc'
if METFORCE:
    FCST_CLIM_INFILE = f'{FCST_DIR}/raw/Climatology/{FCST_INIT_DATE}/{FCST_VAR}_fcst_clim.nc'
    if FCST_VAR in ['LWGAB', 'SWGDN', 'QV2M']:
        alternate_name = {
            'LWGAB': 'LWS',
            'SWGDN': 'SLRSF',
            'QV2M': 'Q2M',
        }
        if not os.path.exists(FCST_CLIM_INFILE):
            FCST_CLIM_INFILE = f'{FCST_DIR}/raw/Climatology/{FCST_INIT_DATE}/'\
            f'{alternate_name[FCST_VAR]}_fcst_clim.nc'
    FCST_INFILE_TEMPLATE = '{}/raw/Monthly/{}/{:04d}/ens{:01d}/{}.{}.{:04d}{:02d}.nc'
    FCST_MODEL = config['BCSD']['metforce_source']
else:
    FCST_CLIM_INFILE = \
        f'{FCST_DIR}/raw/Climatology/{FCST_INIT_DATE}/{MODEL_NAME}/{FCST_VAR}_fcst_clim.nc'
    FCST_INFILE_TEMPLATE = '{}/raw/Monthly/{}/{}/{:04d}/ens{:01d}/{}.nmme.monthly.{:04d}{:02d}.nc'
    FCST_MODEL = MODEL_NAME

CF2VAR = {
    'PRECTOT': 'PRECTOT',
    'LWGAB': 'LWGAB',
    'SWGDN': 'SWGDN',
    'PS': 'PS',
    'QV2M':'QV2M',
    'T2M': 'T2M',
    'U10M': 'WIND',
    'WIND10M': 'WIND',
    }

def align_coordinates(arrays_dict):
    """Ensure all arrays have identical coordinate values"""    

    # Use LAND_MASK as reference coordinates
    reference = arrays_dict['LAND_MASK']
    ref_lat = reference.latitude
    ref_lon = reference.longitude

    aligned_arrays = {}
    for name, arr in arrays_dict.items():
        if name == 'LAND_MASK':
            aligned_arrays[name] = \
                arr.astype(np.float32).chunk({'latitude': LAT_CHUNK, 'longitude': LON_CHUNK})
            continue

        # Check if reindexing is needed
        lat_diff = np.abs(arr.latitude.values - ref_lat.values).max()
        lon_diff = np.abs(arr.longitude.values - ref_lon.values).max()

        if lat_diff > 1e-10 or lon_diff > 1e-10:
            aligned_arrays[name] = \
                arr.reindex(latitude=ref_lat,
                            longitude=ref_lon,
                            method='nearest')\
                   .astype(np.float32).chunk({'latitude': LAT_CHUNK, 'longitude': LON_CHUNK})
        else:
            aligned_arrays[name] = \
                arr.astype(np.float32).chunk({'latitude': LAT_CHUNK, 'longitude': LON_CHUNK})

    return aligned_arrays


def process_point(land_mask, obs_clim, fcst_clim, fcst_monthly):
    ''' function processes a land point '''
    correct_fcst_coarse = np.ones((ENS_NUM), dtype=np.float32)*-9999.
    if land_mask > 0.:
        # Reading observed and forcast time series; Note the 1st column is quantile time series
        obs_quant_ts, obs_clim_ts = obs_clim[0, :], obs_clim[1, :]
        fcst_quant_ts, fcst_clim_ts = fcst_clim[0, :], fcst_clim[1, :]

        # Calculating mean, standard deviation and skew of observed and forecast time series
        obs_mean, obs_sd, obs_skew = calc_stats(obs_clim_ts, TINY)
        fcst_mean, fcst_sd, fcst_skew = calc_stats(fcst_clim_ts, TINY)

        # Note bias correction is done seprately for each ensemble member
        for ens_num in range (0, ENS_NUM):
            target_fcst_val = fcst_monthly[ens_num]
            # First determine the quantile for given target forecast value
            target_fcst_quant = lookup(target_fcst_val, fcst_clim_ts, fcst_quant_ts, \
                                       len(fcst_clim_ts), BC_VAR, 'QUAN', fcst_mean, fcst_sd, \
                                       fcst_skew, TINY)

            # Also note that QUAN helps the the function lookup determine if we are trying to
            # convert a value to quantile or VICE versa. For converting a value to quantile
            # use 'QUAN'. For converting quantile to value use 'DATA'

            ## Now using quantile above, determine the value from observed climatology
            bias_corrected_value = lookup(target_fcst_quant, obs_quant_ts, obs_clim_ts, \
                                          len(obs_clim_ts), BC_VAR, 'DATA', obs_mean, obs_sd, \
                                          obs_skew, TINY)
            if METFORCE:
                correct_fcst_coarse[ens_num] = \
                    limits.clip_array(
                        np.atleast_1d(bias_corrected_value), var_name=CF2VAR.get(FCST_VAR)).flat[0]
            else:
                correct_fcst_coarse[ens_num] = \
                    limits.clip_array(np.atleast_1d(bias_corrected_value), var_name="PRECTOT", \
                                      max_val=0.004, precip=True).flat[0]
    return correct_fcst_coarse

def process_month(logger):
    ''' process the month '''
    logger.info(f"Processing month {LEAD_NUM}", subtask=SUBTASK)
    fcst_date = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1, 6) + \
        relativedelta(months=LEAD_NUM)
    fcst_year, fcst_month = fcst_date.year, fcst_date.month
    obs_clim_file = OBS_CLIM_FILE_TEMPLATE.format(OBS_CLIM_INDIR, OBS_VAR)
    logger.info(f"Reading observed climatology file: {obs_clim_file}", subtask=SUBTASK)
    obs_clim_full = load_ncdata(obs_clim_file, [logger, SUBTASK], var_name='clim')
    obs_clim_array = obs_clim_full.isel(DIST=[0, fcst_month]).rename({'time': 'time_obs'}) \
                                                             .assign_coords(DIST=[0, 1]) \
                                                             .astype(np.float32)

    logger.info(f"Reading forecast climatology file: {FCST_CLIM_INFILE}", subtask=SUBTASK)
    fcst_clim_full = load_ncdata(FCST_CLIM_INFILE, [logger, SUBTASK], var_name='clim',
                                 chunks={'DIST': 1, 'time': 50})
    fcst_clim_array = fcst_clim_full.isel(DIST=[0, LEAD_NUM+1]) \
                                .assign_coords(DIST=[0, 1]) \
                                .astype(np.float32)

    input_ens = []
    for ens in range(ENS_NUM):
        fcst_date = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1) + \
            relativedelta(months=LEAD_NUM)
        if METFORCE:
            infile = FCST_INFILE_TEMPLATE.format(FCST_DIR, FCST_INIT_DATE, \
                                                 INIT_FCST_YEAR, ens+1, FCST_INIT_DATE,
                                                 FCST_MODEL.lower(), fcst_year, fcst_month)
        else:
            infile = FCST_INFILE_TEMPLATE.format(FCST_DIR, FCST_INIT_DATE, MODEL_NAME, \
                                                 INIT_FCST_YEAR, ens+1, FCST_INIT_DATE,
                                                 fcst_year, fcst_month)

        logger.info(f"Reading forecast file: {infile}", subtask=SUBTASK)
        input_ens.append(load_ncdata(infile, [logger, SUBTASK], var_name=FCST_VAR)\
                         .astype(np.float32))

    fcst_monthly = xr.concat(input_ens, dim='Ens').rename({'lat': 'latitude', 'lon': 'longitude'})
    if 'time' in fcst_monthly.dims:
        # drop time dimension
        fcst_monthly = fcst_monthly.squeeze('time')

    #aligned_arrays = align_coordinates({
    #aligned_arrays = {
    #    'LAND_MASK': land_mask,
    #    'OBS_CLIM_ARRAY': obs_clim_array,
    #    'FCST_CLIM_ARRAY': fcst_clim_array,
    #    'FCST_MONTHLY': fcst_monthly}
    # Dask paralellized bias correction
    result = xr.apply_ufunc(
        process_point,
        land_mask_ldt.chunk({'latitude': LAT_CHUNK, 'longitude': LON_CHUNK}),
        obs_clim_array.chunk({'latitude': LAT_CHUNK, 'longitude': LON_CHUNK}),
        fcst_clim_array.chunk({'latitude': LAT_CHUNK, 'longitude': LON_CHUNK}),
        fcst_monthly.chunk({'latitude': LAT_CHUNK, 'longitude': LON_CHUNK}),
        input_core_dims=[[],['DIST','time_obs'],['DIST','time'],['Ens']],
        output_core_dims=[['Ens']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[np.float32],
        dask_gufunc_kwargs={
            'output_sizes': {'Ens': ENS_NUM},
            'allow_rechunk': True  
        }
    ).transpose('Ens', 'latitude', 'longitude')

    # write file
    outfile = OUTFILE_TEMPLATE.format(OUTDIR, FCST_VAR, FCST_MODEL, FCST_INIT_DATE, \
                                      fcst_year, fcst_month)
    out_xr = xr.Dataset(
        {
            FCST_VAR: (['Ens', 'latitude', 'longitude'], result.data)
        },
        coords={
            'Ens': np.arange(ENS_NUM),
            'latitude': LATS,
            'longitude': LONS
        }
    )

    encoding = {
        FCST_VAR: {
            'dtype': 'float32',
            'zlib': True,
            'complevel': 6,
            'shuffle': True,
            '_FillValue': -9999.0,
            'chunksizes': (ENS_NUM, LAT_CHUNK, LON_CHUNK) 
        }
    }
    del result, obs_clim_array, fcst_clim_array, fcst_monthly
    gc.collect()
    logger.info(f"Writing outpufile: {outfile}", subtask=SUBTASK)
    write_ncfile(out_xr, outfile, encoding, [logger, SUBTASK])

if __name__ == '__main__':
    # logging
    task_name = os.environ.get('SCRIPT_NAME')
    LOGGER = TaskLogger(task_name,
                        os.getcwd(),
                        'bcsd/bcsd_library/bias_correction_module.py processing ' +\
                        f'{MODEL_NAME} {FCST_VAR} for {FCST_INIT_DATE}')

    process_month(LOGGER)
