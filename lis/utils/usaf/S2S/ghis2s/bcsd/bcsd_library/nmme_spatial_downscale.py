import os
import sys
from time import ctime as t_ctime
from time import time as t_time
from datetime import datetime
import numpy as np
import xarray as xr
import xesmf as xe
import yaml
# pylint: disable=import-error
from ghis2s.shared.utils import get_domain_info, load_ncdata, write_ncfile
from ghis2s.bcsd.bcsd_library.bcsd_function import VarLimits as lim
from ghis2s.bcsd.bcsd_library.bcsd_function import apply_regridding_with_mask
from ghis2s.shared.logging_utils import TaskLogger
# pylint: enable=import-error

limits = lim()

CMDARGS = str(sys.argv)
CMN = int(sys.argv[1])
MM = CMN - 1

if len(sys.argv) >= 6:  
    CYR = int(sys.argv[2])
    NMME_OUTPUT_DIR = str(sys.argv[3])
    NMME_MODEL = str(sys.argv[4])
    CONFIGFILE = str(sys.argv[5])
else:  
    NMME_OUTPUT_DIR = str(sys.argv[2])
    NMME_MODEL = str(sys.argv[3])
    CONFIGFILE = str(sys.argv[4])

with open(CONFIGFILE, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)

LEAD_MONS = config['EXP']['lead_months']
DATATYPE = config['SETUP']['DATATYPE'] 
SUPPLEMENTARY_DIR = config['SETUP']['supplementarydir'] + '/bcsd_fcst/'
ensemble_sizes = config['EXP']['ensemble_sizes'][0]
ENS_NUM = ensemble_sizes[NMME_MODEL]
NMME_DOWNLOAD_DIR = config['BCSD']['nmme_download_dir']

# Set up variables based on DATATYPE
if DATATYPE == 'hindcast':
    YEAR0 = 1982
    year_begin = YEAR0
    year_end = YEAR0 + 40
    LEAD_MONS = 9
    INFILE_TEMP = '{}/{}/prec.{}.mon_{}_{:04d}_{:04d}.nc'
    log_msg = f'bcsd/bcsd_library/nmme_spatial_downscale.py processing {NMME_MODEL} hindcast for month {CMN:02d}'
    
else:
    year_begin = CYR
    year_end = CYR+1
    INFILE_TEMP = '{}/{}/prec.{}.mon_{}.{:04d}.nc'
    log_msg = f'bcsd/bcsd_library/nmme_spatial_downscale.py processing {NMME_MODEL} for {CYR:04d}{CMN:02d} forecast'

OUTDIR_TEMPLATE = '{}/{}/{}/{:04d}/ens{}/'
OUTFILE_TEMPLATE = '{}/{}.nmme.monthly.{:04d}{:02d}.nc'

if not os.path.exists(NMME_OUTPUT_DIR):
    os.makedirs(NMME_OUTPUT_DIR, exist_ok=True)

nmme_path_dict = {
    'CFSv2': ['NCEP-CFSv2','NCEP-CFSv2'],
    'GEOSv2': ['NASA-GEOSS2S','NASA-GEOSS2S'],
    'CCM4': ['CanSIPS-IC3','CanSIPS-IC3'],
    'GNEMO5': ['CanSIPS-IC3','CanSIPS-IC3'],
    'CanESM5': ['CanSIPS-IC4','CanESM5'],
    'GNEMO52': ['CanSIPS-IC4', 'GEM5.2-NEMO'],
    'CCSM4': ['COLA-RSMAS-CCSM4', 'COLA-RSMAS-CCSM4'],
    'CESM1': ['COLA-RSMAS-CESM1', 'COLA-RSMAS-CESM1'],
    'GFDL': ['GFDL-SPEAR', 'GFDL-SPEAR'],
}

fcast_ens_index = {
    'CFSv2': [0, ENS_NUM],
    'GEOSv2': [0, ENS_NUM],
    'CCM4': [10, 20],
    'GNEMO5': [0, ENS_NUM],
    'CanESM5': [0, ENS_NUM],
    'GNEMO52': [0, ENS_NUM],
    'CCSM4': [0, ENS_NUM],
    'CESM1': [0, ENS_NUM],
    'GFDL': [0, ENS_NUM],
}

hcast_p1 = {
    'CFSv2': [1982, 2010],
    'GEOSv2': [0, ENS_NUM],
    'CCM4': [1991, 2020],
    'GNEMO5': [1991, 2020],
    'CanESM5': [1991, 2020],
    'GNEMO52': [1991, 2020],
    'CCSM4': [1982, 2021],
    'CESM1': [1991, 2021],
    'GFDL': [1991, 2020],
}   
    
MON = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 
       'Sep', 'Oct', 'Nov', 'Dec']
MONTH = ['jan01', 'feb01', 'mar01', 'apr01', 'may01', 'jun01', 'jul01',
          'aug01', 'sep01', 'oct01', 'nov01', 'dec01']
MONTHN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

# Logger
task_name = os.environ.get('SCRIPT_NAME')
subtask = f'{NMME_MODEL}'
logger = TaskLogger(task_name, os.getcwd(), log_msg)

LEADS1 = np.zeros((12, 12), dtype=int)
LDYR = np.zeros((12, 12), dtype=int)

for i in range(0, 12):
    for j in range(0, 12):
        k = MONTHN[i]+j
        KY = 0
        if k >= 13:
            k = k-12
            KY = 1
        LEADS1[i, j] = k
        LDYR[i, j] = KY

# Read domain info and land mask
LATS, LONS = get_domain_info(CONFIGFILE, coord=True)
resol = round((LATS[1] - LATS[0])*100)
resol = f'{resol}km'
weightdir = config['SETUP']['supplementarydir'] + '/bcsd_fcst/'
land_mask = load_ncdata(weightdir + f'NMME_{resol}_landmask.nc4', [logger, subtask])

if DATATYPE == 'forecast':
    nmme_path = nmme_path_dict[NMME_MODEL]
    ens_index = fcast_ens_index[NMME_MODEL]
    INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, nmme_path[0], nmme_path[1], MON[MM], CYR)
    
    logger.info(f"Reading: {INFILE}", subtask=subtask)
    nmme_da = load_ncdata(INFILE.strip(), [logger, subtask], var_name="prec", **dict(decode_times=False))
    
    if NMME_MODEL in ['CCM4', 'GNEMO5', 'CanESM5', 'GNEMO52']:
        nmme_da = nmme_da.transpose('S', 'L', 'M', 'Y', 'X')

    XPREC = np.array(nmme_da.values[:,0:LEAD_MONS,ens_index[0]:ens_index[1],:,:])
    
else:  # hindcast mode 
    XPREC = np.empty([40, LEAD_MONS, ENS_NUM, 181, 360])

    if NMME_MODEL in ['CFSv2', 'CCM4', 'GNEMO5', 'CanESM5', 'GNEMO52', 'CCSM4', 'CESM1', 'GFDL']:
        nmme_path = nmme_path_dict[NMME_MODEL]
        ens_index = fcast_ens_index[NMME_MODEL]
        p1 = hcast_p1[NMME_MODEL]
        y1 = p1[0] - YEAR0
        y2 = p1[1] - YEAR0 + 1
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, nmme_path[0], nmme_path[1], MON[MM], p1[0], p1[1])
        logger.info(f"Reading: {INFILE}", subtask=subtask)
        nmme_da = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
        if NMME_MODEL in ['CCM4', 'GNEMO5', 'CanESM5', 'GNEMO52']:
            nmme_da = nmme_da.transpose('S', 'L', 'M', 'Y', 'X')
        XPREC[y1:y2,:,:,:,:] = np.array(nmme_da.values[:, 0:LEAD_MONS, ens_index[0]:ens_index[1], :, :])
    
    if NMME_MODEL == 'CFSv2':
        if MON[MM] == 'Jan' or MON[MM] == 'Feb':
            SYR2 = 2011
            EYR2 = 2011
            INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
            logger.info(f"Reading: {INFILE}", subtask=subtask)
            nmme_da2 = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
            XPREC[29,:,:,:,:] = np.array(nmme_da2.values[:, 0:LEAD_MONS, 0:ENS_NUM, :, :])

            SYR3 = 2011
            EYR3 = 2021
            INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR3, EYR3)
            logger.info(f"Reading: {INFILE}", subtask=subtask)
            nmme_da3 = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
            XPREC[30:40,:,:,:,:] = np.array(nmme_da3.values[:, 0:LEAD_MONS, 0:ENS_NUM, :, :])
        else:
            SYR2 = 2011
            EYR2 = 2021
            INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
            logger.info(f"Reading: {INFILE}", subtask=subtask)
            nmme_da2 = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
            XPREC[29:40,:,:,:,:] = np.array(nmme_da2.values[:, 0:LEAD_MONS, 0:ENS_NUM, :, :])

    elif NMME_MODEL == 'GEOSv2':
        MODEL = 'NASA-GEOSS2S'
        if MM == 0:
            SYR1 = 1982
            EYR1 = 2017
            INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
            logger.info(f"Reading: {INFILE}", subtask=subtask)
            nmme_da = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
            XPREC[0:36,:,:,:,:] = np.array(nmme_da.values[:, 0:LEAD_MONS, 0:ENS_NUM, :, :])

            SYR2 = 2018
            EYR2 = 2021
            INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
            logger.info(f"Reading: {INFILE}", subtask=subtask)
            nmme_da2 = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
            XPREC[36:40,:,:,:,:] = np.array(nmme_da2.values[:, 0:LEAD_MONS, 0:ENS_NUM, :, :])
        else:
            SYR1 = 1982
            EYR1 = 2016
            INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR1, EYR1)
            logger.info(f"Reading: {INFILE}", subtask=subtask)
            nmme_da = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
            XPREC[0:35,:,:,:,:] = np.array(nmme_da.values[:, 0:LEAD_MONS, 0:ENS_NUM, :, :])

            SYR2 = 2017
            EYR2 = 2021
            INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
            logger.info(f"Reading: {INFILE}", subtask=subtask)
            nmme_da2 = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
            XPREC[35:40,:,:,:,:] = np.array(nmme_da2.values[:, 0:LEAD_MONS, 0:ENS_NUM, :, :])

    elif NMME_MODEL == 'GFDL':
        SYR2 = 2021
        EYR2 = 2021
        INFILE = INFILE_TEMP.format(NMME_DOWNLOAD_DIR, MODEL, MODEL, MON[MM], SYR2, EYR2)
        logger.info(f"Reading: {INFILE}", subtask=subtask)
        nmme_da2 = load_ncdata(INFILE, [logger, subtask], var_name="prec", **dict(decode_times=False))
        XPREC[39:40,:,:,:,:] = np.array(nmme_da2.values[:, 0:LEAD_MONS, 0:ENS_NUM, :, :])
    else:
        logger.error(f"Invalid argument for NMME_MODEL! Received {NMME_MODEL}", subtask=subtask)
        sys.exit(1)

LATI = np.array(nmme_da.Y)
LONI = np.array(nmme_da.X)

# convert mm/day to kg/m^2/s
XPREC = XPREC/86400. 

# Create xarray datasets for regridding
ds_in = xr.Dataset({
    "lat": (["lat"], LATI),
    "lon": (["lon"], LONI),
})

if DATATYPE == 'forecast':
    ds_in["XPREC"] = xr.DataArray(
        data=np.array(XPREC[0,0:LEAD_MONS,:,:,:]),
        dims=["mon","ens", "lat", "lon"],
        coords=dict(
            mon=(["mon"], np.arange(LEAD_MONS)),
            ens=(["ens"], np.arange(ENS_NUM)),
            lat=(["lat"], LATI),
            lon=(["lon"], LONI))
    )
else:  # hindcast
    ds_in["XPREC"] = xr.DataArray(
        data=np.array(XPREC[:,0:LEAD_MONS,:,:,:]),
        dims=["year", "mon","ens", "lat", "lon"],
        coords=dict(
            year=(["year"], np.arange(40)),
            mon=(["mon"], np.arange(LEAD_MONS)),
            ens=(["ens"], np.arange(ENS_NUM)),
            lat=(["lat"], LATI),
            lon=(["lon"], LONI))
    )

ds_out = xr.Dataset({
    "lat": (["lat"], LATS),
    "lon": (["lon"], LONS),
})

# now regridding
weight_file = weightdir + f'NMME_{resol}_conservative.nc'
regridder = xe.Regridder(ds_in, ds_out, "conservative", periodic=True,
                         reuse_weights=True, filename=weight_file)

ds_in2 = ds_in.rename_dims({"lat":"latitude", "lon":"longitude"})
ds_in = ds_in2.rename_vars({"lat":"latitude", "lon":"longitude"})
result = apply_regridding_with_mask(ds_in, regridder, land_mask)
ds_out["XPREC"] = result["XPREC"]

# Output writing
encoding = {'PRECTOT': {'dtype': 'float32', 'zlib': True, '_FillValue': -9999.}}
XPRECI = np.empty([1, LATS.size, LONS.size])
YR = year_begin - 1  

for y in range(year_end - year_begin):
    YR = YR + 1
    for m in range(0, ENS_NUM):
        for l in range(0, LEAD_MONS):
            if DATATYPE == 'forecast':
                XPRECI[0, :, :] = ds_out["XPREC"].values[l,m,:,:]
            else:  # hindcast
                XPRECI[0, :, :] = ds_out["XPREC"].values[y,l,m,:,:]
            
            jy = YR + LDYR[MM, l]
            l1 = LEADS1[MM, l]
            OUTDIR = OUTDIR_TEMPLATE.format(NMME_OUTPUT_DIR, MONTH[MM], NMME_MODEL, YR, m+1)
            OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, MONTH[MM], jy, l1)

            if not os.path.exists(OUTDIR):
                os.makedirs(OUTDIR)

            XPRECI = np.nan_to_num(XPRECI, nan=-9999.)
            XPRECI = limits.clip_array(XPRECI, var_name="PRECTOT", max_val=0.004, precip=True)

            SDATE = datetime(YR, CMN, 1)
            string_date = datetime.strftime(SDATE, "%Y-%m-%d")
            out_xr = xr.Dataset({
                'PRECTOT': xr.DataArray(
                    data=XPRECI,
                    dims=['time', 'lat', 'lon'],
                    coords={
                        'time': xr.DataArray([0], dims=['time'],
                                             attrs={'units': f'days since {string_date}',
                                                    'calendar': 'gregorian'}),
                        'lat': xr.DataArray(LATS, dims=['lat'],
                                            attrs={'units': 'degrees_north'}),
                        'lon': xr.DataArray(LONS, dims=['lon'],
                                            attrs={'units': 'degrees_east'})},
                    attrs={'units': 'kg m-2 s-1'})},
                            attrs={
                                'description': f'Downscaled to {resol}',
                                'history': f'Created {t_ctime(t_time())}',
                                'source': 'Raw NMME at 1deg'})

            logger.info(f"Writing: {OUTFILE}", subtask=subtask)
            write_ncfile(out_xr, OUTFILE, encoding, [logger, subtask])
