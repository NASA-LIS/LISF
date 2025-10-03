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
    
    geos_file = args["forcedir"] + yyyymm_ic + ensNN + f'geos_s2s_v3.{yyyymm_fc}.nc'
    
    rad_file = args["forcedir"] + yyyymm_ic + ensNN + f'rad_geos_s2s_v3.{yyyymm_fc}.nc'
    
    print(rank, ' ',geos_file)
    sys.exit()

    fcst_init['date'] = f"{fcst_init['year']}{fcst_init['monthday']}"
    fcst_init['month'] = int(fcst_init['monthday'][0:2])
    fcst_init['day'] = int(fcst_init['monthday'][2:4])
    fcst_init['timestring'] = f"{fcst_init['date']}{fcst_init['hour']}"
    print(rank, fcst_init)
    sys.exit()
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

    # Load GEOSv3 variables:
    
    
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
