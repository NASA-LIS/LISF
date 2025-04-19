from datetime import datetime, date
import glob
import os
import sys
import yaml
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr

# pylint: disable=import-error
from metricslib import sel_var, compute_anomaly, compute_sanomaly
# pylint: enable=import-error
# pylint: disable=consider-using-f-string

FCST_INIT_MON = int(sys.argv[1])
TARGET_YEAR = int(sys.argv[2])
NMME_MODEL = sys.argv[3]
CONFIGFILE = sys.argv[4]
BASEOUTDIR = sys.argv[5]
ANOM = sys.argv[6]

# Load CONFIG file
with open(CONFIGFILE, 'r', encoding="utf-8") as file:
    CONFIG = yaml.safe_load(file)
HYD_MODEL = CONFIG["EXP"]["lsmdir"]
DOMAIN_NAME = CONFIG["EXP"]["DOMAIN"]
CLIM_SYR = int(CONFIG["BCSD"]["clim_start_year"])
CLIM_EYR = int(CONFIG["BCSD"]["clim_end_year"])
BASEDIR = BASEOUTDIR + "/DYN_" + ANOM + "/"
METRIC_VARS = CONFIG["POST"]["metric_vars"]
WEEKLY_VARS = CONFIG["POST"]["weekly_vars"]
HINDCASTS = CONFIG["SETUP"]["E2ESDIR"] + '/hindcast/s2spost/' + '{:02d}/'.format(FCST_INIT_MON)
FORECASTS = "./s2spost/"

OUTDIR = BASEDIR + '/' + HYD_MODEL
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR, exist_ok=True)

OUTFILE_TEMPLATE = '{}/{}_{}_{}_init_weeks1-6_{:02d}_{:04d}.nc'
TARGET_INFILE_TEMPLATE = \
    '{}/{:04d}{:02d}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.{}_' \
    'PA.ALL_DD.{:04d}{:02d}01_DT.0000_FD.{}_DT.0000_DF.NC'

CLIM_INFILE_TEMPLATE = \
    '{}/????{:02d}/{}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.{}_' \
    'PA.ALL_DD.*{:02d}01_DT.0000_FD.*{:02d}{:02d}_DT.0000_DF.NC'

LEAD_WEEKS = 6
for var_name in WEEKLY_VARS:
    OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, NMME_MODEL, \
                                      var_name, ANOM, FCST_INIT_MON, TARGET_YEAR)

    CURRENTDATE = date(TARGET_YEAR, FCST_INIT_MON, 2)
    for lead in range(LEAD_WEEKS):
        fcast_list = []
        clim_list = []
        clim_std = []
        print(f"[INFO] Computing {var_name} forecast anomaly for week {lead}")
        for count_days in range(7):
            # processing climatology
            INFILE = glob.glob(CLIM_INFILE_TEMPLATE.format(HINDCASTS, \
                                                 FCST_INIT_MON, NMME_MODEL, \
                                                 NMME_MODEL.upper(), DOMAIN_NAME, \
                                                 FCST_INIT_MON, \
                                                 CURRENTDATE.month, CURRENTDATE.day))
            
            
            day_xr = xr.open_mfdataset(INFILE, combine='by_coords')
            sel_cim_data = day_xr.sel(time= \
                       (day_xr.coords['time.year'] >= \
                       CLIM_SYR) & (day_xr.coords['time.year'] <= \
                       CLIM_EYR))
            clim_list.append(sel_var(sel_cim_data, var_name, HYD_MODEL).mean(dim = ['time','ensemble'], skipna = True))
            clim_std.append(sel_var(sel_cim_data, var_name, HYD_MODEL))
            
            # reading forecast
            INFILE = TARGET_INFILE_TEMPLATE.format(FORECASTS, \
                                                   TARGET_YEAR, FCST_INIT_MON, NMME_MODEL, \
                                                   NMME_MODEL.upper(), DOMAIN_NAME, \
                                                   TARGET_YEAR, FCST_INIT_MON, \
                                                   CURRENTDATE.strftime("%Y%m%d"))
            print(f"[INFO] {INFILE}")
            fcast_list.append(INFILE)
            CURRENTDATE += relativedelta(days=1)
            
        weekly_clim = xr.concat(clim_list, dim='day').mean(dim='day')
        weekly_std = xr.concat(clim_std, dim='day').mean(dim='day')
        weekly_std = weekly_std.std(dim = ['time','ensemble'])

        fcst_xr = xr.open_mfdataset(fcast_list, combine='by_coords')
        fcst_da = sel_var(fcst_xr, var_name, HYD_MODEL).mean(dim = ['time'], skipna = True)
 
        # Step-3 loop through each grid cell and convert data into anomaly
        # Defining array to store anomaly data
        lat_count, lon_count, ens_count = \
            len(fcst_xr.coords['lat']), \
            len(fcst_xr.coords['lon']), \
            len(fcst_xr.coords['ensemble'])

        # 4 members in anomaly output and so on.
        if lead == 0:
            all_anom = np.ones((ens_count, LEAD_WEEKS, lat_count, lon_count))*-9999

        print('[INFO] Converting data into anomaly')
        if (not np.array_equal(weekly_clim.lat.values, fcst_da.lat.values)) or \
           (not np.array_equal(weekly_clim.lon.values, fcst_da.lon.values)):
            weekly_clim = weekly_clim.assign_coords({"lon": fcst_da.lon.values,
                                                         "lat": fcst_da.lat.values})
            weekly_std = weekly_std.assign_coords({"lon": fcst_da.lon.values,
                                                         "lat": fcst_da.lat.values})
            
        if ANOM == 'ANOM':
            this_anom = xr.apply_ufunc(
                compute_anomaly,
                fcst_da.chunk({"lat": "auto", "lon": "auto"}).compute(),
                weekly_clim.chunk({"lat": "auto", "lon": "auto"}).compute(),
                input_core_dims=[['ensemble',],[]],
                exclude_dims=set(('ensemble',)),
                output_core_dims=[['ensemble',]],
                vectorize=True,  # loop over non-core dims
                dask="forbidden",
                output_dtypes=[np.float64])

        if ANOM == 'SANOM':
            this_anom = xr.apply_ufunc(
                compute_sanomaly,
                fcst_da.chunk({"lat": "auto", "lon": "auto"}).compute(),
                weekly_clim.chunk({"lat": "auto", "lon": "auto"}).compute(),
                weekly_std.chunk({"lat": "auto", "lon": "auto"}).compute(),
                input_core_dims=[['ensemble',],[],[]],
                exclude_dims=set(('ensemble',)),
                output_core_dims=[['ensemble',]],
                vectorize=True,  # loop over non-core dims
                dask="forbidden",
                output_dtypes=[np.float64])
            
        for ens in range(ens_count):
            all_anom[ens, lead, :, :] = this_anom [:,:,ens]

    
    ### Step-4 Writing output file
    all_anom = np.ma.masked_array(all_anom, mask=(all_anom == -9999))

    ## Creating an latitude and longitude array based on locations of corners
    lats = np.arange(fcst_xr.attrs['SOUTH_WEST_CORNER_LAT'], \
                     fcst_xr.attrs['SOUTH_WEST_CORNER_LAT'] + \
                     (lat_count*0.25), 0.25)
    lons = np.arange(fcst_xr.attrs['SOUTH_WEST_CORNER_LON'], \
                     fcst_xr.attrs['SOUTH_WEST_CORNER_LON'] + \
                     (lon_count*0.25), 0.25)
    anom_xr = xr.Dataset()
    anom_xr['anom'] = (('ens', 'time', 'latitude', 'longitude'), all_anom)
    anom_xr.coords['latitude'] = (('latitude'), lats)
    anom_xr.coords['longitude'] = (('longitude'), lons)
    anom_xr.coords['time'] = (('time'), np.arange(0, LEAD_WEEKS, dtype=int))
    anom_xr.coords['ens'] = (('ens'), np.arange(0, ens_count, dtype=int))
    print(f"[INFO] Writing {OUTFILE}")
    anom_xr.to_netcdf(OUTFILE, format="NETCDF4",
                      encoding = {'anom': {"zlib":True, "complevel":6, "shuffle":True, "missing_value": -9999., "_FillValue": -9999.}})
    
        
        
