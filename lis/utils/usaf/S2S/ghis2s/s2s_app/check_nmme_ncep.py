#!/usr/bin/env python

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

'''
    This is good to plot NMME monthly precip anomaly with that of NOAA-NMME anom
'''
import os
import sys
import calendar
import argparse
from datetime import datetime as dt
import xarray as xr
import xesmf as xe
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4
# pylint: enable=no-name-in-module
import numpy as np
import yaml
#pylint: disable=import-error

NOAA_NCEP = {
    'CCM4': 'CanCM4i',
    'CCSM4':'NCAR_CCSM4',
    'CFSv2':'CFSv2',
    'GEOSv2':'NASA_GEOS5v2',
    'GFDL':'GFDL_SPEAR',
    'GNEMO5':'GEM5_NEMO',}

NOAA_NMME_PATH = '/discover/nobackup/projects/usaf_lis/GHI_S2S/NOAA_NMME/netcdf/realtime_anom/'
NOAA_NMME_TEMPLATE = '{}/{}/{}.prate.{:04d}{:02d}.anom.nc'

NAFPA_PATH = '/discover/nobackup/projects/ghilis/S2S/GLOBAL/DA_Run_Hist/output/SURFACEMODEL/'
NAFPA_FILE_TEMPLATE = '{}/{:04d}{:02d}/LIS_HIST_{:04d}{:02d}010000.d01.nc'

NMME_CLIM_TEMPLATE_RAW = '{}/raw/Climatology/{}/{}/PRECTOT_fcst_clim.nc'
NMME_MONTHLY_TEMPLATE_RAW = '{}/raw/Monthly/{}/{}/{:04d}/{}/{}.nmme.monthly.{:04d}{:02d}.nc'

NMME_CLIM_TEMPLATE_BCSD = 'hindcast/bcsd_fcst/NMME/bcsd/Monthly/{}/PRECTOT.{}.{}_*_*.nc'
NMME_MONTHLY_TEMPLATE_BCSD = '{}/bcsd/Monthly/{}/PRECTOT.{}.{}_{:04d}_{:04d}.nc'

rainf_levels = [0,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,7,9,11,13,15,17,20,25,30,50]

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--init_month', required=True, help='init month 1-12')
    parser.add_argument('-l', '--lead_month', required=True, help='lead month 0-9')
    parser.add_argument('-y', '--year', required=True, help='forcing start year')
    parser.add_argument('-c', '--config_file', required=True, help='config file')
    parser.add_argument('-m', '--model', required=True,
                        help='NMME model [CCM4, CCSM4, CFSv2, GEOSv2, GFDL, GNEMO5]')

    args = parser.parse_args()
    ic_month = int(args.init_month)
    lead_month = int(args.lead_month)
    year = int(args.year)
    model = args.model
    IC_YEAR = int(args.year)
    CONFIGFILE = args.config_file
    OUT_PATH =  os.getcwd() + '/plots/' + calendar.month_abbr[ic_month].upper() + '01/NMME/'

    with open(CONFIGFILE, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    NENS = cfg['EXP']['ensemble_sizes'][0]
    clim_syr = int(cfg["BCSD"]["clim_start_year"])
    clim_eyr = int(cfg["BCSD"]["clim_end_year"])
    sys.path.append(cfg['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.s2splots import plot_utils

    NMME_PATH = 'bcsd_fcst/NMME/'
    if cfg['SETUP']['DATATYPE'] == 'hindcast':
        NMME_PATH = 'hindcast/bcsd_fcst/NMME/'

    mmm = calendar.month_abbr[ic_month].lower() + '01'
    fcast_month = ic_month + lead_month
    if fcast_month > 12:
        fcast_month -= 12
        year = year+1

    lis_month = fcast_month + 1
    lis_year = year
    if lis_month > 12:
        lis_month -= 12
        lis_year = lis_year+1
        clim_syr = clim_syr + 1
        clim_eyr = clim_eyr + 1

    fm_label = str(year) + '-' + calendar.month_abbr[fcast_month] + '_NCEP'
    print (fm_label)

    # python check_nmme_bcsd.py -m CCM4 -l 4 -y 2000
    # Jan 2001
    # compute MODEL mean : 1) raw, 2) bcsd
    nmme_clim_file_raw = NMME_CLIM_TEMPLATE_RAW.format(NMME_PATH,mmm,model)
    print(nmme_clim_file_raw)

    nmme_clim_xr = xr.open_dataset(nmme_clim_file_raw)
    nmme_clim_prcp = np.mean(nmme_clim_xr['clim'].values[lead_month+1,:], axis=0)
    nmme_clim_xr.close()

    nmme_clim_files_bcsd = NMME_CLIM_TEMPLATE_BCSD.format(mmm,model,mmm)
    print(nmme_clim_files_bcsd)

    nmme_clim_bcsd_xr = \
        xr.open_mfdataset(nmme_clim_files_bcsd,concat_dim = 'time', combine='nested')
    nmme_clim_bcsd = \
        nmme_clim_bcsd_xr['PRECTOT'].isel(Lead=lead_month).mean(dim=['time', 'Ens']).values
    nmme_clim_bcsd_xr.close()

    # compute NAFPA mean
    nafpa_clim_files = \
        [NAFPA_FILE_TEMPLATE.format(NAFPA_PATH,cyear, lis_month,cyear,lis_month)
         for cyear in range(clim_syr, clim_eyr + 1)]
    print(nafpa_clim_files)

    nafpa_clim_xr = xr.open_mfdataset(nafpa_clim_files,concat_dim = 'time', combine='nested')
    nafpa_clim = nafpa_clim_xr.mean(dim = 'time')
    nafpa_clim_prcp = nafpa_clim['TotalPrecip_acc'].values

    nafpa_fcst_file = \
        NAFPA_FILE_TEMPLATE.format(NAFPA_PATH,lis_year,lis_month,lis_year,lis_month)
    print(nafpa_fcst_file)

    nafpa_mon_xr = xr.open_dataset(nafpa_fcst_file)
    nafpa_mon_prcp = nafpa_mon_xr['TotalPrecip_acc'].values # [1,2]
    nafpa_anom_prcp = nafpa_mon_prcp - nafpa_clim_prcp  # [2,2]

    # ensemble anomaly
    nmme_monthly_file_bcsd = \
        NMME_MONTHLY_TEMPLATE_BCSD.format(NMME_PATH,mmm,model,mmm,IC_YEAR,IC_YEAR)
    print(nmme_monthly_file_bcsd)

    nmme_mon_bcsd_xr = xr.open_dataset(nmme_monthly_file_bcsd)

    # NOAA-NMME anomaly
    ncep_model = NOAA_NCEP.get(model)
    icdate = '{:04d}/{:02d}/01'.format(IC_YEAR, ic_month)
    if model == 'GNEMO5':
        ddays = (dt.strptime(icdate, "%Y/%m/%d") - dt.strptime('2021/12/01', "%Y/%m/%d")).days
        if ddays < 0:
            ncep_model = 'GEM_NEMO'

    if model == 'GFDL':
        ddays = (dt.strptime(icdate, "%Y/%m/%d") - dt.strptime('2021/02/01', "%Y/%m/%d")).days
        if ddays < 0:
            ncep_model = 'GFDL_FLOR'

    noaa_anom_file = \
        NOAA_NMME_TEMPLATE.format(NOAA_NMME_PATH, ncep_model, ncep_model, IC_YEAR, ic_month)
    print(noaa_anom_file)
    noaa_anom_xr = nc4(noaa_anom_file)
    lati = np.array(noaa_anom_xr.variables['lat'][:])
    loni = np.array(noaa_anom_xr.variables['lon'][:])
    ncep_anom = np.array(noaa_anom_xr.variables['fcst'][:])
    print (ncep_anom[0:NENS.get(model),].shape)
    # regrid 1-deg to 1/4
    ds_in = xr.Dataset(
        {
            "lat": (["lat"], lati),
            "lon": (["lon"], loni),
        })

    ds_in["slice"] = xr.DataArray(
        data = ncep_anom[0,0,:,:],
        dims=["lat", "lon"],
        coords=dict(
            lat=(["lat"], lati),
            lon=(["lon"], loni))
        )

    ds_in["fcst"] = xr.DataArray(
        data = ncep_anom[0:NENS.get(model),],
        dims=["ens","mon", "lat", "lon"],
        coords=dict(
            ens=(["ens"], np.arange(NENS.get(model))),
            mon =(["mon"], np.arange(9)),
            lat=(["lat"], lati),
            lon=(["lon"], loni))
        )

    ds_out = xr.Dataset(
        {
            "lat": (["lat"], nmme_mon_bcsd_xr['latitude'].values),
            "lon": (["lon"], nmme_mon_bcsd_xr['longitude'].values),
        })

    regridder = xe.Regridder(ds_in, ds_out, "bilinear", periodic=True)
    ds_out = regridder(ds_in)

    unit_conv = calendar.monthrange(year,fcast_month)[1]*86400./25.4
    nrows = 2
    ncols = 2
    domain = plot_utils.dicts('boundary', 'GLOBAL')
    load_table = '14WPR'
    under_over = plot_utils.dicts('lowhigh', load_table)
    levels = plot_utils.dicts('anom_levels','Precip_AF')
    plot_title = ['NMME',  'NOAA-NMME', 'BCSD', 'AF10km']
    var_name = 'Precip_AF'
    stitle = var_name + ' Forecast'
    clabel = 'Anomaly (' + plot_utils.dicts('units', var_name) + ')'
    clabel2 = 'Monthly Precip (units mm/d)'

    for ens in range (1, NENS.get(model) +1):
        plot_arr = np.zeros([4,720,1440],dtype=float)
        eee = 'ens' + str(ens)
        nmme_monthly_file_raw = \
            NMME_MONTHLY_TEMPLATE_RAW.format(NMME_PATH,mmm,model,IC_YEAR,eee,mmm,year,fcast_month)
        print(nmme_monthly_file_raw)

        # NMME Raw
        nmme_mon_xr = xr.open_dataset(nmme_monthly_file_raw)
        nmme_mon_prcp = nmme_mon_xr['PRECTOT'].values
        nmme_anom_prcp = nmme_mon_prcp - nmme_clim_prcp

        plot_arr[0,] = nmme_anom_prcp*unit_conv
        plot_arr[1,] = ds_out['fcst'].values[ens-1, lead_month,:,:]*unit_conv
        plot_arr[3,] = nafpa_anom_prcp/25.4

        # NMME BCSD
        bcsd_mon_prcp = \
            nmme_mon_bcsd_xr['PRECTOT'].isel(Lead=lead_month, Ens=ens-1,  time= 0).values
        bcsd_anom_prcp = bcsd_mon_prcp - nmme_clim_bcsd
        plot_arr[2,] = bcsd_anom_prcp*unit_conv

        plot_dir = OUT_PATH + model + '/'
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        # plot anom
        anom_file = plot_dir + \
            '{:04}-{}'.format(IC_YEAR,mmm) + '_' + fm_label + '_' + eee + '_anom.png'
        print (anom_file)
        plot_utils.contours (nmme_mon_xr['lon'].values, nmme_mon_xr['lat'].values, nrows,
                             ncols, plot_arr, load_table, plot_title, domain, anom_file, under_over,
                             fscale=1.1, stitle=stitle, clabel=clabel, levels=levels)

        nmme_mon_xr.close()

    # plot NOAA_NCEP anom on the native grid
    plot_title = ['ENS1',  'ENS2', 'ENS3', 'ENS4']
    loni2 = np.roll(np.where (loni <=180., loni, loni -360.),179)
    ncep_anom2 = np.roll(ncep_anom, 179, axis =3)

    plot_arr = np.zeros([4,181,360],dtype=float)
    plot_arr[0,] = ncep_anom2[0,lead_month,:,:]*unit_conv
    plot_arr[1,] = ncep_anom2[1,lead_month,:,:]*unit_conv
    plot_arr[2,] = ncep_anom2[2,lead_month,:,:]*unit_conv
    plot_arr[3,] = ncep_anom2[3,lead_month,:,:]*unit_conv

    anom_file = plot_dir + '{:04}-{}'.format(IC_YEAR,mmm) + '_' + fm_label + '_RAW_anom.png'
    plot_utils.contours (loni2, lati, nrows,
                         ncols, plot_arr, load_table, plot_title, domain, anom_file, under_over,
                         fscale=1.1, stitle=stitle, clabel=clabel, levels=levels)
