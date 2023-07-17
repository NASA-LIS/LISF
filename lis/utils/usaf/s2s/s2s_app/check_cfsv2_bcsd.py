#!/usr/bin/env python

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

'''
    This is good to plot CFSv2 forcing standerdized anomalies
'''
import os
import sys
import calendar
import argparse
import xarray as xr
import numpy as np
import yaml
#pylint: disable=import-error

NAFPA_PATH = '/discover/nobackup/projects/ghilis/S2S/GLOBAL/DA_Run_Hist/output/SURFACEMODEL/'
NAFPA_FILE_TEMPLATE = '{}/{:04d}{:02d}/LIS_HIST_{:04d}{:02d}010000.d01.nc'

CFSv2_CLIM_TEMPLATE_RAW = '{}/raw/Climatology/{}/{}_fcst_clim.nc'
CFSv2_MONTHLY_TEMPLATE_RAW = '{}/raw/Monthly/{}/{:04d}/{}/{}.cfsv2.{:04d}{:02d}.nc'

CFSv2_CLIM_TEMPLATE_BCSD = 'hindcast/bcsd_fcst/CFSv2_25km/bcsd/Monthly/{}/{}.CFSv2.{}_*_*.nc'
CFSv2_MONTHLY_TEMPLATE_BCSD = '{}/bcsd/Monthly/{}/{}.CFSv2.{}_{:04d}_{:04d}.nc'

lis_name = {
        'PRECTOT': 'TotalPrecip_acc','PS': 'Psurf_f_tavg','T2M': 'Tair_f_tavg',
        'LWS': 'LWdown_f_tavg','SLRSF': 'SWdown_f_tavg','Q2M': 'Qair_f_tavg',
        'WIND10M': 'Wind_f_tavg',}

mean_lwval = {
    'SWdown_f_tavg': 0.,
    'LWdown_f_tavg': 50.,
    'Psurf_f_tavg': 50000.,
    'Qair_f_tavg': 0.,
    'Tair_f_tavg': 200.,
    'Wind_f_tavg': 0.3,
    'TotalPrecip_acc': 0
    }

mean_upval = {
    'SWdown_f_tavg': 450,
    'LWdown_f_tavg': 450,
    'Psurf_f_tavg': 103000,
    'Qair_f_tavg': 2500,
    'Tair_f_tavg': 320.,
    'Wind_f_tavg': 20,
    'TotalPrecip_acc': 50
    }

color_table = {
        'T2M': '14WT2M'
        }
unit_conv = {
       'T2M': 9./5.
        }

NENS = 12

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--init_month', required=True, help='init month 1-12')
    parser.add_argument('-l', '--lead_month', required=True, help='lead month 0-9')
    parser.add_argument('-y', '--year', required=True, help='forcing start year')
    parser.add_argument('-c', '--config_file', required=True, help='config file')
    parser.add_argument('-v', '--variable', required=True,
                        help='variable [LWS PRECTOT PS Q2M SLRSF T2M WIND10M]')

    args = parser.parse_args()
    ic_month = int(args.init_month)
    lead_month = int(args.lead_month)
    year = int(args.year)
    variable = args.variable
    IC_YEAR = int(args.year)
    CONFIGFILE = args.config_file

    OUT_PATH =  os.getcwd() + '/plots/' + calendar.month_abbr[ic_month].upper() + '01/CFSv2/'

    with open(CONFIGFILE, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)

    sys.path.append(cfg['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.s2splots import plot_utils
    clim_syr = int(cfg["BCSD"]["clim_start_year"])
    clim_eyr = int(cfg["BCSD"]["clim_end_year"])
    CFSv2_PATH = 'bcsd_fcst/CFSv2_25km/'
    if cfg['SETUP']['DATATYPE'] == 'hindcast':
        CFSv2_PATH = 'hindcast/bcsd_fcst/CFSv2_25km/'

    var_levels = np.linspace(mean_lwval.get(lis_name.get(variable)),
                             mean_upval.get(lis_name.get(variable)), 21)
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

    fm_label = str(year) + '-' + calendar.month_abbr[fcast_month]
    print (fm_label)

    # python check_nmme_bcsd.py -m CCM4 -l 4 -y 2000
    # Jan 2001
    # compute VARIABLE mean : 1) raw, 2) bcsd
    varname_clim_file_raw = CFSv2_CLIM_TEMPLATE_RAW.format(CFSv2_PATH,mmm,variable)
    print(varname_clim_file_raw)

    varname_clim_xr = xr.open_dataset(varname_clim_file_raw)
    varname_clim_var = np.mean(varname_clim_xr['clim'].values[lead_month+1,:], axis=0)
    varname_clim_std = np.std(varname_clim_xr['clim'].values[lead_month+1,:], axis=0)
    varname_clim_xr.close()

    varname_clim_files_bcsd = CFSv2_CLIM_TEMPLATE_BCSD.format(mmm,variable,mmm)
    print(varname_clim_files_bcsd)

    varname_clim_bcsd_xr = \
        xr.open_mfdataset(varname_clim_files_bcsd,concat_dim = 'time', combine='nested')
    varname_clim_bcsd = \
        varname_clim_bcsd_xr[variable].isel(Lead=lead_month).mean(dim=['time', 'Ens']).values
    varname_clim_bcsd_std = \
        varname_clim_bcsd_xr[variable].isel(Lead=lead_month).std(dim=['time', 'Ens']).values
    varname_clim_bcsd_xr.close()

    # compute NAFPA mean
    nafpa_clim_files = \
        [NAFPA_FILE_TEMPLATE.format(NAFPA_PATH,cyear, lis_month,cyear,lis_month)
         for cyear in range(clim_syr, clim_eyr + 1)]
    print(nafpa_clim_files)

    nafpa_clim_xr = xr.open_mfdataset(nafpa_clim_files,concat_dim = 'time', combine='nested')
    nafpa_clim = nafpa_clim_xr.mean(dim = 'time')
    nafpa_std = nafpa_clim_xr.std(dim = 'time')
    nafpa_clim_var = nafpa_clim[lis_name.get(variable)].values
    nafpa_std_var = nafpa_std[lis_name.get(variable)].values

    nafpa_fcst_file = NAFPA_FILE_TEMPLATE.format(NAFPA_PATH,lis_year,lis_month,lis_year,lis_month)
    print(nafpa_fcst_file)

    nafpa_mon_xr = xr.open_dataset(nafpa_fcst_file)
    nafpa_mon_var = nafpa_mon_xr[lis_name.get(variable)].values # [1,2]
    nafpa_anom_var = (nafpa_mon_var - nafpa_clim_var)/ nafpa_std_var  # [2,2]

    # ensemble anomaly
    varname_monthly_file_bcsd = \
        CFSv2_MONTHLY_TEMPLATE_BCSD.format(CFSv2_PATH,mmm,variable,mmm,IC_YEAR,IC_YEAR)
    print(varname_monthly_file_bcsd)

    varname_mon_bcsd_xr = xr.open_dataset(varname_monthly_file_bcsd)
    #unit_conv = calendar.monthrange(year,fcast_month)[1]*86400./25.4
    domain = plot_utils.dicts('boundary', 'GLOBAL')
    load_table = 'CB11W_'
    under_over = plot_utils.dicts('lowhigh', load_table)
    levels = levels = plot_utils.dicts('anom_levels', 'standardized')
    plot_title = ['CFSv2', 'BCSD', 'AF10km']
    stitle = variable + ' Forecast'
    clabel = 'Standardized Anomaly'
    clabel2 = 'Monthly ' + variable

    for ens in range (1, NENS +1):
        plot_arr = np.zeros([3,720,1440],dtype=float)
        eee = 'ens' + str(ens)
        varname_monthly_file_raw = \
            CFSv2_MONTHLY_TEMPLATE_RAW.format(CFSv2_PATH,mmm,IC_YEAR,eee,mmm,year,fcast_month)
        print(varname_monthly_file_raw)

        # CFSv2 Raw
        varname_mon_xr = xr.open_dataset(varname_monthly_file_raw)
        varname_mon_var = varname_mon_xr[variable].values    # [1,1]
        varname_anom_var = (varname_mon_var - varname_clim_var)/varname_clim_std

        plot_arr[0,] = varname_anom_var
        plot_arr[2,] = nafpa_anom_var

        # CFSv2 BCSD
        bcsd_mon_var = \
            varname_mon_bcsd_xr[variable].isel(Lead=lead_month, Ens=ens-1,  time= 0).values
        bcsd_anom_var = (bcsd_mon_var - varname_clim_bcsd)/varname_clim_bcsd_std
        plot_arr[1,] = bcsd_anom_var

        plot_dir = OUT_PATH + variable + '/'
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        # plot anom
        anom_file = plot_dir + \
            '{:04d}-{}'.format(IC_YEAR,mmm) + '_' + fm_label + '_' + eee + '_anom.png'
        print (anom_file)
        plot_utils.contours (varname_mon_xr['lon'].values, varname_mon_xr['lat'].values, 3,
                             1, plot_arr, load_table, plot_title, domain, anom_file, under_over,
                             fscale=1.1, stitle=stitle, clabel=clabel, levels=levels)
        # plot monthly
        mon_file = plot_dir + '{:04d}-{}'.format(IC_YEAR,mmm) + '_' + fm_label + '_' + eee + '.png'
        plot_arr[0,] = varname_mon_var
        plot_arr[2,] = nafpa_mon_var
        plot_arr[1,] = bcsd_mon_var
        plot_utils.contours (varname_mon_xr['lon'].values, varname_mon_xr['lat'].values, 3,
                             1, plot_arr, 'L21', plot_title, domain, mon_file, under_over,
                             fscale=1.1, stitle=stitle, clabel=clabel2, levels=var_levels)
        varname_mon_xr.close()
