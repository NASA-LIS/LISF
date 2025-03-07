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
    This is good to check Monthly to 6-hourly disaggregation
'''
import os
import sys
import calendar
import argparse
import xarray as xr
import numpy as np
import yaml
#pylint: disable=import-error

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

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--init_month', required=True, help='init month 1-12')
    parser.add_argument('-y', '--year', required=True, help='forcing start year')
    parser.add_argument('-f', '--forcing', required=True, help='forcing (NMME, CFSv2')
    parser.add_argument('-c', '--config_file', required=True, help='config file')

    args = parser.parse_args()
    ic_month = int(args.init_month)
    year = int(args.year)
    forcing = args.forcing
    CONFIGFILE = args.config_file

    IC_YEAR = int(args.year)
    mmm = calendar.month_abbr[ic_month].lower() + '01'
    OUT_PATH =  os.getcwd() + '/plots/' + \
        calendar.month_abbr[ic_month].upper() + '01/DISAGG/' + forcing + '/'
    with open(CONFIGFILE, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)

    sys.path.append(cfg['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.s2splots import plot_utils

    BCSD_PATH =  'bcsd_fcst/'
    if cfg['SETUP']['DATATYPE'] == 'hindcast':
        BCSD_PATH =  'hindcast/bcsd_fcst/'
    LOAD_TABLE = 'L21'
    under_over = plot_utils.dicts('lowhigh', LOAD_TABLE)
    plot_title = ['Monthly', 'Monthly mean from 6-hourly']
    domain = plot_utils.dicts('boundary', 'GLOBAL')

    if forcing == 'NMME':
        rainf_levels = \
            [0,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,7,9,11,13,15,17,20,25,30,50]
        CLABEL2 = 'Monthly Precip (units mm/d)'
        STITLE = 'PRECTOT'
        NENS = cfg['EXP']['ensemble_sizes'][0]
        for nmme_model in  cfg['EXP']['NMME_models']:
            plot_dir = OUT_PATH + nmme_model + '/'
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)

            # Read Monthly
            mon_file = BCSD_PATH + '/NMME/bcsd/Monthly/' +\
                '{}/PRECTOT.{}.{}_{:04d}_{:04d}.nc'.format(mmm, nmme_model,mmm, IC_YEAR, IC_YEAR)
            monthly_xr = xr.open_dataset(mon_file)
            print ('NMME Monthly File : ', mon_file)

            # Read 6-hourly
            for lead_month in range(0,9):
                year = IC_YEAR
                fcast_month = ic_month + lead_month
                if fcast_month > 12:
                    fcast_month -= 12
                    year = year+1
                for ens in range (1, NENS.get(nmme_model) +1):
                    plot_arr = np.zeros([2,720,1440],dtype=float)
                    mon_arr = monthly_xr['PRECTOT'].isel(Lead=lead_month, Ens = ens -1, time =0)
                    EEE = 'ens' + str(ens)
                    six_file = BCSD_PATH + '/NMME/bcsd/6-Hourly/' + '{}/{}/{:04d}/{}/'\
                        'PRECTOT.{:04d}{:02d}.nc4'.format(
                            mmm,nmme_model,IC_YEAR,EEE,year,fcast_month)
                    print ('Comparing ' + nmme_model, ': L=',lead_month, ' E=',ens)
                    six_xr = xr.open_dataset(six_file)
                    six_xarr = six_xr['PRECTOT'] # .to_array.mean(dim=['time']).to_numpy()
                    six_arr = six_xarr.mean(dim=['time'])

                    figname = plot_dir + \
                        'PRECTOT.{:04d}_{}_{}.{:04d}{:02d}.png'.format(
                            IC_YEAR, mmm, EEE, year,fcast_month)
                    plot_arr[0,] = mon_arr.to_numpy()*86400.
                    plot_arr[1,] = six_arr.to_numpy()*86400.
                    plot_utils.contours (
                        monthly_xr['longitude'].values, monthly_xr['latitude'].values, 2,
                        1, plot_arr, LOAD_TABLE, plot_title, domain, figname, under_over,
                        fscale=1.1, stitle=STITLE, clabel=CLABEL2, levels=rainf_levels)
                    six_xr.close()
            monthly_xr.close()

    else:
        VARIABLE = 'T2M'
        CLABEL2 = 'Monthly ' + VARIABLE + '(K)'
        var_levels = np.linspace(mean_lwval.get(lis_name.get(VARIABLE)),
                                 mean_upval.get(lis_name.get(VARIABLE)), 21)
        STITLE = VARIABLE
        NENS = 12

        # Read Monthly
        mon_file =  BCSD_PATH + '/CFSv2_25km/bcsd/Monthly/' + \
            '{}/{}.CFSv2.{}_{:04d}_{:04d}.nc'.format(mmm, VARIABLE,mmm, IC_YEAR, IC_YEAR)
        monthly_xr = xr.open_dataset(mon_file)
        print ('CFSv2 Monthly File : ', mon_file)
        plot_dir = OUT_PATH
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        # Read 6-hourly
        for lead_month in range(0,9):
            year = IC_YEAR
            fcast_month = ic_month + lead_month
            if fcast_month > 12:
                fcast_month -= 12
                year = year+1
            for ens in range (1, NENS +1):
                plot_arr = np.zeros([2,720,1440],dtype=float)
                mon_arr = monthly_xr[VARIABLE].isel(Lead=lead_month, Ens = ens -1, time =0)
                EEE = 'ens' + str(ens)
                six_file = BCSD_PATH + '/CFSv2_25km/bcsd/6-Hourly/' + \
                    '{}/{:04d}/{}/{}.{:04d}{:02d}.nc4'.format(
                        mmm,IC_YEAR,EEE,VARIABLE,year,fcast_month)
                print ('Comparing CFSv2 ' + VARIABLE + ': L=',lead_month, ' E=',ens)
                six_xr = xr.open_dataset(six_file)
                six_xarr = six_xr[VARIABLE] # .to_array.mean(dim=['time']).to_numpy()
                six_arr = six_xarr.mean(dim=['time'])

                figname = plot_dir + VARIABLE + \
                    '.{:04d}_{}_{}.{:04d}{:02d}.png'.format(IC_YEAR, mmm, EEE, year,fcast_month)
                plot_arr[0,] = mon_arr.to_numpy()
                plot_arr[1,] = six_arr.to_numpy()
                plot_utils.contours (
                    monthly_xr['longitude'].values, monthly_xr['latitude'].values, 2,
                    1, plot_arr, 'L21', plot_title, domain, figname, under_over,
                    fscale=1.1, stitle=STITLE, clabel=CLABEL2, levels=var_levels)
                six_xr.close()
        monthly_xr.close()
