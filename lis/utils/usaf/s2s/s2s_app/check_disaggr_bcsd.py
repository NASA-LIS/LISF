import os
import sys
import numpy as np
import calendar
import xarray as xr
import yaml
import argparse
import yaml

# NMME
# /gpfsm/dnb06/projects/p204/S2S/GLOBAL/E2ES_557ww-7.5/hindcast/bcsd_fcst/NMME/bcsd/Monthly/sep01/PRECTOT.GNEMO5.sep01_1991_1991.nc
# float PRECTOT(time, Lead, Ens, latitude, longitude) ;
#                PRECTOT:_FillValue = 9.96921e+36f ;
#                PRECTOT:units = "kg/m^2/s" ;
# /gpfsm/dnb06/projects/p204/S2S/GLOBAL/E2ES_557ww-7.5/hindcast/bcsd_fcst/NMME/bcsd/6-Hourly/sep01/GNEMO5/1991/ens[1-10]/PRECTOT.199109.nc4
# PRECTOT.199109.nc4
# PRECTOT.199110.nc4
# PRECTOT.199111.nc4
# PRECTOT.199112.nc4
# PRECTOT.199201.nc4
# PRECTOT.199202.nc4
# PRECTOT.199203.nc4
# PRECTOT.199204.nc4
# PRECTOT.199205.nc4
#        float PRECTOT(time, latitude, longitude) ;
#                PRECTOT:_FillValue = -9999.f ;
#                PRECTOT:least_significant_digit = 5 ;
#                PRECTOT:units = "kg/m^2/s" ;
#                PRECTOT:standard_name = "PRECTOT" ;
#        time = UNLIMITED ; // (120 currently)
#        longitude = 1440 ;
#        latitude = 720 ;

# CFSv2
# /gpfsm/dnb06/projects/p204/S2S/GLOBAL/E2ES_557ww-7.5/hindcast/bcsd_fcst/CFSv2_25km/bcsd/Monthly/sep01/T2M.CFSv2.sep01_1991_1991.nc
#         longitude = 1440 ;
#        latitude = 720 ;
#        Lead = 9 ;
#        Ens = 12 ;
#        time = UNLIMITED ; // (1 currently)
#        float T2M(time, Lead, Ens, latitude, longitude) ;
#                T2M:_FillValue = 9.96921e+36f ;
#                T2M:units = "K" ;

# /gpfsm/dnb06/projects/p204/S2S/GLOBAL/E2ES_557ww-7.5/hindcast/bcsd_fcst/CFSv2_25km/bcsd/6-Hourly/sep01/1991/ens1/T2M.199109.nc4
#  float T2M(time, lat, lon) ;
#                T2M:_FillValue = -9999.f ;
#                T2M:least_significant_digit = 5 ;
#                T2M:units = "K" ;
#                T2M:standard_name = "T2M" ;

lis_name = {
        'PRECTOT': 'Rainf_f_tavg','PS': 'Psurf_f_tavg','T2M': 'Tair_f_tavg',
        'LWS': 'LWdown_f_tavg','SLRSF': 'SWdown_f_tavg','Q2M': 'Qair_f_tavg',
        'WIND10M': 'Wind_f_tavg',}

mean_lwval = {
    'SWdown_f_tavg': 0.,
    'LWdown_f_tavg': 50.,
    'Psurf_f_tavg': 50000.,
    'Qair_f_tavg': 0.,
    'Tair_f_tavg': 200.,
    'Wind_f_tavg': 0.3,
    'Rainf_f_tavg': 0
    }

mean_upval = {
    'SWdown_f_tavg': 450,
    'LWdown_f_tavg': 450,
    'Psurf_f_tavg': 103000,
    'Qair_f_tavg': 2500,
    'Tair_f_tavg': 320.,
    'Wind_f_tavg': 20,
    'Rainf_f_tavg': 50
    }

if __name__ == "__main__":
    
    """Main driver."""
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
    OUT_PATH =  os.getcwd() + '/plots/DISAGG/' + forcing + '/'
    with open(CONFIGFILE, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)

    sys.path.append(cfg['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.s2splots import plot_utils
    
    under_over = ['black', '#B404AE']
    plot_title = ['Monthly', 'Monthly mean from 6-hourly']
    domain = plot_utils.dicts('boundary', 'GLOBAL')
    
    if forcing == 'NMME':
        rainf_levels = [0,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,7,9,11,13,15,17,20,25,30,50]        
        clabel2 = 'Monthly Precip (units mm/d)'
        stitle = 'PRECTOT'
        NENS = cfg['EXP']['ensemble_sizes'][0]
        for nmme_model in  cfg['EXP']['NMME_models']:
            plot_dir = OUT_PATH + nmme_model + '/'
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
                
            # Read Monthly
            mon_file = 'bcsd_fcst/NMME/bcsd/Monthly/' + '{}/PRECTOT.{}.{}_{:04d}_{:04d}.nc'.format(mmm, nmme_model,mmm, IC_YEAR, IC_YEAR)
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
                    eee = 'ens' + str(ens)
                    six_file = 'bcsd_fcst/NMME/bcsd/6-Hourly/' + '{}/{}/{:04d}/{}/PRECTOT.{:04d}{:02d}.nc4'.format(mmm,nmme_model,IC_YEAR,eee,year,fcast_month)
                    print ('Comparing ' + nmme_model, ': L=',lead_month, ' E=',ens)
                    six_xr = xr.open_dataset(six_file)
                    six_xarr = six_xr['PRECTOT'] # .to_array.mean(dim=['time']).to_numpy()
                    six_arr = six_xarr.mean(dim=['time'])

                    figname = plot_dir + 'PRECTOT.{:04d}_{}_{}.{:04d}{:02d}.png'.format(IC_YEAR, mmm, eee, year,fcast_month)
                    plot_arr[0,] = mon_arr.to_numpy()*86400.
                    plot_arr[1,] = six_arr.to_numpy()*86400.
                    plot_utils.contours (monthly_xr['longitude'].values, monthly_xr['latitude'].values, 2,
                                         1, plot_arr, 'L21', plot_title, domain, figname, under_over,
                                         fscale=1.1, stitle=stitle, clabel=clabel2, levels=rainf_levels)
                    six_xr.close()            
            monthly_xr.close()
        
    else:
        variable = 'T2M'
        clabel2 = 'Monthly ' + variable + '(K)'
        var_levels = np.linspace(mean_lwval.get(lis_name.get(variable)), mean_upval.get(lis_name.get(variable)), 21)
        stitle = variable
        NENS = 12
        
        # Read Monthly
        mon_file = 'bcsd_fcst/CFSv2_25km/bcsd/Monthly/' + '{}/{}.CFSv2.{}_{:04d}_{:04d}.nc'.format(mmm, variable,mmm, IC_YEAR, IC_YEAR)
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
                mon_arr = monthly_xr[variable].isel(Lead=lead_month, Ens = ens -1, time =0)
                eee = 'ens' + str(ens)
                six_file = 'bcsd_fcst/CFSv2_25km/bcsd/6-Hourly/' + '{}/{:04d}/{}/{}.{:04d}{:02d}.nc4'.format(mmm,IC_YEAR,eee,variable,year,fcast_month)
                print ('Comparing CFASv2 ' + variable + ': L=',lead_month, ' E=',ens)
                six_xr = xr.open_dataset(six_file)
                six_xarr = six_xr[variable] # .to_array.mean(dim=['time']).to_numpy()
                six_arr = six_xarr.mean(dim=['time'])
                
                figname = plot_dir + variable + '.{:04d}_{}_{}.{:04d}{:02d}.png'.format(IC_YEAR, mmm, eee, year,fcast_month)
                plot_arr[0,] = mon_arr.to_numpy()
                plot_arr[1,] = six_arr.to_numpy()
                plot_utils.contours (monthly_xr['longitude'].values, monthly_xr['latitude'].values, 2,
                                    1, plot_arr, 'L21', plot_title, domain, figname, under_over,
                                     fscale=1.1, stitle=stitle, clabel=clabel2, levels=var_levels)
                six_xr.close()            
        monthly_xr.close()
