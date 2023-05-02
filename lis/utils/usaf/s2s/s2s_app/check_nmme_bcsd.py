#!/usr/bin/env python

import os
import sys
import numpy as np
import calendar
import xarray as xr
import yaml
import argparse
import yaml

# Observed clim: ryan/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/PRECTOT_obs_clim.nc (DIST, time, latitude, longitude) - DIST 1-13 Jan-Dec
#  /discover/nobackup/projects/usaf_lis/razamora/NAFPAMod_Aggregate_Monthly/2005/LIS_HIST_200503010000.d01.nc
# 6 panel plots per ens, per model, per month
# NMME monthly        NMME anom
# AF10KM monthly      AF10KM anom
# BCSD monthly        BCSD anom


# NMME
# /discover/nobackup/projects/usaf_lis/razamora/NAFPAMod_Aggregate_Monthly/1991-2021/LIS_HIST_{199111}010000.d01.nc
# BCSD
#

NAFPA_PATH = '/discover/nobackup/projects/ghilis/S2S/GLOBAL/Forcing_Merge/M2CH2BC_USAFNAFPAMod_S2S_Mon/SURFACEMODEL/'
NMME_PATH = os.getcwd() + '/bcsd_fcst/NMME/'
OUT_PATH =  os.getcwd() + '/plots/NMME/'

NAFPA_CLIM_TEMPLATE = '{}/*{:02d}/LIS_HIST_*{:02d}010000.d01.nc'
NAFPA_FILE_TEMPLATE = '{}/{:04d}{:02d}/LIS_HIST_{:04d}{:02d}010000.d01.nc'

NMME_CLIM_TEMPLATE_RAW = '{}/raw/Climatology/{}/{}/PRECTOT_fcst_clim.nc'
NMME_MONTHLY_TEMPLATE_RAW = '{}/raw/Monthly/{}/{}/{:04d}/{}/{}.nmme.monthly.{:04}{:02d}.nc'

NMME_CLIM_TEMPLATE_BCSD = '{}/bcsd/Monthly/{}/PRECTOT.{}.{}_*_*.nc'
NMME_MONTHLY_TEMPLATE_BCSD = '{}/bcsd/Monthly/{}/PRECTOT.{}.{}_{:04d}_{:04d}.nc'

rainf_levels = [0,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,7,9,11,13,15,17,20,25,30,50]

if __name__ == "__main__":
    
    """Main driver."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--init_month', required=True, help='init month 1-12')
    parser.add_argument('-l', '--lead_month', required=True, help='lead month 0-9')
    parser.add_argument('-y', '--year', required=True, help='forcing start year')
    parser.add_argument('-c', '--config_file', required=True, help='config file')    
    parser.add_argument('-m', '--model', required=True, help='NMME model [CCM4, CCSM4, CFSv2, GEOSv2, GFDL, GNEMO5]')
 
    args = parser.parse_args()    
    ic_month = int(args.init_month)
    lead_month = int(args.lead_month)
    year = int(args.year)
    model = args.model
    IC_YEAR = int(args.year)
    CONFIGFILE = args.config_file
    
    with open(CONFIGFILE, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    NENS = cfg['EXP']['ensemble_sizes'][0]
    sys.path.append(cfg['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.s2splots import plot_utils

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
    fm_label = str(year) + '-' + calendar.month_abbr[fcast_month]
    print (fm_label)

    # python check_nmme_bcsd.py -m CCM4 -l 4 -y 2000
    # Jan 2001
    # compute MODEL mean : 1) raw, 2) bcsd
    nmme_clim_file_raw = NMME_CLIM_TEMPLATE_RAW.format(NMME_PATH,mmm,model)
    print(nmme_clim_file_raw)
    # /discover/nobackup/projects/ghilis/S2S/GLOBAL/E2ES_557ww-7.5/hindcast/bcsd_fcst/NMME//raw/Climatology/sep01/CCM4/PRECTOT_fcst_clim.nc
    nmme_clim_xr = xr.open_dataset(nmme_clim_file_raw)
    nmme_clim_prcp = np.mean(nmme_clim_xr['clim'].values[lead_month+1,:], axis=0)
    nmme_clim_xr.close()

    nmme_clim_files_bcsd = NMME_CLIM_TEMPLATE_BCSD.format(NMME_PATH,mmm,model,mmm)
    print(nmme_clim_files_bcsd)
    # /discover/nobackup/projects/ghilis/S2S/GLOBAL/E2ES_557ww-7.5/hindcast/bcsd_fcst/NMME//bcsd/Monthly/sep01/PRECTOT.CCM4.sep01_*_*.nc
    nmme_clim_bcsd_xr = xr.open_mfdataset(nmme_clim_files_bcsd,concat_dim = 'time', combine='nested')
    nmme_clim_bcsd = nmme_clim_bcsd_xr['PRECTOT'].isel(Lead=lead_month).mean(dim=['time', 'Ens']).values
    nmme_clim_bcsd_xr.close()
        
    # compute NAFPA mean
    nafpa_clim_files = NAFPA_CLIM_TEMPLATE.format(NAFPA_PATH,lis_month,lis_month)
    print(nafpa_clim_files)
    # /discover/nobackup/projects/ghilis/S2S/GLOBAL/Forcing_Merge/M2CH2BC_USAFNAFPAMod_S2S_Mon/SURFACEMODEL//*02/LIS_HIST_*02010000.d01.nc
    nafpa_clim_xr = xr.open_mfdataset(nafpa_clim_files,concat_dim = 'time', combine='nested')
    nafpa_clim = nafpa_clim_xr.mean(dim = 'time')
    nafpa_clim_prcp = nafpa_clim['Rainf_f_tavg'].values
    
    nafpa_fcst_file = NAFPA_FILE_TEMPLATE.format(NAFPA_PATH,lis_year,lis_month,lis_year,lis_month) 
    print(nafpa_fcst_file)
    # /discover/nobackup/projects/ghilis/S2S/GLOBAL/Forcing_Merge/M2CH2BC_USAFNAFPAMod_S2S_Mon/SURFACEMODEL//200102/LIS_HIST_200102010000.d01.nc
    nafpa_mon_xr = xr.open_dataset(nafpa_fcst_file)
    nafpa_mon_prcp = nafpa_mon_xr['Rainf_f_tavg'].values # [1,2]
    nafpa_anom_prcp = nafpa_mon_prcp - nafpa_clim_prcp   # [2,2]

    # ensemble anomaly
    nmme_monthly_file_bcsd = NMME_MONTHLY_TEMPLATE_BCSD.format(NMME_PATH,mmm,model,mmm,IC_YEAR,IC_YEAR)
    print(nmme_monthly_file_bcsd)
    # /discover/nobackup/projects/ghilis/S2S/GLOBAL/E2ES_557ww-7.5/hindcast/bcsd_fcst/NMME//bcsd/Monthly/sep01/PRECTOT.CCM4.sep01_2000_2000.nc
    nmme_mon_bcsd_xr = xr.open_dataset(nmme_monthly_file_bcsd)

    unit_conv = calendar.monthrange(year,fcast_month)[1]*86400./25.4
    nrows = 3
    ncols = 1
    domain = plot_utils.dicts('boundary', 'GLOBAL')
    load_table = '14WPR'
    levels = plot_utils.dicts('anom_levels','Precip_AF')
    plot_title = ['NMME', 'BCSD', 'AF10km']
    under_over = ['black', '#B404AE']
    var_name = 'Precip_AF'
    stitle = var_name + ' Forecast'
    clabel = 'Anomaly (' + plot_utils.dicts('units', var_name) + ')'
    clabel2 = 'Monthly Precip (units mm/d)'
    
    for ens in range (1, NENS.get(model) +1):
        plot_arr = np.zeros([3,720,1440],dtype=float)
        eee = 'ens' + str(ens)
        nmme_monthly_file_raw = NMME_MONTHLY_TEMPLATE_RAW.format(NMME_PATH,mmm,model,IC_YEAR,eee,mmm,year,fcast_month)
        print(nmme_monthly_file_raw)
        # /discover/nobackup/projects/ghilis/S2S/GLOBAL/E2ES_557ww-7.5/hindcast/bcsd_fcst/NMME//raw/Monthly/sep01/CCM4/2000/ens1/sep01.nmme.monthly.200101.nc

        # NMME Raw
        nmme_mon_xr = xr.open_dataset(nmme_monthly_file_raw)
        nmme_mon_prcp = nmme_mon_xr['PRECTOT'].values    # [1,1]
        nmme_anom_prcp = nmme_mon_prcp - nmme_clim_prcp  # [1,2]

        plot_arr[0,] = nmme_anom_prcp*unit_conv
        plot_arr[2,] = nafpa_anom_prcp*unit_conv
        
        # NMME BCSD
        bcsd_mon_prcp = nmme_mon_bcsd_xr['PRECTOT'].isel(Lead=lead_month, Ens=ens-1,  time= 0).values # [1,3]
        bcsd_anom_prcp = bcsd_mon_prcp - nmme_clim_bcsd                                               # [2,3]
        plot_arr[1,] = bcsd_anom_prcp*unit_conv
        
        plot_dir = OUT_PATH + model + '/'
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        # plot anom
        anom_file = plot_dir + '{:04}-{}'.format(IC_YEAR,mmm) + '_' + fm_label + '_' + eee + '_anom.png'
        print (anom_file)
        plot_utils.contours (nmme_mon_xr['lon'].values, nmme_mon_xr['lat'].values, nrows,
                             ncols, plot_arr, load_table, plot_title, domain, anom_file, under_over,
                             fscale=1.1, stitle=stitle, clabel=clabel, levels=levels)
        # plot monthly
        mon_file = plot_dir + '{:04}-{}'.format(IC_YEAR,mmm) + '_' + fm_label + '_' + eee + '.png'
        plot_arr[0,] = nmme_mon_prcp*86400.
        plot_arr[2,] = nafpa_mon_prcp*86400.
        plot_arr[1,] = bcsd_mon_prcp*86400.
        plot_utils.contours (nmme_mon_xr['lon'].values, nmme_mon_xr['lat'].values, nrows,
                             ncols, plot_arr, 'L21', plot_title, domain, mon_file, under_over,
                             fscale=1.1, stitle=stitle, clabel=clabel2, levels=rainf_levels)
        
        
        
        nmme_mon_xr.close()
