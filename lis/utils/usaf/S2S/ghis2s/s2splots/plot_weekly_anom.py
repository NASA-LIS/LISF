import os
import sys
import calendar
from datetime import datetime, date
from dateutil.relativedelta import relativedelta
import argparse
import xarray as xr
import numpy as np
import yaml
# pylint: disable=import-error
import plot_utils
from ghis2s.s2smetric.metricslib import get_anom
# pylint: enable=import-error

def plot_anoms(fcst_year, fcst_mon, cwd, config, region, anom_type):
    plotdir_template = cwd + '/s2splots/{:04d}{:02d}/' 
    plotdir = plotdir_template.format(fcst_year, fcst_mon)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir, exist_ok=True)

    figure_template = '{}/NMME_plot_{}_{}_weeks1-6_{}.png'
    
    lead_week = [0, 1, 2, 3, 4, 5]
    # Universal setup of plots:
    nrows = 2
    ncols = 3
    domain = plot_utils.dicts('boundary', region)

    data_dir = cwd + f'/s2smetric/{fcst_year:04d}{fcst_mon:02d}/'
    cartopy_dir = config['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/'
    
    for var_name in config["POST"]["weekly_vars"]:
        if anom_type == 'ANOM':
            clabel = 'Anomaly (' + plot_utils.dicts('units', var_name) + ')'
            load_table = 'clim_reanaly'
        else:
            if var_name == 'RZSM':
                continue
            clabel = 'Standardized Anomaly'            
            if var_name == 'TOP40ST':
                load_table = 'CB11W_'
            else:
                load_table = 'CB11W'

        under_over = plot_utils.dicts('lowhigh', load_table)               
        # READ ANOMALIES
        anom = get_anom(data_dir, var_name, anom_type, weekly=True)
        anom_crop = plot_utils.crop(domain, anom)
        median_anom = np.median(anom_crop.anom.values, axis=0)
        plot_arr = median_anom[lead_week, ]
        if anom_type == 'SANOM':
            levels=plot_utils.dicts('anom_levels', 'standardized')
            plot_arr = np.where(plot_arr < np.max(levels), plot_arr, np.nan)
            plot_arr = np.where(plot_arr > np.min(levels), plot_arr, np.nan)
            
        BEGDATE = date(fcst_year, fcst_mon, 2)
        titles = []
        for lead in lead_week:
            ENDDATE = BEGDATE + relativedelta(days=6)
            titles.append(var_name + ' '+  BEGDATE.strftime("%Y%m%d") + '-' + ENDDATE.strftime("%Y%m%d"))
            BEGDATE += relativedelta(days=7)
            
        maxloc = np.unravel_index(np.nanargmax(plot_arr), plot_arr.shape)
        figure = figure_template.format(plotdir, region, var_name, anom_type.lower())
        stitle = var_name + ' Forecast'
        if anom_type == 'ANOM':
            anom_minmax = plot_utils.dicts('anom_minmax', var_name) 
            plot_utils.contours (anom_crop.lon.values, anom_crop.lat.values, nrows,
                                 ncols, plot_arr, load_table, titles, domain,
                                 figure, under_over,
                                 fscale=1.2, stitle=stitle, clabel=clabel, min_val=anom_minmax[0], max_val=anom_minmax[1],
                                 cartopy_datadir=cartopy_dir, projection=['polar', 90.])
        else:
            plot_utils.contours (anom_crop.lon.values, anom_crop.lat.values, nrows,
                                 ncols, plot_arr, load_table, titles, domain,
                                 figure, under_over,
                                 fscale=1.2, stitle=stitle, clabel=clabel, levels=levels,
                                 cartopy_datadir=cartopy_dir, projection=['polar', 90.])
        #sys.exit()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--fcst_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--fcst_mon', required=True, help= 'forecast end year')
    parser.add_argument('-c', '--configfile', required=True, help='config file name')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    
    args = parser.parse_args()
    configfile = args.configfile
    fcst_year = int(args.fcst_year)
    fcst_mon = int(args.fcst_mon)
    cwd = args.cwd
    
    # load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    plot_anoms(fcst_year, fcst_mon, cwd, config, 'ARCTIC', 'ANOM')
    plot_anoms(fcst_year, fcst_mon, cwd, config, 'ARCTIC', 'SANOM')
