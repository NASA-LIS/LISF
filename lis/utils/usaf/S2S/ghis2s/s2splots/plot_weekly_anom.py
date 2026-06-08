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
New s2splots weekly scripts (S. Mahanama; Jul-2025)
'''
# pylint: disable=no-value-for-parameter
# pylint: disable=f-string-without-interpolation,too-many-positional-arguments
# pylint: disable=too-many-arguments,too-many-locals,consider-using-f-string,too-many-statements

import os
from datetime import date
import argparse
from concurrent.futures import ProcessPoolExecutor
from dateutil.relativedelta import relativedelta
import numpy as np
import yaml
# pylint: disable=import-error
import plot_utils
from ghis2s.s2smetric.metricslib import get_anom
from ghis2s.shared.logging_utils import TaskLogger
# pylint: enable=import-error

task_name = os.environ.get('SCRIPT_NAME')
logger = TaskLogger(task_name,
                    os.getcwd(),
                    f'Running s2splots/plot_weekly_anom.py')

def plot_anoms(fcst_year, fcst_mon, cwd, config, region, anom_type):
    '''
    Plot the anomalies.
    '''
    plotdir_template = cwd + '/s2splots/{:04d}{:02d}/'
    plotdir = plotdir_template.format(fcst_year, fcst_mon)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir, exist_ok=True)

    figure_template = '{}/NMME_plot_{}_{}_weeks1-6_{}.png'

    lead_week = list(range(min(6, config["EXP"]["lead_months"])))
    # Universal setup of plots:
    nrows = 2
    ncols = 3
    domain = plot_utils.dicts('boundary', region)
    subtask = region + ': ' + anom_type
    data_dir = cwd + f'/s2smetric/{fcst_year:04d}{fcst_mon:02d}/'
    cartopy_dir = config['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/'
    if region == 'ARCTIC':
        var_list = ['SWE', 'SnowDepth', 'RZSM', 'TOP40RELSM', 'TOP40ST']
    else:
        var_list = config["POST"]["weekly_vars"]

    for var_name in var_list:
        if anom_type == 'ANOM':
            clabel = 'Anomaly (' + plot_utils.dicts('units', var_name) + ')'
            load_table = 'clim_reanaly'
        else:
            #if var_name == 'RZSM':
            #    continue
            clabel = 'Standardized Anomaly'
            if var_name == 'TOP40ST':
                load_table = 'CB11W_'
            else:
                load_table = 'CB11W'

        under_over = plot_utils.dicts('lowhigh', load_table)
        # READ ANOMALIES
        anom_crop = get_anom(data_dir, var_name, anom_type, domain, weekly=True)
        median_anom = np.median(anom_crop.anom.values, axis=0)
        plot_arr = median_anom[lead_week, ]
        if anom_type == 'SANOM':
            levels=plot_utils.dicts('anom_levels', 'standardized')
            plot_arr = np.where(plot_arr < np.max(levels), plot_arr, np.nan)
            plot_arr = np.where(plot_arr > np.min(levels), plot_arr, np.nan)

        begdate = date(fcst_year, fcst_mon, 2)
        titles = []
        for lead in lead_week:
            enddate = begdate + relativedelta(days=6)
            titles.append(
                var_name + ' '+  begdate.strftime("%Y%m%d") + '-' + enddate.strftime("%Y%m%d")
            )
            begdate += relativedelta(days=7)

        figure = figure_template.format(plotdir, region, var_name, anom_type.lower())
        logger.info(f"Plotting {figure}", subtask=subtask)
        stitle = var_name + ' Forecast'
        if anom_type == 'ANOM':
            anom_minmax = plot_utils.dicts('anom_minmax', var_name)
            plot_utils.contours (anom_crop.lon.values, anom_crop.lat.values, nrows,
                                 ncols, plot_arr, load_table, titles, domain,
                                 figure, under_over,
                                 fscale=1.2, stitle=stitle, clabel=clabel,
                                 min_val=anom_minmax[0], max_val=anom_minmax[1],
                                 cartopy_datadir=cartopy_dir, projection=['polar', 90.])
        else:
            if region == 'ARCTIC':
                plot_utils.contours (anom_crop.lon.values, anom_crop.lat.values, nrows,
                                     ncols, plot_arr, load_table, titles, domain,
                                     figure, under_over,
                                     fscale=1.2, stitle=stitle, clabel=clabel, levels=levels,
                                     cartopy_datadir=cartopy_dir, projection=['polar', 90.])
            else:
                plot_utils.contours (anom_crop.lon.values, anom_crop.lat.values, nrows,
                                     ncols, plot_arr, load_table, titles, domain, figure,
                                     under_over, fscale=0.8, stitle=stitle, clabel=clabel,
                                     levels=levels, cartopy_datadir=cartopy_dir)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--fcst_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--fcst_mon', required=True, help= 'forecast end year')
    parser.add_argument('-c', '--configfile', required=True, help='config file name')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-s', '--sanom', required=False, default=None,
                        help='plot SANOM all variables')

    args = parser.parse_args()
    FCST_YEAR = int(args.fcst_year)
    FCST_MON = int(args.fcst_mon)
    CWD = args.cwd

    # load config file
    with open(args.configfile, 'r', encoding="utf-8") as file:
        config_ = yaml.safe_load(file)

    if args.sanom is None:
        plot_anoms(FCST_YEAR, FCST_MON, CWD, config_, 'ARCTIC', 'ANOM')
        plot_anoms(FCST_YEAR, FCST_MON, CWD, config_, 'ARCTIC', 'SANOM')
    else:
        num_workers = int(os.environ.get('NUM_WORKERS', 7))
        regions = ['GLOBAL', 'AFRICA', 'EUROPE', 'CENTRAL_ASIA', 'SOUTH_EAST_ASIA', 'NORTH_AMERICA',
                   'SOUTH_AMERICA']
        with ProcessPoolExecutor(max_workers=7) as executor:
            futures = []
            for region_ in regions:
                futures.append(executor.submit(plot_anoms, FCST_YEAR, FCST_MON, CWD, config_,
                                               region_, 'SANOM'))

            for future in futures:
                result = future.result()
