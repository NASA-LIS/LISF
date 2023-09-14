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
This script plots OL anomalies at lead times 0,1,2,3,4 months
for a given forecast start month and year. The script consolidated
Abheera Hazra's two scripts Plot_real-time_OUTPUT_AFRICOM_NMME_RT_FCST_anom.py and
Plot_real-time_OUTPUT_AFRICOM_NMME_RT_FCST_sanom.py into a single script.
'''
# pylint: disable=no-value-for-parameter

import os
import calendar
import argparse
import xarray as xr
# pylint: disable=no-name-in-module
from netCDF4 import Dataset
# pylint: enable=no-name-in-module
import numpy as np
import yaml
# pylint: disable=import-error
import plot_utils
# pylint: enable=import-error

def plot_anoms(syear, smonth, cwd_, config_, region, standardized_anomaly = None):
    '''
    This function processes arguments and make plots.
    '''
    # Input file template for straight anomalies:
    infile_template = '{}/{}_ANOM_init_monthly_{:02d}_{:04d}.nc'
    figure_template = '{}/NMME_plot_{}_{}_FCST_anom_g.png'
    if standardized_anomaly == 'Y':
        infile_template = '{}/{}_SANOM_init_monthly_{:02d}_{:04d}.nc'
        figure_template = '{}/NMME_plot_{}_{}_FCST_sanom_g.png'

    plotdir_template = cwd_ + '/s2splots/{:04d}{:02d}/' + '/' + config_["EXP"]["lsmdir"] + '/'
    plotdir = plotdir_template.format(fcst_year, fcst_mon)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    data_dir_template = cwd_ + '/s2smetric/{:04d}{:02d}/metrics_cf/' + config_["EXP"]["lsmdir"] + '/'
    data_dir = data_dir_template.format(fcst_year, fcst_mon)

    lead_month = [0, 1, 2, 3, 4, 5]

    nrows = 2
    ncols = 3
    boundary = plot_utils.dicts('boundary', region)
    vmask = (((upstream_lon >= boundary[2]) & (upstream_lon <= boundary[3])) &
             ((upstream_lat >= boundary[0]) & (upstream_lat <= boundary[1])))

    dlon = downstream_lon[vmask]
    dlat = downstream_lat[vmask]
    ulon = upstream_lon[vmask]
    ulat = upstream_lat[vmask]
    carea = cum_area[vmask]
    var_name = 'Streamflow'

    levels = plot_utils.dicts('anom_levels', var_name)
    if standardized_anomaly == 'Y':
        levels = plot_utils.dicts('anom_levels', 'standardized')
    under_over = ['gray', 'blue']

    infile = infile_template.format(data_dir, '*_' + var_name, smonth, syear)
    print("Reading infile {}".format(infile))

    print(infile)
    anom = xr.open_mfdataset(infile, concat_dim='ens',
                             preprocess=plot_utils.preproc, combine='nested')
    anom_crop = plot_utils.crop(boundary, anom.latitude, anom.longitude, anom)
    median_anom = np.nanmedian(anom_crop.anom.values, axis=0)
    plot_arr = median_anom[lead_month, ]
    figure = figure_template.format(plotdir, region, var_name)
    print(figure)
    titles = []
    for lead in lead_month:
        fcast_month = smonth+lead
        fcast_year = syear
        if fcast_month > 12:
            fcast_month -= 12
            fcast_year = syear+1
        titles.append(calendar.month_abbr[fcast_month] + ', ' + str(fcast_year))
    stitle = var_name + ' Forecast'
    clabel = 'Anomaly (' + plot_utils.dicts('units', var_name) + ')'
    if standardized_anomaly == 'Y':
        clabel = 'Standardized Anomaly'
    plot_utils.google_map(anom_crop.longitude.values, anom_crop.latitude.values, nrows,
                          ncols, plot_arr, 'DROUGHT_INV', titles, boundary, figure, under_over,
                          dlat, dlon, ulat, ulon, carea, fscale=0.8, stitle=stitle,
                          clabel=clabel, levels=levels)
    del anom
    del anom_crop

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

    rnetwork =  Dataset (config['SETUP']['supplementarydir'] + '/s2splots/RiverNetwork_information.nc4', mode='r')
    downstream_lon = np.array (rnetwork.variables['DownStream_lon'][:])
    downstream_lat = np.array (rnetwork.variables['DownStream_lat'][:])
    upstream_lon   = np.array (rnetwork.variables['UpStream_lon'][:])
    upstream_lat   = np.array (rnetwork.variables['UpStream_lat'][:])
    cum_area = np.array (rnetwork.variables['CUM_AREA'][:])

    if config ["EXP"]["DOMAIN"] == 'AFRICOM':
        FONT_SIZE1 = 30
        FONT_SIZE2 = 40
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'FAME')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'FAME', standardized_anomaly = 'Y')

    if config ["EXP"]["DOMAIN"] == 'GLOBAL':
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'AFRICA')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'EUROPE')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'CENTRAL_ASIA')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SOUTH_EAST_ASIA')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'NORTH_AMERICA')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SOUTH_AMERICA')

        plot_anoms(fcst_year, fcst_mon, cwd, config, 'AFRICA', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'EUROPE', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'CENTRAL_ASIA', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SOUTH_EAST_ASIA', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'NORTH_AMERICA', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SOUTH_AMERICA', standardized_anomaly = 'Y')
