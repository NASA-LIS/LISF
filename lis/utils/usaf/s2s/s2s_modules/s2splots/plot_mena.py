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


USAF_COLORS = True

def plot_anoms_basin(syear, smonth, cwd, config, dlon, dlat, ulon, ulat,
                     carea, boundary, region, google_path, hymap_mask):
    '''
    This function processes arguments and make plots.
    '''

    infile_template = '{}/{}_SANOM_init_monthly_{:02d}_{:04d}.nc'
    figure_template = '{}/NMME_plot_{}_{}_basins_sanom.png'

    plotdir_template = cwd + '/s2splots/{:04d}{:02d}/' + '/' + config["EXP"]["lsmdir"] + '/'
    plotdir = plotdir_template.format(syear, smonth)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    data_dir_template = cwd + '/s2smetric/{:04d}{:02d}/metrics_cf/' + config["EXP"]["lsmdir"] + '/'
    data_dir = data_dir_template.format(syear, smonth)

    lead_month = [0, 1, 2, 3, 4, 5]
    nrows = 2
    ncols = 3

    var_name = "Streamflow"
    levels = plot_utils.dicts('anom_levels', 'standardized')
    under_over = ['gray', 'blue']

    infile = infile_template.format(data_dir, '*_' + var_name, smonth, syear)
    print("Reading infile {}".format(infile))

    anom = xr.open_mfdataset(infile, concat_dim='ens',
                             preprocess=plot_utils.preproc, combine='nested')
    anom_crop = plot_utils.crop(boundary, anom.latitude, anom.longitude, anom)
    median_anom = np.nanmedian(anom_crop.anom.values, axis=0)
    plot_arr = median_anom[lead_month, ]
    for i in range (0, len(lead_month)):
        plot_arr[i,:,:] = np.where(hymap_mask >= 1.0e9, plot_arr[i,:,:],-9999.)

    figure = figure_template.format(plotdir, region, var_name)

    titles = []
    for lead in lead_month:
        fcast_month = smonth+lead
        fcast_year = syear
        if fcast_month > 12:
            fcast_month -= 12
            fcast_year = syear+1
        titles.append(calendar.month_abbr[fcast_month] + ', ' + str(fcast_year))
    stitle = var_name + ' Forecast'
    clabel = 'Standardized Anomaly'

    cartopy_dir = config['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/'
    plot_utils.google_map(anom_crop.longitude.values, anom_crop.latitude.values, nrows,
                          ncols, plot_arr, 'CB11W_', titles, boundary, figure, under_over,
                          dlat, dlon, ulat, ulon, carea, google_path, fscale=0.8, stitle=stitle,
                          clabel=clabel, levels=levels, cartopy_datadir=cartopy_dir)
    del anom
    del anom_crop

def plot_anoms(syear, smonth, cwd, config, region, standardized_anomaly = None):
    '''
    This function processes arguments and make plots.
    '''
    #def preproc(ds):
    #    ds = ds.isel(ens=0)
    #    return ds

    # Input file template for straight anomalies:
    infile_template = '{}/{}_ANOM_init_monthly_{:02d}_{:04d}.nc'
    # Figure file template for straight anomalies:
    figure_template = '{}/NMME_plot_{}_{}_FCST_anom.png'

    # Standardized anomaly input and plot file templates (option):
    if standardized_anomaly == 'Y':
        infile_template = '{}/{}_SANOM_init_monthly_{:02d}_{:04d}.nc'
        figure_template = '{}/NMME_plot_{}_{}_FCST_sanom.png'

    plotdir_template = cwd + '/s2splots/{:04d}{:02d}/' + config["EXP"]["lsmdir"] + '/'
    plotdir = plotdir_template.format(syear, smonth)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    # Where input data files are located ("s2smetric dir - CF-convention nc files):
    data_dir_template = cwd + '/s2smetric/{:04d}{:02d}/metrics_cf/' + \
        config["EXP"]["lsmdir"] + '/'
    data_dir = data_dir_template.format(syear, smonth)

    lead_month = [0, 1, 2, 3, 4, 5]

    # Universal setup of plots:
    nrows = 2
    ncols = 3
    domain = plot_utils.dicts('boundary', region)

    for var_name in ['RZSM', 'SFCSM', 'Precip', 'AirT']:
        # Streamflow specifics
        if var_name == 'Streamflow':
            ldtfile = config['SETUP']['supplementarydir'] + '/lis_darun/' + \
                config['SETUP']['ldtinputfile']
            ldt = xr.open_mfdataset(ldtfile)
            ldt_crop = plot_utils.crop(domain, ldt.lat, ldt.lon, ldt)
            mask = ldt_crop.HYMAP_drain_area.values

        levels = plot_utils.dicts('anom_levels', var_name)
        if standardized_anomaly == 'Y':
            levels = plot_utils.dicts('anom_levels', 'standardized')

        # colors defualts
        load_table = plot_utils.dicts('anom_tables', var_name)
        if USAF_COLORS and standardized_anomaly is None:
            if var_name in {'AirT', 'Air-T', 'Air_T'}:
                load_table = '14WT2M'
                levels = plot_utils.dicts('anom_levels', 'Air_T_AF')

            elif var_name == 'Precip':
                load_table = '14WPR'
                levels = plot_utils.dicts('anom_levels','Precip_AF')

        under_over = plot_utils.dicts('lowhigh', load_table)

        # READ ANOMALIES
        infile = infile_template.format(data_dir, '*_' + var_name, smonth, syear)
        print("Reading infile {}".format(infile))
        if var_name == 'Streamflow':
            anom = xr.open_mfdataset(infile, concat_dim='ens',
                                     preprocess=plot_utils.preproc, combine='nested')
        else:
            anom = xr.open_mfdataset(infile, concat_dim='ens', combine='nested')
        anom_crop = plot_utils.crop(domain, anom.latitude, anom.longitude, anom)
        median_anom = np.median(anom_crop.anom.values, axis=0)

        if (var_name in {'AirT', 'Air-T', 'Air_T'}) and \
           USAF_COLORS and standardized_anomaly is None:
            median_anom = median_anom*9./5.
        if var_name == 'Precip' and USAF_COLORS and standardized_anomaly is None:
            median_anom = median_anom/25.4

        if var_name in {'Total-Runoff', 'ET'}:
            if standardized_anomaly is None:
                median_anom = median_anom * 86400.

        plot_arr = median_anom[lead_month, ]
        figure = figure_template.format(plotdir, region, var_name)
        titles = []
        for lead in lead_month:
            if var_name == 'Streamflow':
                plot_arr[lead,] = np.ma.masked_where(mask<=2.0e9, plot_arr[lead,])
            fcast_month = smonth+lead
            fcast_year = syear
            if fcast_month > 12:
                fcast_month -= 12
                fcast_year = syear+1
            titles.append(calendar.month_abbr[fcast_month] + ', ' + str(fcast_year))
        stitle = var_name + ' Forecast'
        clabel = 'Anomaly (' + plot_utils.dicts('units', var_name) + ')'

        if USAF_COLORS and standardized_anomaly is None:
            if var_name in {'AirT', 'Air-T', 'Air_T'}:
                clabel = 'Anomaly (' + plot_utils.dicts('units', 'Air_T_AF') + ')'
            elif var_name == 'Precip':
                clabel = 'Anomaly (' + plot_utils.dicts('units', 'Precip_AF') + ')'

        if standardized_anomaly == 'Y':
            clabel = 'Standardized Anomaly'
        cartopy_dir = config['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/'
        plot_utils.contours (anom_crop.longitude.values, anom_crop.latitude.values, nrows,
                             ncols, plot_arr, load_table, titles, domain, figure, under_over,
                             fscale=0.8, stitle=stitle, clabel=clabel, levels=levels,
                             cartopy_datadir=cartopy_dir)
        del anom
        del anom_crop

def process_domain (fcst_year, fcst_mon, cwd, config, rnetwork, region):
    ''' processes a single domain GLOBAL or USAF COM '''
    downstream_lon = np.array (rnetwork.variables['DownStream_lon'][:])
    downstream_lat = np.array (rnetwork.variables['DownStream_lat'][:])
    upstream_lon   = np.array (rnetwork.variables['UpStream_lon'][:])
    upstream_lat   = np.array (rnetwork.variables['UpStream_lat'][:])
    cum_area = np.array (rnetwork.variables['CUM_AREA'][:])

    google_path = config['SETUP']['supplementarydir'] + '/s2splots/'

    boundary = plot_utils.dicts('boundary', region)
    vmask = (((upstream_lon >= boundary[2]) & (upstream_lon <= boundary[3])) &
         ((upstream_lat >= boundary[0]) & (upstream_lat <= boundary[1])))

    ldtfile = config['SETUP']['supplementarydir'] + '/lis_darun/' + \
        config['SETUP']['ldtinputfile']
    ldt = xr.open_mfdataset(ldtfile)
    ldt_crop = plot_utils.crop(boundary, ldt.lat, ldt.lon, ldt)
    hymap_mask = ldt_crop.HYMAP_drain_area.values

    plot_anoms_basin(fcst_year, fcst_mon, cwd, config, downstream_lon[vmask],
               downstream_lat[vmask], upstream_lon[vmask],upstream_lat[vmask],
                     cum_area[vmask], boundary, region, google_path, hymap_mask)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--fcst_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--fcst_mon', required=True, help= 'forecast end year')
    parser.add_argument('-c', '--configfile', required=True, help='config file name')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')

    args = parser.parse_args()
    configfile = args.configfile
    fcst_year_ = int(args.fcst_year)
    fcst_mon_ = int(args.fcst_mon)
    cwd_ = args.cwd

    # load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        config_ = yaml.safe_load(file)

    rnetwork_ =  Dataset (config_['SETUP']['supplementarydir'] + \
                          '/s2splots/RiverNetwork_information.nc4', mode='r')

    for region_ in ['TUNISIA', 'ME_CRES']:
        process_domain (fcst_year_, fcst_mon_, cwd_, config_, rnetwork_, region_)
        plot_anoms(fcst_year_, fcst_mon_, cwd_, config_, region_)
        plot_anoms(fcst_year_, fcst_mon_, cwd_, config_, region_, standardized_anomaly = 'Y')
