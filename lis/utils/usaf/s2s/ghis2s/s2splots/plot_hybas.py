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
This script plots streamflow anomalies along river pathways while using the Google map as a canvas.
- Sarith Mahanama (2023-01-13
'''
# pylint: disable=no-value-for-parameter

import os
import calendar
import argparse
import math
import xarray as xr
# pylint: disable=no-name-in-module
from netCDF4 import Dataset
# pylint: enable=no-name-in-module
import numpy as np
import yaml
# pylint: disable=import-error
import plot_utils
# pylint: enable=import-error

STANDARDIZED_ANOMALY = 'Y'

DEFCOMS = ['INDOPACOM', 'CENTCOM', 'AFRICOM', 'EUCOM', 'SOUTHCOM']

def plot_anoms(syear, smonth, cwd, config, dlon, dlat, ulon, ulat,
               carea, boundary, region, google_path, hybas_mask):
    '''
    This function processes arguments and make plots.
    '''

    infile_template = '{}/{}_ANOM_init_monthly_{:02d}_{:04d}.nc'
    figure_template = '{}/NMME_plot_{}_{}_FCST_anom.png'
    if STANDARDIZED_ANOMALY == 'Y':
        infile_template = '{}/{}_SANOM_init_monthly_{:02d}_{:04d}.nc'
        figure_template = '{}/NMME_plot_{}_{}_FCST_sanom.png'

    plotdir_template = cwd + '/s2splots/{:04d}{:02d}/' + '/' + config["EXP"]["lsmdir"] + '/'
    plotdir = plotdir_template.format(syear, smonth)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    data_dir_template = cwd + '/s2smetric/{:04d}{:02d}/metrics_cf/' + config["EXP"]["lsmdir"] + '/'
    data_dir = data_dir_template.format(syear, smonth)

    lead_month = [0, 1, 2]
    nrows = 1
    ncols = 3

    var_name = "Streamflow"
    if STANDARDIZED_ANOMALY == 'Y':
        levels = plot_utils.dicts('anom_levels', 'standardized')

    infile = infile_template.format(data_dir, '*_' + var_name, smonth, syear)
    print("Reading infile {}".format(infile))

    anom = xr.open_mfdataset(infile, concat_dim='ens',
                             preprocess=plot_utils.preproc, combine='nested')
    anom_crop = plot_utils.crop(boundary, anom.latitude, anom.longitude, anom)
    median_anom = np.nanmedian(anom_crop.anom.values, axis=0)
    plot_arr = median_anom[lead_month, ]
    for i in range (0, len(lead_month)):
        plot_arr[i,:,:] = np.where(hybas_mask > 0, plot_arr[i,:,:],-9999.)

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
    clabel = 'Anomaly (' + plot_utils.dicts('units', var_name) + ')'
    if STANDARDIZED_ANOMALY == 'Y':
        clabel = 'Standardized Anomaly'

    under_over = plot_utils.dicts('lowhigh', 'CB11W')
    cartopy_dir = config['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/'
    plot_utils.google_map(anom_crop.longitude.values, anom_crop.latitude.values, nrows,
                          ncols, plot_arr, 'CB11W', titles, boundary, figure, under_over,
                          dlat, dlon, ulat, ulon, carea, google_path, fscale=0.8, stitle=stitle,
                          clabel=clabel, levels=levels, cartopy_datadir=cartopy_dir)
    del anom
    del anom_crop

def process_domain (fcst_year, fcst_mon, cwd, config, rnetwork, plot_domain):
    ''' processes a single domain GLOBAL or USAF COM '''
    downstream_lon = np.array (rnetwork.variables['DownStream_lon'][:])
    downstream_lat = np.array (rnetwork.variables['DownStream_lat'][:])
    upstream_lon   = np.array (rnetwork.variables['UpStream_lon'][:])
    upstream_lat   = np.array (rnetwork.variables['UpStream_lat'][:])
    cum_area = np.array (rnetwork.variables['CUM_AREA'][:])
    bmask = xr.open_dataset(config['SETUP']['supplementarydir'] + \
                            '/s2splots/HYBAS_' + plot_domain + '.nc4')
    lonv = bmask.lon.values
    latv = bmask.lat.values
    #nr, nc = (latv.size, lonv.size)
    lons, lats = np.meshgrid(lonv, latv)
    bas = 0
    google_path = config['SETUP']['supplementarydir'] + '/s2splots/'
    for bid in bmask.BASIN_ID.values:
        hybas_mask = bmask.basin_mask.values[bas,:,:]
        tx_ = np.ma.compressed(np.ma.masked_where(hybas_mask == 0, lons))
        ty_ = np.ma.compressed(np.ma.masked_where(hybas_mask == 0, lats))
        boundary = [math.floor(ty_.min()), math.ceil(ty_.max()),
                    math.floor(tx_.min()), math.ceil(tx_.max())]
        vmask = (((upstream_lon >= boundary[2]) & (upstream_lon <= boundary[3])) &
             ((upstream_lat >= boundary[0]) & (upstream_lat <= boundary[1])))
        sub_mask = plot_utils.crop(boundary, bmask.lat, bmask.lon, bmask.basin_mask)
        region = "{:10d}".format(bid)
        plot_anoms(fcst_year, fcst_mon, cwd, config, downstream_lon[vmask],
                   downstream_lat[vmask], upstream_lon[vmask],upstream_lat[vmask],
                   cum_area[vmask], boundary, region, google_path, sub_mask.values[bas,:,:])
        bas += 1

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
    exp_domain = config_["EXP"]["DOMAIN"]
    if exp_domain == 'GLOBAL':
        for plot_domain_ in DEFCOMS:
            process_domain (fcst_year_, fcst_mon_, cwd_, config_, rnetwork_, plot_domain_)
    else:
        process_domain (fcst_year_, fcst_mon_, cwd_, config_, rnetwork_, exp_domain)
