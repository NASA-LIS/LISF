#!/usr/bin/env python
'''
This script plots OL anomalies at lead times 0,1,2,3,4 months
for a given forecast start month and year. The script consolidated
Abheera Hazra's two scripts Plot_real-time_OUTPUT_AFRICOM_NMME_RT_FCST_anom.py and
Plot_real-time_OUTPUT_AFRICOM_NMME_RT_FCST_sanom.py into a single script.
'''
# pylint: disable=no-value-for-parameter

from __future__ import division
import os
import calendar
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import xarray as xr
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from load_colors import load_table
import yaml
import argparse

USAF_COLORS = True
FONT_SIZE1 = 16
FONT_SIZE2 = 24

def plot_anoms(syear, smonth, cwd, config, region, standardized_anomaly = None):
    '''
    This function processes arguments and make plots.
    '''
    def preproc(ds):
        ds = ds.isel(ens=0) 
        return ds
    
    # Input file template for straight anomalies:
    infile_template = '{}/{}_ANOM_init_monthly_{:02d}_{:04d}.nc'
    # Figure file template for straight anomalies:
    figure_template = '{}/NMME_plot_{}_{}_FCST_anom.png'

    # Standardized anomaly input and plot file templates (option):
    if standardized_anomaly == 'Y':
        infile_template = '{}/{}_SANOM_init_monthly_{:02d}_{:04d}.nc'
        figure_template = '{}/NMME_plot_{}_{}_FCST_sanom.png'

    plotdir_template = cwd + '/s2splots/{:04d}{:02d}/' + config["EXP"]["lsmdir"] + '/'
    plotdir = plotdir_template.format(fcst_year, fcst_mon)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    # Where input data files are located ("s2smetric dir - CF-convention nc files):
    data_dir_template = cwd + '/s2smetric/{:04d}{:02d}/metrics_cf/' + config["EXP"]["lsmdir"] + '/'
    data_dir = data_dir_template.format(fcst_year, fcst_mon)
    
    get_boundary = {
        'EA': (22, 55, -12, 23),
        'WA': (-19, 26, -5, 25),
        'SA': (8, 52, -37, 6),
        'SA1':(24, 33, -31, -24),
        'FAME':(-20, 55, -40, 40)
        'EUROPE': (-10, 50, 35, 70),
        'CENTRAL_ASIA': (35, 90, 14, 55),
        'SOUTH_EAST_ASIA': (75, 140, -15, 30),
        'SOUTH_AMERICA':(-85, -35, -55, 10),
        'NORTH_AMERICA':(-165, -58, 9, 72),
        'AFRICA':(-20, 55, -40, 40),
        'GLOBAL':(-179, 179, -89, 89)        
    }    
    # Anomaly unit ranges for each variable:
    get_levels = {
        'Streamflow': [-10000, -900, -600, -300, -100, -25, 25, 100, 300, 600, 900, 10000],
        'Total-Runoff': [-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        'Precip_AF': [-10, -6, -4, -2, -1, -0.5, -0.25,0.25, 0.5, 1., 2., 4., 6., 10.],
        'Precip': [-10, -6, -4, -2, -1, -0.25, 0.25, 1., 2., 4., 6., 10.],
        'Air-T': [-4., -3., -2., -1., -0.5, -0.25, 0.25, 0.5, 1., 2., 3., 4.],
        'Air_T': [-4., -3., -2., -1., -0.5, -0.25, 0.25, 0.5, 1., 2., 3., 4.],
        #'Air_T_AF': [-15., -12., -9., -6., -3., 3., 6., 9., 12., 15.],
        'Air_T_AF': [-5., -4., -3., -2., -1., -0.5, 0.5, 1., 2., 3., 4., 5.],
        'TWS': np.array([-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])*500.,
        'ET': np.array([-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])*2.
    }
    get_units = {
        'Streamflow': 'm^3/s',
        'Precip': 'mm/d',
        'Precip_AF': 'in/mon',
        'Air-T': 'K',
        'Air_T': 'K',
        'Air_T_AF': 'F',
        'TWS': 'mm',
        'ET': 'mm/d'
    }
    lead_month = [0, 1, 2, 3, 4, 5]
    
    # Plotting Parameters
    add_land = True
    add_rivers = True
    resol = '50m'  # use data at this scale
    bodr = cartopy.feature.NaturalEarthFeature(category='cultural',
                                               name='admin_0_boundary_lines_land',
                                               scale=resol, facecolor='none', alpha=0.7)
    coastlines = cartopy.feature.NaturalEarthFeature('physical', 'coastline',
                                                     scale=resol, edgecolor='black', facecolor='none')
    land = cartopy.feature.NaturalEarthFeature('physical', 'land', scale=resol, edgecolor='k',
                                               facecolor=cfeature.COLORS['land'])
    ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='none',
                                               facecolor=cfeature.COLORS['water'])
    lakes = cartopy.feature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='b',
                                                facecolor=cfeature.COLORS['water'])
    rivers = cartopy.feature.NaturalEarthFeature('physical', 'rivers_lake_centerlines',
                                                 scale=resol, edgecolor='b', facecolor='none')
    # Universal setup of plots:
    nrows = 2
    ncols = 3
    figwidth = 25
    default_levels = [-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
    standardized_levels = [-4.0, -3.0, -2, -1, -0.5, -0.25, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0]    

    mpl.style.use('bmh')
    cbar_axes = [0.2, 0.04, 0.6, 0.03]

    mask_file100 = config["PLOTS"]["mask_file100"]
    cmap, col_under, col_higher, extend = plt.cm.RdYlGn, 'black', '#B404AE', 'both'
    
    for var_name in config["POST"]["metric_vars"]:
        #if var_name == 'Streamflow':
        #    continue
        levels = get_levels.get(var_name, default_levels)
        if standardized_anomaly == 'Y':
            levels = standardized_levels

        if var_name == 'Air-T' or var_name == 'Air_T':
            #print ('CMAPR', var_name)
            style_color = load_table('L11W')
            if USAF_COLORS and standardized_anomaly is None:
                style_color = load_table('14WT2M')
            color_arr = []
            for color in style_color:
                rgb = [float(value) / 255 for value in color]
                color_arr.append(rgb)

            # normalize bound values
            if USAF_COLORS and standardized_anomaly is None:
                norm = mpl.colors.BoundaryNorm(get_levels.get('Air_T_AF'), ncolors=256)
                levels = get_levels.get('Air_T_AF', default_levels)
            else:
                norm = mpl.colors.BoundaryNorm(get_levels.get('Air_T'), ncolors=256)

            # create a colormap
            cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
            cmap.set_under(col_under)
            cmap.set_over(col_higher)
            if USAF_COLORS and standardized_anomaly is None:
                cmap.set_over('gray')

        elif var_name == 'Precip' and USAF_COLORS and standardized_anomaly is None:
            style_color = load_table('14WPR')
            color_arr = []
            for color in style_color:
                rgb = [float(value) / 255 for value in color]
                color_arr.append(rgb)
            norm = mpl.colors.BoundaryNorm(get_levels.get('Precip_AF'), ncolors=256)
            levels = get_levels.get('Precip_AF', default_levels)
            cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
            cmap.set_under(col_under)
            cmap.set_over(col_higher)
            if USAF_COLORS and standardized_anomaly is None:
                cmap.set_under('gray')
                cmap.set_over('blue')
        
        else:
            #print ('CMAP')
            cmap = mpl.cm.get_cmap("RdYlGn").copy()
            extend = 'both'
            col_under = 'black'
            col_higher = '#B404AE'
            cmap.set_under(col_under)
            cmap.set_over(col_higher)
    

        count_plot = 0

        # First plotting Initial conditions
        # ---------------------------------

        fig = plt.figure(figsize=(figwidth,
                                  figwidth*(nrows*(get_boundary.get(region)[3]
                                                   - get_boundary.get(region)[2]))
                                  /(ncols*(get_boundary.get(region)[1]
                                           - get_boundary.get(region)[0]))))
        gs_ = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)
        ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
        infile = infile_template.format(data_dir, '*_' + var_name, smonth, syear)
        print("Reading infile {}".format(infile))
        
        if var_name == 'Streamflow':
            print(infile)
            anom = xr.open_mfdataset(infile, concat_dim='ens',preprocess=preproc, combine='nested')

            median_anom = np.nanmedian(anom.anom.values, axis=0)
            #if standardized_anomaly is None:
            mask1 = xr.open_mfdataset(mask_file100)
            mask100=mask1.SFMASK.values
        else:
            anom = xr.open_mfdataset(infile, concat_dim='ens', combine='nested')
            median_anom = np.median(anom.anom.values, axis=0)

        if var_name == 'Air_T' and USAF_COLORS and standardized_anomaly is None:
            median_anom = median_anom*9./5.
        if var_name == 'Precip' and USAF_COLORS and standardized_anomaly is None:
            median_anom = median_anom*30./25.4
            
        if var_name == 'Total-Runoff' or var_name == 'ET':
            if standardized_anomaly is None:
                median_anom = median_anom * 86400.

        for lead in lead_month:
            ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
            count_plot += 1
            fcast_month = smonth+lead
            fcast_year = syear
            if fcast_month > 12:
                fcast_month -= 12
                fcast_year = syear+1
            ax_.set_extent(get_boundary.get(region), ccrs.PlateCarree())

            # Streamflow variable:
            if var_name == 'Streamflow':
                if standardized_anomaly == 'Y':
                    sanom = median_anom[lead, ]
                    sanom = np.ma.masked_where(mask100==0, sanom)
                    cs_ = plt.pcolormesh(anom.longitude.values, anom.latitude.values, sanom,
                                     norm=colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False), cmap=cmap,zorder=3)
                else:                        
                    cs_ = plt.pcolormesh(anom.longitude.values, anom.latitude.values, median_anom[lead, ],
                                     norm=colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False), cmap=cmap,zorder=3)                
            else:
                if var_name == 'Air-T' or var_name == 'Air_T':
                    # Air-T reversed colors red (warm) green (cold)
                    cs_ = plt.contourf(anom.longitude.values, anom.latitude.values, median_anom[lead, ],
                                   levels, cmap=cmap, extend=extend, transform=ccrs.PlateCarree())
                else:
                    cs_ = plt.contourf(anom.longitude.values, anom.latitude.values, median_anom[lead, ],
                                   levels, cmap=cmap, extend=extend, transform=ccrs.PlateCarree())

            gl_ = ax_.gridlines(draw_labels=True)
            # Updated for Python 3.9
            gl_.top_xlabels = False
            gl_.bottom_xlabels = False
            gl_.right_labels = False
            gl_.left_labels = False

            # First row of plot panels (lead-months 0-2)
            if lead == 0:
                plt.text(-0.15, 0.5, var_name + ' Forecast', verticalalignment='center',
                         horizontalalignment='center', transform=ax_.transAxes,
                         color='Blue', fontsize=FONT_SIZE1, rotation=90)
                gl_.left_labels = True
            plt.title(calendar.month_abbr[fcast_month] + ', ' + str(fcast_year), fontsize=FONT_SIZE2)
            if lead >= 3:
                gl_.bottom_labels = True
            if lead == 3:
                plt.text(-0.15, 0.5, var_name + ' Forecast', verticalalignment='center',
                         horizontalalignment='center', transform=ax_.transAxes,
                         color='Blue', fontsize=FONT_SIZE1, rotation=90)
                gl_.left_labels = True
                
            ax_.coastlines()
            ax_.add_feature(cfeature.BORDERS)
            ax_.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')
            cax = fig.add_axes(cbar_axes)

            if var_name == 'Streamflow' and standardized_anomaly == 'Y':
                ax_.add_feature(land, facecolor='lightgrey',zorder=2)
            cbar = fig.colorbar(cs_, cax=cax, orientation='horizontal', ticks=levels)

            if USAF_COLORS and standardized_anomaly is None:
                if var_name == 'Air_T':
                    cbar.set_label('Anomaly (' + get_units.get('Air_T_AF', 'm^3/m^3') + ')', fontsize=FONT_SIZE2)
                elif var_name == 'Precip':
                    cbar.set_label('Anomaly (' + get_units.get('Precip_AF', 'm^3/m^3') + ')', fontsize=FONT_SIZE2)
            else:
                cbar.set_label('Anomaly (' + get_units.get(var_name, 'm^3/m^3') + ')', fontsize=FONT_SIZE2)
            if standardized_anomaly == 'Y':
                cbar.set_label('Standardized Anomaly', fontsize=FONT_SIZE2)
            cbar.ax.tick_params(labelsize=20)

            figure = figure_template.format(plotdir, region, var_name)
            print(figure)
            plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')
            #plt.show()

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
    with open(configfile, 'r') as file:
        config = yaml.safe_load(file)

    if config ["EXP"]["domain"] == 'AFRICOM':
        FONT_SIZE1 = 30
        FONT_SIZE2 = 40
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'FAME')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'WA'  )
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'EA'  )
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SA'  )
        
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'FAME', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'WA'  , standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'EA'  , standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SA'  , standardized_anomaly = 'Y')

    if config ["EXP"]["domain"] == 'GLOBAL':

        plot_anoms(fcst_year, fcst_mon, cwd, config, 'GLOBAL')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'AFRICA')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'EUROPE')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'CENTRAL_ASIA')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SOUTH_EAST_ASIA')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'NORTH_AMERICA')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SOUTH_AMERICA')
        
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'GLOBAL', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'AFRICA', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'EUROPE', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'CENTRAL_ASIA', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SOUTH_EAST_ASIA', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'NORTH_AMERICA', standardized_anomaly = 'Y')
        plot_anoms(fcst_year, fcst_mon, cwd, config, 'SOUTH_AMERICA', standardized_anomaly = 'Y')

    # set correct permission to directories and files 
    E2ESDIR = config["SETUP"]["E2ESDIR"]
    SPCODE = config["SETUP"]["SPCODE"]
    os.chdir(E2ESDIR)
    command = "chgrp -R " + SPCODE + ' ' + E2ESDIR
    res = os.system(command)
    command = "find . -type d -exec chmod 0775 {} \;"
    res = os.system(command)
    command = 'find ' + E2ESDIR + '/. -name "*.nc" -type f -exec chmod 0664 {} \;'
    res = os.system(command)
    command = 'find ' + E2ESDIR + '/. -name "*.NC" -type f -exec chmod 0664 {} \;'
    res = os.system(command)
    command = 'find ' + E2ESDIR + '/. -name "*.TIF" -type f -exec chmod 0664 {} \;'
    res = os.system(command)
    command = 'find ' + E2ESDIR + '/. -name "*.png" -type f -exec chmod 0664 {} \;'
    res = os.system(command)    
    os.chdir(cwd)
