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
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
import cartopy.io.img_tiles as cimgt
import types
from load_colors import load_table
import math
import requests
import PIL
import yaml
import argparse
import sys

USAF_COLORS = True
boundary = [-20, 55, -40, 40]
standardized_anomaly = 'Y'
region = ''
mask_100 = 0


class CachedTiler(object):
    def __init__(self, tiler):
        self.tiler = tiler

    def __getattr__(self, name):
        attr = getattr(self.tiler, name, None)
        if isinstance(attr, types.MethodType):
            attr = types.MethodType(attr.__func__, self)
        return attr

    def get_image(self, tile):
        tileset_name = '{}'.format(self.tiler.__class__.__name__.lower())
        cache_dir = os.path.expanduser(os.path.join('~/', 'image_tiles', tileset_name))
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        tile_fname = os.path.join(cache_dir, '_'.join(str(v) for v in tile) + '.png')
        if not os.path.exists(tile_fname):
            response = requests.get(self._image_url(tile),
                                    stream=True)

            with open(tile_fname, "wb") as fh:
                for chunk in response:
                    fh.write(chunk)
        with open(tile_fname, 'rb') as fh:
            img = PIL.Image.open(fh)
            img = img.convert(self.desired_tile_form)     
        return img, self.tileextent(tile), 'lower'
def plot_anoms(syear, smonth, cwd, config, dlon, dlat, ulon, ulat):
    '''
    This function processes arguments and make plots.
    '''
    def preproc(ds):
        ds = ds.isel(ens=0) 
        return ds

    def getclosest_ij(lats,lons,latpt,lonpt):
        # find squared distance of every point on grid
        dist_sq = (lats-latpt)**2 + (lons-lonpt)**2
        minindex_flattened = dist_sq.argmin()
        return np.unravel_index(minindex_flattened, lats.shape)

    def map2pfaf (anom, lats, lons, levels):
        carr = []
        xx1 = []
        xx2 = []
        yy1 = []
        yy2 = []

        for i, value in enumerate(dlon):
            iy = min(range(len(lats)), key=lambda j: abs(lats[j] - ulat[i]))
            ix = min(range(len(lons)), key=lambda j: abs(lons[j] - ulon[i]))
            au = anom[iy,ix]
            iy = min(range(len(lats)), key=lambda j: abs(lats[j] - dlat[i]))
            ix = min(range(len(lons)), key=lambda j: abs(lons[j] - dlon[i]))
            ad = anom[iy,ix]

            if au and ad:
                this_val = 0.5*(au+ad)
                this_val = min (this_val, max(levels) - 0.01)
                this_val = max (this_val, min(levels) + 0.01)
                k = min(range(len(levels)), key=lambda i: abs(levels[i]-this_val))
                carr.append(k)
                xx1.append(dlon[i])
                xx2.append(ulon[i])
                yy1.append(dlat[i])
                yy2.append(ulat[i])
        return carr, xx1, xx2, yy1, yy2
            
    infile_template = '{}/{}_ANOM_init_monthly_{:02d}_{:04d}.nc'
    figure_template = '{}/NMME_plot_{}_{}_FCST_anom.png'
    if standardized_anomaly == 'Y':
        infile_template = '{}/{}_SANOM_init_monthly_{:02d}_{:04d}.nc'
        figure_template = '{}/NMME_plot_{}_{}_FCST_sanom.png'

    plotdir_template = cwd + '/{:04d}{:02d}/' + '/' + config["EXP"]["lsmdir"] + '/'
    plotdir = plotdir_template.format(fcst_year, fcst_mon)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    data_dir_template = cwd + '/s2smetric/{:04d}{:02d}/metrics_cf/' + '/' + config["EXP"]["lsmdir"] + '/'
    data_dir = data_dir_template.format(fcst_year, fcst_mon)
    
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
    nrows = 1
    ncols = 3
    figwidth = 25
    default_levels = [-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
    standard_levels = [-4.0, -3.0, -2, -1, -0.5, -0.25, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0]    

    mpl.style.use('bmh')
    cbar_axes = [0.2, 0.04, 0.6, 0.03]

    cmap, col_under, col_higher, extend = plt.cm.RdYlGn, 'black', '#B404AE', 'both'

    var_name = "Streamflow"
    levels = get_levels.get(var_name, default_levels)
    if standardized_anomaly == 'Y':
        levels = standard_levels        
        
    style_color = load_table('DROUGHT_INV')
    color_arr = []
    for color in style_color:
        rgb = [float(value) / 255 for value in color]
        color_arr.append(rgb)

    norm = mpl.colors.BoundaryNorm(levels, ncolors=256)
    cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
    cmap.set_under('gray')
    cmap.set_over('blue')
        
    count_plot = 0

    # Initialize our Google Maps tiles
    actual_tiler = cimgt.GoogleTiles()
    imagery = CachedTiler(actual_tiler)
        
    # First plotting Initial conditions
    # ---------------------------------

    fig = plt.figure(figsize=(figwidth,
                              figwidth*(nrows*(boundary[3]
                                               - boundary[2]))
                              /(ncols*(boundary[1]
                                       - boundary[0]))), facecolor='white')
    gs_ = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)
    ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
    infile = infile_template.format(data_dir, '*_' + var_name, smonth, syear)
    print("Reading infile {}".format(infile))
       
    anom = xr.open_mfdataset(infile, concat_dim='ens',preprocess=preproc)
    median_anom = np.nanmedian(anom.anom.values, axis=0)
            
    #for lead in lead_month:
    for lead in [0,1,2]:
        ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
        
        count_plot += 1
        fcast_month = smonth+lead
        fcast_year = syear
        if fcast_month > 12:
            fcast_month -= 12
            fcast_year = syear+1

        ax_.set_extent([boundary[1],boundary[0],boundary[2],boundary[3]], crs=ccrs.Geodetic())
        ax_.add_image(imagery, 9)
                        
        if standardized_anomaly == 'Y':
            sanom = median_anom[lead, ]
            sanom = np.ma.masked_where(mask_100==0, sanom)
            # loop through vector values and overplot
            carr, xx1, xx2, yy1, yy2 = map2pfaf(sanom, anom.latitude.values, anom.longitude.values, levels)
            
            pfaf_cnt = 0
            for this_col in carr:
                rgb = np.array(style_color[:][this_col])/255.
                track = sgeom.LineString(zip([xx1[pfaf_cnt],xx2[pfaf_cnt]], [yy1[pfaf_cnt], yy2[pfaf_cnt]]))
                ax_.add_geometries([track], ccrs.PlateCarree(),facecolor='none', edgecolor=rgb, linewidth=2)
                pfaf_cnt += 1

        else:                        
            cs_ = plt.pcolormesh(anom.longitude.values, anom.latitude.values, median_anom[lead, ],
                                 norm=colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False), cmap=cmap,zorder=3)                

        gl_ = ax_.gridlines(draw_labels=True)
        gl_.xlabels_top = False
        gl_.xlabels_bottom = False
        gl_.ylabels_right = False
        gl_.ylabels_left = False
        if lead == 0:
            plt.text(-0.15, 0.5, var_name + ' Forecast', verticalalignment='center',
                     horizontalalignment='center', transform=ax_.transAxes,
                     color='Blue', fontsize=30, rotation=90)
            gl_.ylabels_left = True
        plt.title(calendar.month_abbr[fcast_month] + ', ' + str(fcast_year), fontsize=40)
        if lead >= 3:
            gl_.xlabels_bottom = True
        if lead == 3:
            plt.text(-0.15, 0.5, var_name + ' Forecast', verticalalignment='center',
                     horizontalalignment='center', transform=ax_.transAxes,
                     color='Blue', fontsize=30, rotation=90)
            gl_.ylabels_left = True
        ax_.coastlines()
        ax_.add_feature(cfeature.BORDERS)
        ax_.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')
        cax = fig.add_axes(cbar_axes)
        
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False), cmap=cmap), cax=cax, orientation='horizontal', ticks=levels)

        if standardized_anomaly == 'Y':
            cbar.set_label('Standardized Anomaly', fontsize=30)
        cbar.ax.tick_params(labelsize=20)
        
        figure = figure_template.format(plotdir, region, var_name)
        print(figure)
        plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')

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

    RNetWork =  Dataset ('/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/SRTM-TopoData/RiverNetwork_information.nc4',mode='r')
    downstream_lon = np.array (RNetWork.variables['DownStream_lon'][:])
    downstream_lat = np.array (RNetWork.variables['DownStream_lat'][:])
    upstream_lon   = np.array (RNetWork.variables['UpStream_lon'][:])
    upstream_lat   = np.array (RNetWork.variables['UpStream_lat'][:])    

    bmask = xr.open_dataset('/discover/nobackup/projects/usaf_lis/smahanam/IDL/AFRICOM.nc4')
    lons = bmask.longitude.values
    lats = bmask.latitude.values
    
    b = 0
    for bid in bmask.HYBAS_ID.values:
        mask_100 = bmask.basin_mask.values[b,]
        print (mask_100.sum())
        tx = np.ma.compressed(np.ma.masked_where(mask_100 == 0, lons))
        ty = np.ma.compressed(np.ma.masked_where(mask_100 == 0, lats))
        boundary = [math.floor(tx.min()), math.ceil(tx.max()), math.floor(ty.min()), math.ceil(ty.max())]
        vmask = (((upstream_lon >= boundary[0]) & (upstream_lon <= boundary[1])) &
             ((upstream_lat >= boundary[2]) & (upstream_lat <= boundary[3])))
        region = bid.decode("utf-8")
        plot_anoms(fcst_year, fcst_mon, cwd, config, downstream_lon[vmask], downstream_lat[vmask], upstream_lon[vmask], upstream_lat[vmask])
        b += 1

