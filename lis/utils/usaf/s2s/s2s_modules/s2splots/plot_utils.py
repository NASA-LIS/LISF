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
plotting functions:
(1) dicts contains all domain boundaries, colorbar levels, untis for all plotting variables
(2) figure_size: sets plotting figure size depends on figure width, number of rows and colummns,
(3) load_table: custom color tables using RGB values
(4) compute_daius: computes the radius of the marker to plot stations on the map
(5) preproc is used to read only the 1st ensemble member for streamflow anomaly processing
(6) crop: crops a xarray dataset
(7) map2pfaf: maps LIS grid cells to watersheds to plot river pathway
(8) CachedTiler an object to pull Google tiles for Cartopy to use as the canvas
(9) contours : plots contours
(10) google_map : plots river path ways on the Google map
(11) stations : plots values using colors at each station depcting a small circle on a map.
Sarith Mahanama 2023-01-13
'''

import os
import types
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib import colors
from matplotlib import patches
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import shapely.geometry as sgeom
import requests
import PIL
import numpy as np
mpl.use('pdf')

FONT_SIZE1 = 16
FONT_SIZE2 = 18
FONT_SCALE = 1.

# Plotting Parameters
ADD_LAND = True
ADD_RIVERS = True
RESOL = '50m'  # use data at this scale
FIGWIDTH = 25
cbar_axes = [0.15, 0.04, 0.7, 0.02]

mpl.use('pdf')
mpl.style.use('bmh')
COL_UNDER = 'black'
COL_HIGHER = '#B404AE'
EXTEND = 'both'
PLOT_CIRCLE_RADIUS = 0.5

def dicts(dic_name, key):
    ''' shares dictionaries for all plotting programs '''
    boundary = {
        'EA': (-12, 23, 22, 55),
        'WA': (-5, 25, -19, 26),
        'SA': (-37, 6, 8, 52),
        'SA1':(-31, -24, 24, 33),
        'FAME':(-40, 40, -20, 55),
        'EUROPE': (35, 70, -10, 50),
        'CENTRAL_ASIA': (5, 55, 33, 95),
        'SOUTH_EAST_ASIA': (-10, 50, 90, 150),
        'SOUTH_AMERICA':(-58, 16, -85, -33),
        'NORTH_AMERICA':(9, 72, -165, -58),
        'AFRICA':(-40, 40, -20, 60),
        'GLOBAL':(-89, 89, -179, 179),
        'TUNISIA':(30, 38, 7, 12),
        'ME_CRES':(21, 39, 24, 62),
        'KOREA':(30, 45, 115, 140),
        'ROMANIA':(40, 55, 18, 45)
    }
    # Anomaly unit ranges for each variable:
    anom_levels = {
        'Streamflow': [-10000, -900, -600, -300, -100, -25, 25, 100, 300, 600, 900, 10000],
        'Total-Runoff': [-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        'Precip_AF': [-10, -6, -4, -2, -1, -0.5, -0.25,0.25, 0.5, 1., 2., 4., 6., 10.],
        'Precip': [-10, -6, -4, -2, -1, -0.25, 0.25, 1., 2., 4., 6., 10.],
        'Air-T': [-4., -3., -2., -1., -0.5, -0.25, 0.25, 0.5, 1., 2., 3., 4.],
        'Air_T': [-4., -3., -2., -1., -0.5, -0.25, 0.25, 0.5, 1., 2., 3., 4.],
        'AirT': [-4., -3., -2., -1., -0.5, -0.25, 0.25, 0.5, 1., 2., 3., 4.],
        'Air_T_AF': [-5., -4., -3., -2., -1., -0.5, 0.5, 1., 2., 3., 4., 5.],
        'TWS': np.array([-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])*500.,
        'ET': np.array([-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])*2.,
        'standardized': [-4.0, -3.0, -2, -1, -0.5, -0.25, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0]
    }
    units = {
        'Streamflow': 'm^3/s',
        'Precip': 'mm/d',
        'Precip_AF': 'in/mon',
        'Air-T': 'K',
        'Air_T': 'K',
        'AirT': 'K',
        'Air_T_AF': 'F',
        'TWS': 'mm',
        'ET': 'mm/d'
    }
    anom_tables = {
        'Streamflow': 'CB11W',
        'Precip': 'CB11W',
        'Precip_AF': 'CB11W',
        'Air-T': 'CB11W_',
        'Air_T': 'CB11W_',
        'AirT': 'CB11W_',
        'Air_T_AF': 'CB11W',
        'TWS': 'CB11W',
        'ET': 'CB11W'
    }
    lowhigh = {
        'CB11W': ['black','indigo'],
        'CB11W_': ['indigo','black'],
        '14WT2M': ['#68228B','black'],     # darkorchid4','black']
        '14WPR': ['#8B4500','blueviolet'], # darkorange4','blueviolet']
        'L21': ['blueviolet','#8B4500']
    }

    default_levels = [-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
    default_units = 'm^3/m^3'
    if dic_name == 'boundary':
        ret = boundary.get(key)
    elif dic_name == 'anom_levels':
        ret = anom_levels.get(key, default_levels)
    elif dic_name == 'units':
        ret = units.get(key, default_units)
    elif dic_name == 'anom_tables':
        ret = anom_tables.get(key, 'CB11W')
    elif dic_name == 'lowhigh':
        ret = lowhigh.get(key)
    return ret

def cartopy_dir(data_dir):
    '''
    configuring cartopy
    '''
    global BODR, COASTLINES, land, OCEAN, lakes, rivers
    if data_dir is not None:
        cartopy.config['data_dir'] = data_dir

    BODR = cartopy.feature.NaturalEarthFeature(category='cultural',
                                               name='admin_0_boundary_lines_land',
                                               scale=RESOL, facecolor='none', alpha=0.7)
    COASTLINES = cartopy.feature.NaturalEarthFeature('physical', 'coastline',
                                                     scale=RESOL, edgecolor='black',
                                                     facecolor='none')
    land = cartopy.feature.NaturalEarthFeature('physical', 'land', scale=RESOL, edgecolor='k',
                                               facecolor=cfeature.COLORS['land'])
    OCEAN = cartopy.feature.NaturalEarthFeature('physical', 'ocean', scale=RESOL, edgecolor='none',
                                                facecolor=cfeature.COLORS['water'])
    lakes = cartopy.feature.NaturalEarthFeature('physical', 'lakes', scale=RESOL, edgecolor='b',
                                                facecolor=cfeature.COLORS['water'])
    rivers = cartopy.feature.NaturalEarthFeature('physical', 'rivers_lake_centerlines',
                                             scale=RESOL, edgecolor='b', facecolor='none')


def figure_size(figwidth, domain, nrows, ncols):
    ''' defines plotting figure size'''
    this_size = (figwidth,
                 figwidth*(nrows*(domain[1] - domain[0]))/
                 (ncols*(domain[3] - domain[2])))
    return this_size

def load_table (table_key):
    '''
    loads RGB based color tables
    some good plettes can be found:
    https://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
    '''

    tables = {
        'DROUGHT':[[  0,  0,  0],
                   [  0,115,  0],
                   [  0,159,  0],
                   [  0,210,  0],
                   [ 47,255, 67],
                   [200,255,130],
                   [255,255,255],
                   [255,255,  0],
                   [255,219,  0],
                   [255,157,  0],
                   [249,  0,  0],
                   [197,  0,  0]],
        'DROUGHT_INV':[[167,  0,  0],
                       [197,  0,  0],
                       [249,  0,  0],
                       [255,157,  0],
                       [255,205,  0],
                       [255,255,  0],
                       [247,250,200],
                       [200,255,130],
                       [ 47,255, 67],
                       [  0,210,  0],
                       [  0,159,  0],
                       [  0,115,  0],
                       [ 28, 64, 33]],
        'L21':[[  0,  0,130],
               [  0,  0,200],
               [  0,  0,255],
               [  0, 83,255],
               [  0,115,255],
               [  0,167,255],
               [  0,195,255],
               [  0,227,255],
               [  0,255,255],
               [ 55,255,199],
               [120,255,135],
               [190,255, 67],
               [240,255, 15],
               [255,219,  0],
               [255,187,  0],
               [255,159,  0],
               [255,131,  0],
               [255, 51,  0],
               [233, 23,  0],
               [197,  0,  0],
               [158,  0,  0]],
        'L21W':[[  0,  0,130],
                [  0,  0,200],
                [  0,  0,255],
                [  0, 83,255],
                [  0,115,255],
                [  0,167,255],
                [  0,195,255],
                [  0,227,255],
                [  0,255,255],
                [255,255,255],
                [255,255,255],
                [190,255, 67],
                [240,255, 15],
                [255,219,  0],
                [255,187,  0],
                [255,159,  0],
                [255,131,  0],
                [255, 51,  0],
                [233, 23,  0],
                [197,  0,  0],
                [158,  0,  0]],
        'CB11W':[[103,0,31],
                [178,24,43],
                [214,96,77],
                [244,165,130],
                [253,219,199],
                [247,247,247],
                [209,229,240],
                [146,197,222],
                [67,147,195],
                [33,102,172],
                [5,48,97]],
        'L11W':[[  0,  0,130],
                [  0,  0,255],
                [  0,115,255],
                [  0,195,255],
                [  0,255,255],
                [255,255,255],
                [190,255, 67],
                [255,187,  0],
                [255,131,  0],
                [233, 23,  0],
                [158,  0,  0]],
        'TYPE1':[[255,245,215],
                 [106, 91,154],
                 [202,178,214],
                 [251,154,153],
                 [  0, 85,  0],
                 [ 29,115,  0],
                 [ 77,145,  0],
                 [109,165,  0],
                 [142,185, 13],
                 [233, 23,  0],
                 [255,131,  0],
                 [255,131,200],
                 [255,191,  0],
                 [127, 39,  4],
                 [164, 53,  3],
                 [164, 53,200],
                 [217, 72,  1],
                 [217, 72,200],
                 [204,204,204],
                 [104,104,200],
                 [  0, 70,200]],
        'TYPE2':[[233, 23,  0],
                 [255,131,  0],
                 [255,191,  0],
                 [255,255,178],
                 [210,255,255],
                 [  0,255,255],
                 [  0,155,255],
                 [  0,  0,200],
                 [204,204,204],
                 [170,240,240],
                 [255,255,100],
                 [220,240,100],
                 [205,205,102],
                 [  0,100,  0],
                 [  0,160,  0],
                 [170,200,  0],
                 [  0, 60,  0],
                 [ 40,100,  0],
                 [120,130,  0],
                 [140,160,  0],
                 [190,150,  0],
                 [150,100,  0],
                 [255,180, 50],
                 [255,235,175],
                 [  0,120, 90],
                 [  0,150,120],
                 [  0,220,130],
                 [195, 20,  0],
                 [255,245,215],
                 [  0, 70,200]],
        'GREEN':[[200,255,200],
                 [150,255,150],
                 [ 47,255, 67],
                 [ 60,230, 15],
                 [  0,219,  0],
                 [  0,187,  0],
                 [  0,159,  0],
                 [  0,131,  0]],
        'BLUE':[[55,255,199],
                [ 0,255,255],
                [ 0,227,255],
                [ 0,195,255],
                [ 0,167,255],
                [ 0,115,255],
                [ 0, 83,255],
                [ 0,  0,255],
                [ 0,  0,200],
                [ 0,  0,130]],
        'RED':[[255,255,153],
               [240,255, 15],
               [255,219,  0],
               [255,187,  0],
               [255,159,  0],
               [255,131,  0],
               [255, 51,  0],
               [233, 23,  0],
               [197,  0,  0]],
        'GREY':[[245,245,245],
                [225,225,225],
                [205,205,205],
                [185,185,185],
                [165,165,165],
                [145,145,145],
                [125,125,125],
                [105,105,105],
                [ 85, 85, 85]],
        '14WPR':[[197,  0,  0],
                 [233, 23,  0],
                 [255, 51,  0],
                 [255,131,  0],
                 [255,159,  0],
                 [255,187,  0],
                 [255,255,255],
                 [  0,167,255],
                 [  0,115,255],
                 [  0, 83,255],
                 [  0,  0,255],
                 [  0,  0,200],
                 [  0,  0,130]],
        '14WT2M':[[179, 66,245],
                  [  0,  0,255],
                  [  0,115,255],
                  [  0,195,255],
                  [  0,255,255],
                  [255,255,255],
                  [255,187,  0],
                  [255,131,  0],
                  [255,  0,  0],
                  [128, 48,  9],
                  [196,159,128]],
        'NCV_gebco':[[ 18,  10,  58],
                     [ 23,  49, 110],
                     [ 19,  89, 140],
                     [ 26, 103, 164],
                     [ 30, 114, 178],
                     [ 29, 139, 196],
                     [ 26, 165, 210],
                     [ 27, 184, 223],
                     [ 26, 204, 235],
                     [ 26, 216, 241],
                     [ 38, 223, 241],
                     [ 49, 229, 235],
                     [104, 242, 233],
                     [160, 255, 229],
                     [195, 209,  80],
                     [225, 224, 102],
                     [223, 196,  91],
                     [210, 178,  81],
                     [189, 150,  46],
                     [163, 127,  46],
                     [153, 118,  43],
                     [142, 109,  38],
                     [134, 103,  36],
                     [116,  88,  29]],
        'VERIFY' :[[255,187,  0],
                   [  0,115,255],
                   [255,  0,  0]],
        'srip_reanalysis':[[226,  31,  38],
                           [246, 153, 153],
                           [ 41,  95, 138],
                           [ 95, 152, 198],
                           [175, 203, 227],
                           [114,  59, 122],
                           [173, 113, 181],
                           [214, 184, 218],
                           [245, 126,  32],
                           [253, 191, 110],
                           [236,   0, 140],
                           [247, 153, 209],
                           [  0, 174, 239],
                           [ 96, 200, 232],
                           [ 52, 160,  72],
                           [179,  91,  40],
                           [255, 215,   0],
                           [  0,   0,   0],
                           [119, 119, 119]],
        'cb_9step':[[255,   0,   0],
                    [255, 128,   0],
                    [255, 255,   0],
                    [  0, 255,   0],
                    [  0,   0, 255],
                    [128,   0, 255],
                    [219, 219, 255],
                    [194, 194, 250],
                    [158, 158, 247],
                    [130, 130, 255],
                    [ 97,  97, 255],
                    [ 64,  64, 232],
                    [  0,   0, 194],
                    [  0,   0, 148],
                    [222, 250, 245],
                    [194, 245, 237],
                    [156, 230, 217],
                    [112, 204, 191],
                    [ 43, 184, 163],
                    [  0, 156, 133],
                    [  0, 120, 102],
                    [  0,  92,  79],
                    [219, 255, 219],
                    [186, 245, 186],
                    [140, 235, 140],
                    [ 92, 209,  92],
                    [  0, 184,   0],
                    [  0, 145,   0],
                    [  0, 105,   0],
                    [  0,  77,   0],
                    [235, 204, 255],
                    [222, 176, 255],
                    [199, 148, 237],
                    [186, 112, 237],
                    [171,  77, 237],
                    [138,  51, 199],
                    [107,   0, 186],
                    [ 84,   0, 145],
                    [250, 227, 240],
                    [247, 204, 230],
                    [245, 173, 214],
                    [240, 138, 194],
                    [217,  92, 163],
                    [189,   0, 130],
                    [153,   0, 107],
                    [117,   0,  82],
                    [255, 219, 219],
                    [255, 189, 189],
                    [255, 145, 145],
                    [250,  97,  97],
                    [214,  26,  26],
                    [163,   0,   0],
                    [125,   0,   0],
                    [ 92,   0,   0],
                    [255, 252, 214],
                    [252, 242, 168],
                    [252, 237, 128],
                    [227, 209,   0],
                    [199, 186,  43],
                    [161, 150,   0],
                    [120, 112,   0],
                    [ 84,  82,   0],
                    [255, 222, 199],
                    [252, 199, 161],
                    [250, 176, 125],
                    [232, 143,  79],
                    [209, 105,  31],
                    [186,  77,   0],
                    [153,  64,   0],
                    [115,  48,   0],
                    [240, 240, 240],
                    [222, 222, 222],
                    [199, 199, 199],
                    [171, 171, 171],
                    [145, 145, 145],
                    [120, 120, 120],
                    [ 94,  94,  94],
                    [ 74,  74,  74]],        
        }

    # add a few 24-level monotone colors
    start_brown = np.array([204, 153, 102])
    end_brown = np.array([51, 25, 0]) 
    start_green = np.array([255, 255, 152])
    end_green = np.array([0, 51, 0]) 
    
    num_steps = 24
    gradient_brown = np.linspace(start_brown, end_brown, num_steps, dtype=int)
    gradient_green = np.linspace(start_green, end_green, num_steps, dtype=int)

    rgb_list = []
    for rgb in gradient_brown:
        rgb_list.append([rgb[0],rgb[1], rgb[2]])
    tables['mono_brown'] = rgb_list

    rgb_list = []
    for rgb in gradient_green:
        rgb_list.append([rgb[0],rgb[1], rgb[2]])
    tables['mono_green'] = rgb_list

    if table_key[-1] == '_':
        ct_ = tables[table_key[:-1]]
        ct_.reverse()
    else:
        ct_ = tables[table_key]
    return ct_

def compute_radius(ortho, radius_degrees, lon, lat):
    ''' compute radius '''
    phi1 = lat + radius_degrees if lat <= 0 else lat - radius_degrees
    _, _y1 = ortho.transform_point(lon, phi1, ccrs.PlateCarree())
    return abs(_y1)

def preproc(ds_):
    ''' select only the 1st ensemble for streamflow '''
    ds_ = ds_.isel(ens=0)
    return ds_

def crop (limits, lat, lon, xrin):
    ''' crops a data set'''
    xr_lon = (lon >= limits[2]) & (lon <= limits[3])
    xr_lat = (lat >= limits[0]) & (lat <= limits[1])
    crop_xcm = xrin.where(xr_lon & xr_lat, drop=True)
    return crop_xcm

def getclosest_ij(lats,lons,latpt,lonpt):
    ''' find squared distance of every point on grid '''
    dist_sq = (lats-latpt)**2 + (lons-lonpt)**2
    minindex_flattened = dist_sq.argmin()
    return np.unravel_index(minindex_flattened, lats.shape)

def map2pfaf (anom, lats, lons, dlat, dlon, ulat, ulon, carea, levels):
    ''' maps variable to pfaf channel '''
    carr = []
    xx1 = []
    xx2 = []
    yy1 = []
    yy2 = []
    cua = []

    for i, value in enumerate(dlon):
        iy_ = min(range(len(lats)), key=lambda j: abs(lats[j] - ulat[i]))
        ix_ = min(range(len(lons)), key=lambda j: abs(lons[j] - ulon[i]))
        au_ = anom[iy_,ix_]
        iy_ = min(range(len(lats)), key=lambda j: abs(lats[j] - dlat[i]))
        ix_ = min(range(len(lons)), key=lambda j: abs(lons[j] - dlon[i]))
        ad_ = anom[iy_,ix_]
        if ad_ == -9999. or au_ == -9999.:
            continue

        if au_ and ad_:
            this_val = 0.5*(au_+ad_)
            this_val = min (this_val, max(levels) - 0.01)
            this_val = max (this_val, min(levels) + 0.01)
            k = min(range(len(levels)), key=lambda i: abs(levels[i]-this_val))
            carr.append(k)
            xx1.append(dlon[i])
            xx2.append(ulon[i])
            yy1.append(dlat[i])
            yy2.append(ulat[i])
            cua.append(carea[i])

    return carr, xx1, xx2, yy1, yy2, cua

class CachedTiler(object):
    ''' for Google tiles background '''
    def __init__(self, tiler, google_path):
        self.tiler = tiler
        self.google_path = google_path

    def __getattr__(self, name):
        attr = getattr(self.tiler, name, None)
        if isinstance(attr, types.MethodType):
            attr = types.MethodType(attr.__func__, self)
        return attr

    def get_image(self, tile):
        ''' down load image if not available '''
        tileset_name = '{}'.format(self.tiler.__class__.__name__.lower())
        cache_dir = os.path.expanduser(os.path.join(self.google_path, 'image_tiles', tileset_name))
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        tile_fname = os.path.join(cache_dir, '_'.join(str(v) for v in tile) + '.png')
        if not os.path.exists(tile_fname):
            response = requests.get(self._image_url(tile),
                                    stream=True, timeout=300)

            with open(tile_fname, "wb") as fh_:
                for chunk in response:
                    fh_.write(chunk)
        with open(tile_fname, 'rb') as fh_:
            img = PIL.Image.open(fh_)
            img = img.convert(self.desired_tile_form)
        return img, self.tileextent(tile), 'lower'

def contours (_x, _y, nrows, ncols, var, color_palette, titles, domain, figure, \
              under_over, min_val=None, max_val=None, fscale=None, levels=None, \
              stitle=None, clabel=None, cartopy_datadir=None):
    ''' plot contour maps'''
    cartopy_dir(cartopy_datadir)
    if fscale is None:
        fscale = FONT_SCALE

    if color_palette == 'RdYlGn':
        cmap = mpl.cm.get_cmap(color_palette).copy()
    else:
        style_color = load_table(color_palette)
        color_arr = []
        for color in style_color:
            rgb = [float(value) / 255 for value in color]
            color_arr.append(rgb)
        cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)

    if levels is None:
        if min_val is None:
            min_val = np.floor(np.nanmin(var [var>0.]))
            max_val = np.ceil(np.nanmax(var))
        levels = np.linspace(min_val, max_val, len(style_color))
    cmap.set_under(under_over[0])
    cmap.set_over(under_over[1])

    nplots = len(titles)
    fig = plt.figure(figsize= figure_size(FIGWIDTH, domain, nrows, ncols))
    gs_ = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)
    cax = fig.add_axes(cbar_axes)

    # plot maps
    for count_plot in range(nplots):
        ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
        cs_ = plt.pcolormesh(_x, _y, var[count_plot,],
                             norm=colors.BoundaryNorm(levels,ncolors=cmap.N, clip=False),
                             cmap=cmap,zorder=3, alpha=0.8)
        gl_ = ax_.gridlines(draw_labels=True)
        gl_.top_labels = False
        gl_.bottom_labels = False
        gl_.left_labels = False
        gl_.right_labels = False

        plt.title(titles[count_plot], fontsize=fscale*FONT_SIZE2)

        if np.mod (count_plot, ncols) == 0:
            gl_.left_labels = True

            if stitle is not None:
                plt.text(-0.15, 0.5, stitle, verticalalignment='center',
                         horizontalalignment='center', transform=ax_.transAxes,
                         color='Blue', fontsize=fscale*FONT_SIZE1, rotation=90)
        if (nplots - count_plot -1) < ncols:
            gl_.bottom_labels = True
        ax_.coastlines(alpha=0.1, zorder=3)
        ax_.add_feature(OCEAN, linewidth=0.2, zorder=3 )
        ax_.add_feature(COASTLINES, edgecolor='black', alpha=1, zorder=3)
        ax_.add_feature(BODR, linestyle='--', edgecolor='k', alpha=1, zorder=3)
        ax_.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')
        if (domain[3] - domain[2]) < 180.:
            ax_.add_feature(cfeature.STATES,  linestyle=':',linewidth=0.9,
                            edgecolor='black', facecolor='none')
        if under_over[0] == "white" and under_over[1] == "white":
            cbar = fig.colorbar(cs_, cax=cax, orientation='horizontal', ticks=levels)
        else:
            cbar = fig.colorbar(cs_, cax=cax, orientation='horizontal', ticks=levels,extend=EXTEND)
        cbar.ax.tick_params(labelsize=fscale*20, labelrotation=90)
        if clabel is not None:
            cbar.set_label(clabel, fontsize=fscale*30)
        plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')
    plt.close()

def google_map(_x, _y, nrows, ncols, var, color_palette, titles, domain, figure, \
               under_over, dlat, dlon, ulat, ulon, carea, google_path, min_val=None, \
               max_val=None, fscale=None, levels=None, \
               stitle=None, clabel=None, cartopy_datadir=None):
    ''' plots streams using google map as the background image'''
    cartopy_dir(cartopy_datadir)
    if fscale is None:
        fscale = FONT_SCALE

    style_color = load_table(color_palette)
    color_arr = []
    for color in style_color:
        rgb = [float(value) / 255 for value in color]
        color_arr.append(rgb)
    cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
    cmap.set_under(under_over[0])
    cmap.set_over(under_over[1])

    # Initialize our Google Maps tiles
    actual_tiler = cimgt.GoogleTiles()
    imagery = CachedTiler(actual_tiler, google_path)

    nplots = len(titles)
    fig = plt.figure(figsize= figure_size(FIGWIDTH, domain, nrows, ncols))
    gs_ = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)
    cax = fig.add_axes(cbar_axes)

    # plot maps

    for count_plot in range(nplots):
        ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
        ax_.set_extent([domain[3],domain[2],domain[0], domain[1]], crs=ccrs.Geodetic())
        ax_.add_image(imagery, 9)
        manom = var[count_plot,]
        carr, xx1, xx2, yy1, yy2, cua = map2pfaf(manom, _y, _x, dlat, dlon,
                                                 ulat, ulon, carea, levels)
        maxa = max(cua)
        mina = min(cua)
        pfaf_cnt = 0

        for this_col in carr:
            #if cua[pfaf_cnt] >= 10000.:
            this_col = max(this_col, 0)
            if this_col >= len(style_color):
                this_col = len(style_color) -1
            rgb = np.array(style_color[:][this_col])/255.
            track = sgeom.LineString(zip([xx1[pfaf_cnt],xx2[pfaf_cnt]],
                                         [yy1[pfaf_cnt], yy2[pfaf_cnt]]))
            lw_ = 1.5*(cua[pfaf_cnt] - mina)/(maxa - mina) + 1.
            ax_.add_geometries([track], ccrs.PlateCarree(),facecolor='none',
                               edgecolor=rgb, linewidth=lw_)
            pfaf_cnt += 1
        gl_ = ax_.gridlines(draw_labels=True)
        gl_.top_labels = False
        gl_.bottom_labels = False
        gl_.left_labels = False
        gl_.right_labels = False

        plt.title(titles[count_plot], fontsize=fscale*FONT_SIZE2)
        if np.mod (count_plot, ncols) == 0:
            gl_.left_labels = True

            if stitle is not None:
                plt.text(-0.15, 0.5, stitle, verticalalignment='center',
                         horizontalalignment='center', transform=ax_.transAxes,
                         color='Blue', fontsize=fscale*FONT_SIZE1, rotation=90)
        if (nplots - count_plot -1) < ncols:
            gl_.bottom_labels = True
        ax_.coastlines()
        ax_.add_feature(cfeature.BORDERS)
        ax_.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')

        cbar = fig.colorbar(plt.cm.ScalarMappable
                            (norm=colors.BoundaryNorm(levels, ncolors=cmap.N,clip=False),
                             cmap=cmap),
                            cax=cax, orientation='horizontal', ticks=levels,extend=EXTEND)
        cbar.ax.tick_params(labelsize=fscale*20)
        if clabel is not None:
            cbar.set_label(clabel, fontsize=fscale*30)
        plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')
    plt.close()

def stations(nrows, ncols, var, color_palette, titles, domain, figure,
             min_val=None, max_val=None, fscale=None, levels=None):
    ''' plots stattion on a map'''
    if min_val is None:
        min_val = np.floor(np.nanmin(var [var>0.]))
        max_val = np.ceil(np.nanmax(var))
    if fscale is None:
        fscale = FONT_SCALE

    style_color = load_table(color_palette)
    if levels is None:
        levels = np.linspace(min_val, max_val, len(style_color))

    color_arr = []
    for color in style_color:
        rgb = [float(value) / 255 for value in color]
        color_arr.append(rgb)
    cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
    cmap.set_under(COL_UNDER)
    cmap.set_over(COL_HIGHER)

    nplots = len(titles)
    fig = plt.figure(figsize= figure_size(FIGWIDTH, domain, nrows, ncols))
    gs_ = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)
    cax = fig.add_axes(cbar_axes)

    # plot maps
    count_plot = 0
    for key in titles:
        ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
        ax_.set_extent([domain[3],domain[2],domain[0],domain[1]], crs=ccrs.Geodetic())

        for sta in np.arange(0,len(var)):
            rmse = var[sta]
            _y, _x = rmse['loc'][0], rmse['loc'][1]
            this_val = rmse[key]
            this_val = min (this_val, max(levels) - 0.5)
            this_val = max (this_val, min(levels) + 0.5)
            this_col = min(range(len(levels)), key=lambda i: abs(levels[i]-this_val))
            rgb = np.array(style_color[:][this_col])/255.
            proj = ccrs.Orthographic(central_longitude=_x, central_latitude=_y)
            r_ortho = compute_radius(proj, PLOT_CIRCLE_RADIUS, _x,_y)
            ax_.add_patch(patches.Circle(xy=[_x,_y], radius=r_ortho, color=rgb, alpha=0.3,
                                         transform=proj, zorder=30))

        gl_ = ax_.gridlines(draw_labels=True)
        gl_.top_labels = False
        gl_.bottom_labels = False
        gl_.left_labels = False
        gl_.right_labels = False

        plt.title(titles[count_plot], fontsize=fscale*FONT_SIZE2)

        if np.mod (count_plot, ncols) == 0:
            gl_.left_labels = True

        if (nplots - count_plot -1) < ncols:
            gl_.bottom_labels = True

        ax_.coastlines()
        ax_.add_feature(cfeature.BORDERS)
        ax_.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')

        cbar = fig.colorbar(
            plt.cm.ScalarMappable(norm=colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False),
                                  cmap=cmap),cax=cax, orientation='horizontal', ticks=levels)
        cbar.ax.tick_params(labelsize=10, rotation = 90)
        plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')
        count_plot += 1
    plt.close()
