#!/usr/bin/env python

'''
This script plots GHI_S2S, IMERG, CHIRPS, and USAF precipitation anomalies on their native grids for a given month.

1) GHI_S2S precip anomalies are read from s2smetric outputs {yyyy}{mm}/ metrics_cf/AFRICOM/NOAHMP/{model}_Precip_ANOM_init_monthly_{mm}_{yyyy}.nc files

2) CHIRPS:
/discover/nobackup/projects/lis/MET_FORCING/CHIRPSv2/daily_p05/chirps-v2.0.{yyyy}.days_p05.nc CHIRPS daily precipitatoin data were used to 
derive monthly climatology for the period 2008-2020 while monthly CHIRPS precipitation is computed using daily precip data in below annual files.  
/discover/nobackup/projects/lis/MET_FORCING/CHIRPSv2/daily_p05/prelim/chirps-v2.0.{yyyy}.days_p05.nc  
Variable Name, Dimensions and units:
 precip(time, latitude, longitude), "mm/day"

3) IMERG:
/discover/nobackup/projects/lis/MET_FORCING/IMERG/Monthly_V06B/{yyyy}/3B-MO.MS.MRG.3IMERG.{yyyy}{mm}01-S000000-E235959.{mm}.V06B.HDF5 monthly data
were used to derive monthly climatology for the period 2008-2020 while /discover/nobackup/projects/lis/MET_FORCING/IMERG/Late_V06B/{yyyy}{mm}/ data are
used to compute precip for the month in question.

4) USAF raw 
/discover/nobackup/projects/usaf_lis/smahanam/s2splots/preproc_usaf_precip.csh processed
usaf_lis73rc8_10km precipitation files in /discover/nobackup/projects/usaf_lis/MET_FORCING/usaf_lis73rc8_10km/ and created the 13-year monthly climatology 
file: /discover/nobackup/projects/usaf_lis/smahanam/observed_precip/PS.AFWA_SC.U_DI.C_DC.ANLYS_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA_monthly_clim_2008-2020.nc
Monthly precip for the month in question is computed using precip data in:
/discover/nobackup/projects/usaf_lis/MET_FORCING/usaf_lis73rc8_10km/${yyyy}${mm}{dd}/
PS.AFWA_SC.U_DI.C_DC.ANLYS_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA.03-HR-SUM_DD.{yyyy}{mm}{dd}_DT.0000_DF.GR1
Variable Names, Dimensions and units:
Clim file: APCP(time, lat, lon), lon, lat, "kg m-2" -> mm/3hrs
monthly file: var61 (time, lat, lon), lat, lon

'''

import numpy as np
import xarray as xr
import netCDF4 as nc4
import os
import glob
import sys
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from load_colors import load_table
import yaml
import argparse

obs_precip_path = '/discover/nobackup/projects/usaf_lis/smahanam/observed_precip/'
figure_template = '{}/OBS_precip_anom_{}.png'

clim_file = {
    'IMERG': '3B-MO.MS.MRG.3IMERG.YYYYMM01-S000000-E235959.01.V06B_monthly_clim_2013-2021.nc4',
    'CHIRPS': 'chirps-v2.0.days_p05_monthly_clim_2013-2021.nc4',
    'USAF': 'PS.AFWA_SC.U_DI.C_DC.ANLYS_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA_monthly_clim_2013-2021.nc'
}
metforc_dir = {
    'IMERG': '/discover/nobackup/projects/lis/MET_FORCING/IMERG/Late_V06B/',
    'CHIRPS': '/discover/nobackup/projects/lis/MET_FORCING/CHIRPSv2/daily_p05/prelim/',
    'USAF': obs_precip_path + 'USAF_monthly/'
}

figwidth = 25
nrows = 2
ncols = 5
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
units = 'in/mon'
levels = [-10, -6, -4, -2, -1, -0.5, -0.25,0.25, 0.5, 1., 2., 4., 6., 10.]
#levels = [-10, -6, -4, -2, -1, -0.25, 0.25, 1., 2., 4., 6., 10.]
cmap, col_under, col_higher, extend = plt.cm.RdYlGn, 'black', '#B404AE', 'both'
#cmap.set_under(col_under)
#cmap.set_over(col_higher)
style_color = load_table('14WPR')
color_arr = []
for color in style_color:
    rgb = [float(value) / 255 for value in color]
    color_arr.append(rgb)
norm = mpl.colors.BoundaryNorm(levels, ncolors=256)
cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
cmap.set_under('gray')
cmap.set_over('blue')

mpl.style.use('bmh')
cbar_axes = [0.2, 0.04, 0.6, 0.03]
        
class Climatology (object):

    def chirps_clim ():
        chirp_template = '/discover/nobackup/projects/lis/MET_FORCING/CHIRPSv2/daily_p05/chirps-v2.0.{:04d}.days_p05.nc'
        prcp_a = []
        for year in range (2008,2021):
            # read annual files
            print (chirp_template.format(year))
            chirp = xr.open_dataset(chirp_template.format(year))
            # construct the list of annual datasets
            prcp_a.append (chirp.precip.groupby('time.month').mean(dim='time'))

        prcp_all = xr.concat(prcp_a, dim='month')
        clim = prcp_all.groupby('month').mean(dim='month')
        clim.to_netcdf(obs_precip_path + clim_file("CHIRPS"), format="NETCDF4", engine="netcdf4")
        return

    def imerg_clim ():
        imerg_template = '/discover/nobackup/projects/lis/MET_FORCING/IMERG/Monthly_V06B/{:04d}/'
        prcp_a = []
        for year in range (2008,2021):
            this_year = []
            files = sorted(os.listdir (imerg_template.format(year)))
            for file in files:
                # read monthly files
                print (file)
                ncf = nc4.Dataset(imerg_template.format(year)+ file, format='NETCDF4', diskless=True, persist=False)
                nch = ncf.groups.get('Grid')
                imerg = xr.open_dataset(xr.backends.NetCDF4DataStore(nch))
                # construct the list for the year in question
                this_year.append(imerg.precipitation)
            # construct the list of annual datasets
            prcp_a.append (xr.concat(this_year, dim='time'))
            
        prcp_all =  xr.concat(prcp_a, dim='time')
        clim = prcp_all.groupby('time.month').mean(dim='time')
        clim.to_netcdf(obs_precip_path + clim_file("IMERG"), format="NETCDF4", engine="netcdf4")
        return

class PlotMaps (object):

    def crop (limits,lon, lat, xrin):
        xr_lon = (lon >= limits[0]) & (lon <= limits[1])
        xr_lat = (lat >= limits[2]) & (lat <= limits[3])
        crop_xcm = xrin.where(xr_lon & xr_lat, drop=True)
        return crop_xcm
    
    def subplot (gs_,count_plot,ctitle,limits,plotdir,region,anom,longitude,latitude,figure):
        add_land,add_rivers,resol
        ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())        
        ax_.set_extent(limits, ccrs.PlateCarree())
        cs_ = plt.contourf(longitude, latitude, anom*30./25.4,levels, cmap=cmap,
                           extend=extend, transform=ccrs.PlateCarree())
        gl_ = ax_.gridlines(draw_labels=True)
        gl_.xlabels_top = False
        gl_.xlabels_bottom = False
        gl_.ylabels_right = False
        gl_.ylabels_left = False
        plt.title(ctitle, fontsize=25)
        if count_plot + 1> ncols:
            gl_.xlabels_bottom = True
        if count_plot % ncols == 0:
            gl_.ylabels_left = True
        ax_.coastlines()
        ax_.add_feature(cfeature.BORDERS)
        ax_.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')
        cax = fig.add_axes(cbar_axes)
        cbar = fig.colorbar(cs_, cax=cax, orientation='horizontal', ticks=levels)
        cbar.set_label('Anomaly (in/mon)', fontsize=30)
        cbar.ax.tick_params(labelsize=20)
        print(ctitle)
        print(figure)
        plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')
        
if __name__ == '__main__':

    ''' 
    NOTE:
    To compute CHIRPS and IMERG climatology on their native grids, run below 2 functions:
    Climatology.chirps_clim()
    Climatology.imerg_clim()
    '''
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
    
    proc_month = date(fcst_year, fcst_mon, 1) - relativedelta(months=1)
    yyyy = str(proc_month.year)
    mm = str(proc_month.month).zfill(2)
    year = proc_month.year
    month = proc_month.month
    # load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)
    
    plotdir_template = cwd + '/plots/{:04d}{:02d}/'
    plotdir = plotdir_template.format(fcst_year,fcst_mon)+  config["EXP"]["domain"] + '/' + config["EXP"]["lsmdir"] + '/'
    s2smetric_temp  = config["SETUP"]["E2ESDIR"] + '/s2smetric/output/{:04d}{:02d}/metrics_cf/' + \
        config["EXP"]["domain"] + '/'  + config["EXP"]["lsmdir"] + '/{}_Precip_ANOM_init_monthly_{:02d}_{:04d}.nc'
    get_boundary = {
        'EA': (22, 55, -12, 23),
        'WA': (-19, 26, -5, 25),
        'SA': (8, 52, -37, 6),
        'SA1':(24, 33, -31, -24),
        'FAME':(-20, 55, -40, 40)
    }
    
    # compute yyyymm monthly precip anomaly
    # -------------------------------------

    # generate USAF monthly precip file
    # ----------------------------------
    command = './preproc_usaf_precip.csh ' + yyyy + ' ' + mm
    #res = os.system(command)
    #assert res == 0,"USAF precipitation not ready."
    usaf_monthly = nc4.Dataset(metforc_dir.get('USAF') + yyyy + mm + '.nc')

    # (2) CHIRPS
    chirp = xr.open_dataset(metforc_dir.get('CHIRPS') + 'chirps-v2.0.' + yyyy + '.days_p05.nc')
     
    # Plotting anomalies
    # ------------------

    region = 'FAME'
    figure = figure_template.format(plotdir, yyyy+mm)    
    count_plot = 0

    
    fig = plt.figure(figsize=(figwidth,
                              figwidth*(nrows*(get_boundary.get(region)[3]
                                               - get_boundary.get(region)[2]))
                              /(ncols*(get_boundary.get(region)[1]
                                       - get_boundary.get(region)[0]))))
    gs_ = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)

    # (1) CHIRPS
    # ----------
    cfile = obs_precip_path + clim_file.get('CHIRPS')
    xc = xr.open_dataset(cfile)
    xcm = xc.isel(month=[month-1])
    crop_xcm = PlotMaps.crop (get_boundary.get(region),xcm.longitude, xcm.latitude, xcm)
    
    mfile = metforc_dir.get('CHIRPS') + 'chirps-v2.0.' + yyyy + '.days_p05.nc'    
    xm = xr.open_dataset(mfile)
    xmy = xm.precip.groupby('time.month').mean(dim='time')
    xmm = xmy.isel(month=[month-1])
    crop_xmm = PlotMaps.crop (get_boundary.get(region),xmm.longitude, xmm.latitude, xmm)
    anom = crop_xmm.values -  crop_xcm.precip.values
    PlotMaps.subplot(gs_,count_plot,'CHIRPS: '+ yyyy + mm, get_boundary.get(region),plotdir, region, anom[0,], crop_xcm.longitude, crop_xcm.latitude,figure) 
    count_plot += 1
    
    # (2) USAF raw
    cfile = obs_precip_path + clim_file.get('USAF')
    xc = xr.open_dataset(cfile)
    xcm = xc.isel(time=[month-1])
    crop_xcm = PlotMaps.crop (get_boundary.get(region),xcm.lon, xcm.lat, xcm)

    mfile = metforc_dir.get('USAF') + yyyy + mm + '.nc'
    xmm = xr.open_dataset(mfile)
    crop_xmm = PlotMaps.crop (get_boundary.get(region),xmm.lon, xmm.lat, xmm)
    anom = 8.*(crop_xmm.var61.values -  crop_xcm.APCP.values)  # coverted from mm/3hr to mm/day
    PlotMaps.subplot(gs_,count_plot,'USAF: '+ yyyy + mm, get_boundary.get(region),plotdir, region, anom[0,], crop_xcm.lon, crop_xcm.lat,figure) 
    count_plot += 1
    
    # (3) IMERG
    cfile = obs_precip_path + clim_file.get('IMERG')
    xc = xr.open_dataset(cfile)
    xcm = xc.isel(month=[month-1])
    crop_xcm = PlotMaps.crop (get_boundary.get(region),xcm.lon, xcm.lat, xcm)

    mdir = metforc_dir.get('IMERG') + yyyy + mm + '/' 
    mfiles = sorted(os.listdir (mdir))
    this_month = []
    this_day = []
    file_count = 0

    # generate monthly mean IMERG file
    for file in mfiles:
        # read 30min files
        ncf = nc4.Dataset(mdir + file, format='NETCDF4', diskless=True, persist=False)
        nch = ncf.groups.get('Grid')
        imerg = xr.open_dataset(xr.backends.NetCDF4DataStore(nch))    
        this_day.append(imerg.precipitationCal)
        file_count += 1
        if file_count == 48:
            file_count = 0
            xr_tmp = xr.concat(this_day, dim='time')
            this_month.append(xr_tmp.groupby('time.day').mean(dim='time'))
            this_day = 0
            this_day = []
    prcp_m = xr.concat(this_month, dim='day')
    xtmf = np.mean(prcp_m, axis =0)
    xtmf.to_netcdf(obs_precip_path + 'IMERG_monthly/' + yyyy+ mm+ '.nc4', format="NETCDF4", engine="netcdf4")
    prcp_m =0
    xtmf = 0
    
    xmm = xr.open_dataset(obs_precip_path + 'IMERG_monthly/' + yyyy+ mm+ '.nc4')
    crop_xmm = PlotMaps.crop (get_boundary.get(region),xmm.lon, xmm.lat, xmm)
    anom = 24.*(crop_xmm.precipitationCal.values -  crop_xcm.precipitation.values)  # coverted from mm/hr to mm/day
    PlotMaps.subplot(gs_,count_plot,'IMERG: '+ yyyy + mm, get_boundary.get(region),plotdir, region, np.transpose(anom[0,]), crop_xcm.lon, crop_xcm.lat,figure) 
    count_plot += 1

    # (4) GHI S2S
    infile = s2smetric_temp.format(year,month, '*', month, year)
    anom = xr.open_mfdataset(infile, concat_dim='ens')
    median_anom = np.median(anom.anom.values, axis=0)
    PlotMaps.subplot(gs_,count_plot,'GHI S2S: '+ yyyy + mm, get_boundary.get(region),plotdir, region, median_anom[0, ], anom.longitude, anom.latitude,figure) 
    count_plot += 1
    
    # (5) CCM4
    infile = s2smetric_temp.format(year,month, 'CCM4', month, year)
    anom = xr.open_dataset(infile)
    ens_median = np.median(anom.anom.values, axis=0)
    PlotMaps.subplot(gs_,count_plot,'CCM4: '+ yyyy + mm, get_boundary.get(region),plotdir, region, ens_median[0, ], anom.longitude, anom.latitude,figure) 
    count_plot += 1
    
    # (6) CCSM4
    infile = s2smetric_temp.format(year,month, 'CCSM4', month, year)
    anom = xr.open_dataset(infile)
    ens_median = np.median(anom.anom.values, axis=0)
    PlotMaps.subplot(gs_,count_plot,'CCSM4: '+ yyyy + mm, get_boundary.get(region),plotdir, region, ens_median[0, ], anom.longitude, anom.latitude,figure)    
    count_plot += 1
    
    # (7) CFSv2
    infile = s2smetric_temp.format(year,month, 'CFSv2', month, year)
    anom = xr.open_dataset(infile)
    ens_median = np.median(anom.anom.values, axis=0)
    PlotMaps.subplot(gs_,count_plot,'CFSv2: '+ yyyy + mm, get_boundary.get(region),plotdir, region, ens_median[0, ], anom.longitude, anom.latitude,figure)    
    count_plot += 1
    
    # (8) GEOSv2
    infile = s2smetric_temp.format(year,month, 'GEOSv2', month, year)
    anom = xr.open_dataset(infile)
    ens_median = np.median(anom.anom.values, axis=0)
    PlotMaps.subplot(gs_,count_plot,'GEOSv2: '+ yyyy + mm, get_boundary.get(region),plotdir, region, ens_median[0, ], anom.longitude, anom.latitude,figure)        
    count_plot += 1
    
    # (9) GFDL
    infile = s2smetric_temp.format(year,month, 'GFDL', month, year)
    anom = xr.open_dataset(infile)
    ens_median = np.median(anom.anom.values, axis=0)
    PlotMaps.subplot(gs_,count_plot,'GFDL: '+ yyyy + mm, get_boundary.get(region),plotdir, region, ens_median[0, ], anom.longitude, anom.latitude,figure)        
    count_plot += 1
    
    #(10) GNEMO5
    infile = s2smetric_temp.format(year, month, 'GNEMO5', month, year)
    anom = xr.open_dataset(infile)
    ens_median = np.median(anom.anom.values, axis=0)
    PlotMaps.subplot(gs_,count_plot,'GNEMO5: '+ yyyy + mm, get_boundary.get(region),plotdir, region, ens_median[0, ], anom.longitude, anom.latitude,figure)        
    count_plot += 1

#/gpfsm/dnb06/projects/p159/jjacob/AFWA_EVAL/Python/IMERG/IMERG_reg_all_daily_ops_SL12.py
