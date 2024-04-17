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
This script derives monthly climatology of the METFORCE used in the LIS-DA step.
'''

import os
import sys
import glob
import argparse
from datetime import date
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib import colors
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from dateutil.relativedelta import relativedelta
import xarray as xr
import xesmf as xe
import numpy as np
import yaml
import eccodes
# pylint: disable=invalid-name, consider-using-f-string, import-outside-toplevel, redefined-outer-name

USAF_COLORS = True
FONT_SIZE1 = 16
FONT_SIZE2 = 20
FONT_SCALE = 1.

# Plotting Parameters
ADD_LAND = True
ADD_RIVERS = True
RESOL = '50m'  # use data at this scale
FIGWIDTH = 25
EXTEND = 'both'

mpl.use('pdf')
mpl.style.use('bmh')

OUT_VARS = ['AirT', 'Precip']

usaf_vars = {
    'AirT': 'TMP',
    'Precip': 'APCP'
    }

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
def compute_median (anom, lead):
    ''' compute median anomaly of the specified lad '''
    da_slice=[]
    for da in anom:
        da_slice.append(da.isel(time=lead))
    da_conc = xr.concat(da_slice, dim = 'ens')
    return da_conc.median(dim = 'ens')

def copy_1d_to_2d(data1d, Ni, Nj):
    ''' vector to 2D '''
    data2d = np.ndarray((Nj, Ni), data1d.dtype)
    for j in range(0, Nj):
        start = j*Ni
        end = (j+1)*Ni
        data2d[j,:] = data1d[start:end]
    return data2d

def grib2_xr(gfile, args):
    ''' read GR1 file, regrid to the LIS grid and return a xarray dataset '''
    ftn = open(gfile, 'rb')
    nmsgs = eccodes.codes_count_in_file(ftn)

    for i in range(0, nmsgs):
        gid = eccodes.codes_grib_new_from_file(ftn)

        Ni = eccodes.codes_get(gid, 'Ni')
        Nj = eccodes.codes_get(gid, 'Nj')
        yearOfCentury = eccodes.codes_get(gid, 'yearOfCentury')
        month = eccodes.codes_get(gid, 'month')
        day = eccodes.codes_get(gid, 'day')
        hour = eccodes.codes_get(gid, 'hour')
        minute = eccodes.codes_get(gid, 'minute')
        indicatorOfParameter = \
            eccodes.codes_get(gid, 'indicatorOfParameter')
        data1d = eccodes.codes_get_array(gid, 'values')
        data1d[data1d == 9999.] = np.nan

        np.ma.set_fill_value(data1d, -9999.)
        if indicatorOfParameter == 11:
            t_2m = copy_1d_to_2d(data1d, Ni, Nj)
        elif indicatorOfParameter == 61:
            pcp3hr = copy_1d_to_2d(data1d, Ni, Nj)
        del data1d
        eccodes.codes_release(gid)
    ftn.close()

    lats = np.ndarray(Nj, np.dtype('f4'))
    lons = np.ndarray(Ni, np.dtype('f4'))
    for i in range(0, Ni):
        lons[i] = -179.9297 + (i*0.140625)
    for j in range(0, Nj):
        lats[j] = -89.95312 + (j*0.09375)

    dt = xr.DataArray(t_2m, dims=('lat', 'lon'), coords={'lat': lats, 'lon': lons})
    dp = xr.DataArray(pcp3hr, dims=('lat', 'lon'), coords={'lat': lats, 'lon': lons})
    ds_in = xr.Dataset(
        {
            "lat": (["lat"], lats),
            "lon": (["lon"], lons),
        })
    ds_in['APCP'] = dp
    ds_in['TMP'] = dt

    del t_2m
    del pcp3hr

    # Now regrid to the LIS grid
    if args["run_regrid"] is True:
        # resample to the S2S grid
        args["run_regrid"] = False
        ds_in["slice"] = ds_in["TMP"]
        args["regridder"] = xe.Regridder(ds_in, args["ds_out"], "bilinear", periodic=True)

    ds_out2 = args["regridder"](ds_in)
    ds_out = ds_out2.drop_vars('slice', errors="ignore")

    return ds_out

def comp_monthly_nafpa (config_file, year, month, flabel):
    ''' The function regrid to output grid and compute monthly mean '''

    with open(config_file, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    from s2s_modules.bcsd_fcst.bcsd_library.bcsd_stats_functions import get_domain_info


    # Read all CFSv2 forcings
    met_path = cfg["SETUP"]["METFORC"] + '/' + flabel + '/'

    file_list = []
    file_list.extend(glob.glob(met_path + '{:04d}{:02d}*/*.GR1'.format(year, month)))

    # read in GR1 files and regrid to a xarray Dataset
    ds = []

    lats, lons = get_domain_info(config_file, coord=True)
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], lats),
            "lon": (["lon"], lons),
        }
    )
    args = {'run_regrid' : True,
            'ds_out': ds_out,
            'regridder': 0}

    for gfile in file_list:
        print(gfile)
        ds.append(grib2_xr(gfile, args))

    ds_xr =  xr.concat(ds, dim = 'time')
    ds_mean = ds_xr.mean(dim = 'time')
    return ds_mean

def write_clim(config_file, month, flabel):
    ''' The function regrid to output grid and compute climatological omean across clim years '''

    with open(config_file, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    from s2s_modules.bcsd_fcst.bcsd_library.bcsd_stats_functions import get_domain_info

    clim_years = cfg["POST"]["verification_clim"]
    cyear_beg = clim_years[0]
    cyear_end = clim_years[1]

    outdir = 'hindcast/bcsd_fcst/{}/Climatology_{:04d}-{:04d}/'.format(flabel, cyear_beg, cyear_end)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Read all CFSv2 forcings
    met_path = cfg["SETUP"]["METFORC"] + '/' + flabel + '/'

    file_list = []
    for year in range(cyear_beg, cyear_end + 1):
        file_list.extend(glob.glob(met_path + '{:04d}{:02d}*/*.GR1'.format(year, month)))

    outfile = outdir + ".".join(file_list[0].split("/")[-1].split(".")[:-3]) + '.{:02d}.nc4'.format(month)
    print (outfile)

    # read in GR1 files and regrid to a xarray Dataset
    ds = []

    lats, lons = get_domain_info(config_file, coord=True)
    ds_out = xr.Dataset(
        {
            "lat": (["lat"], lats),
            "lon": (["lon"], lons),
        }
    )
    args = {'run_regrid' : True,
            'ds_out': ds_out,
            'regridder': 0}

    for gfile in file_list:
        print(gfile)
        ds.append(grib2_xr(gfile, args))

    ds_xr =  xr.concat(ds, dim = 'time')
    ds_mean = ds_xr.mean(dim = 'time')
    ds_mean['TMP'].attrs['units'] = 'K'
    ds_mean['APCP'].attrs['units'] = 'mm/3hr'
    ds_mean.to_netcdf(outfile, format="NETCDF4",
            encoding = {'TMP':{"zlib":True, "complevel":6,
                               "shuffle":True, "missing_value":-9999.},
                        'APCP': {"zlib":True, "complevel":6,
                                 "shuffle":True, "missing_value":-9999.}})
    sys.exit()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--fmonth', required=True, help='month 1-12')
    parser.add_argument('-y', '--fyear', required=True, help='year')
    parser.add_argument('-c', '--config_file', required=True, help='config file')
    parser.add_argument('-l', '--lag_mon', required=True, help='lag month')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')

    args = parser.parse_args()
    fmonth = int(args.fmonth)
    fyear = int(args.fyear)
    lag = int(args.lag_mon)
    cwd = args.cwd

    flabel = 'usaf_lis75s2s_gfs2galwem'
    lead = lag - 1
    clim_month = fmonth - lag
    if clim_month <= 0.:
        clim_month = clim_month + 12

    ic_date = date(fyear, fmonth, 1) - relativedelta(months=lag)
    year = ic_date.year
    month = ic_date.month
    fdate = date(year, month, 1)
    vdate_beg = date(year, month, 1) + relativedelta(months=lead)
    vdate_end = vdate_beg + relativedelta(months=1) - relativedelta(days=1)
    fstr = fdate.strftime("%Y-%m-%d")
    vstr1 = vdate_beg.strftime("%d%b%Y")
    vstr2 = vdate_end.strftime("%d%b%Y")

    # (1) write NAFPA climatology
    # write_clim(args.config_file, month, flabel)

    with open(args.config_file, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)

    sys.path.append(cfg['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.s2splots import plot_utils
    plotdir = cwd + '/s2splots/{:04d}{:02d}/'.format(fyear,fmonth) + cfg["EXP"]["lsmdir"] + '/'
    cartopy.config['data_dir'] = cfg['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/'
    ndays = (date(year, month+1, 1) - date(year, month, 1)).days

    # (1) usaf_lis75s2s_gfs2galwem anomaly
    #climdir = cfg["SETUP"]["E2ESDIR"] + \
    #    '/hindcast/bcsd_fcst/{}/Climatology_{:04d}-{:04d}/'.format(flabel, cyear_beg, cyear_end)
    # nafpa_clim_xr = xr.open_dataset(glob.glob(climdir + '*{:02d}.nc4'.format(month))[0])

    nafpa_mon_xr =  comp_monthly_nafpa(args.config_file, year, month, flabel)
    #nafpa_mon_xr.to_netcdf('nafpa_05.nc4', format="NETCDF4",
    #                  encoding = {'TMP':{"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.},
    #                              'APCP': {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.}})
    #nafpa_mon_xr = xr.open_dataset('nafpa_{:02d}.nc4'.format(month))

    # forcecast anomaly
    end_date = date(year, month, 1) + relativedelta(months=int(cfg["EXP"]["lead_months"]))
    s2smdir = cwd + '/s2smetric/{:04d}{:02d}/metrics_cf/'.format(year,month) + \
        cfg["EXP"]["lsmdir"] + '/'
    domain = plot_utils.dicts('boundary', 'GLOBAL')

    levels_dict = {
        'AirT': [-12., -8., -6., -4., -2., -0.5, 0.5, 2., 4., 6., 8., 12.],
        'Precip': plot_utils.dicts('anom_levels','Precip_AF')}
    tables = {'AirT': '14WT2M',
              'Precip': '14WPR'
        }
    conv_factors = {
        'AirT': 9./5.,
        'Precip': ndays/25.4}  # convert to in/month
    clabels = {
        'AirT': 'Anomaly (' + plot_utils.dicts('units', 'Air_T_AF') + ')',
        'Precip': 'Anomaly (' + plot_utils.dicts('units', 'Precip_AF') + ')'}

    titles = ['Forecast:'+fstr + '            Valid: '+vstr1+'-'+vstr2, 'Observed VT: ' +vstr1+'-'+vstr2, 'Forecast:'+fstr + '            Valid: '+vstr1+'-'+vstr2]
    cbar_axes_horizontal = [0.15, 0.04, 0.7, 0.02]
    cbar_axes_vertical = [0.9, 0.37, 0.03, 0.5] # [left, bottom, width, height]

    for var in OUT_VARS:
        # read USAF-LIS7.3rc8_25km 30-year climatology
        if var == 'AirT':
            usaf_clim_xr =  xr.open_dataset(cfg["SETUP"]["E2ESDIR"] + \
                                            '/hindcast/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/T2M_obs_clim.nc')
        if var == 'Precip':
            usaf_clim_xr =  xr.open_dataset(cfg["SETUP"]["E2ESDIR"] + \
                                            '/hindcast/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/PRECTOT_obs_clim.nc')
        usaf_clim_var = np.mean(usaf_clim_xr['clim'].values[clim_month,:], axis=0)
        usaf_clim_xr.close()
        # var specifics
        load_table = tables.get(var)
        levels = levels_dict.get(var)
        under_over = plot_utils.dicts('lowhigh', load_table)
        convf = conv_factors.get(var)
        clabel = clabels.get(var)
        figure = plotdir + var + '_verification_F' + fdate.strftime("%Y%m%d") + '_V' + vdate_beg.strftime("%Y%m%d") + '-' + vdate_end.strftime("%Y%m%d") + '.png'

        if var == 'AirT':
            nafpa_anom = (np.array(nafpa_mon_xr[usaf_vars.get(var)].values) - usaf_clim_var)*convf # nafpa_clim_xr[usaf_vars.get(var)].values)*convf
        if var == 'Precip':
            nafpa_anom = (np.array(nafpa_mon_xr[usaf_vars.get(var)].values)*8. - usaf_clim_var*86400.)*convf
#            nafpa_anom = nafpa_anom * 8. # nafpa anoms are mm/[3hr] converted to per mm/[day]
        plot_arr = np.zeros([3, nafpa_anom.shape[0], nafpa_anom.shape[1]],dtype=float)
        plot_arr[1,:] = nafpa_anom
        tmp_arr = np.zeros([nafpa_anom.shape[0], nafpa_anom.shape[1]],dtype=float)
        # read NC files
        anoms = []
        print (s2smdir + 'PS.557WW_SC.U_DI.C_GP.LIS-S2S-*_GR.C0P25DEG_AR.GLOBAL_PA.S2SMETRICS_DD.{:04d}{:02d}01_FP.{:04d}{:02d}01-{}_DF.NC'.format(year, month, year, month, end_date.strftime("%Y%m%d")))
        for model in cfg["EXP"]["NMME_models"]:
            ncfile = s2smdir + 'PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.GLOBAL_PA.S2SMETRICS_DD.{:04d}{:02d}01_FP.{:04d}{:02d}01-{}_DF.NC'.format(model.upper(), year, month, year, month, end_date.strftime("%Y%m%d"))
            ncdata = xr.open_dataset(ncfile)
            anoms.append(ncdata[var + '_ANOM'])
            del ncdata

        nc_med = compute_median (anoms, lead)
        plot_arr[0,:] = nc_med.values *convf
        if var == 'Precip':
           plot_arr[0,:] = plot_arr[0,:]/ndays 

        # plotting
        style_color = plot_utils.load_table(load_table)
        color_arr = []
        for color in style_color:
            rgb = [float(value) / 255 for value in color]
            color_arr.append(rgb)
        cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
        cmap.set_under(under_over[0])
        cmap.set_over(under_over[1])
        nplots = len(titles)

        nrows = 3
        ncols = 1
        fscale= 1.
        fig = plt.figure(figsize= plot_utils.figure_size(FIGWIDTH, domain, nrows, ncols))
        gs_ = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)
        cax = fig.add_axes(cbar_axes_vertical)

        # plot maps
        for count_plot in range(nplots-1):

            ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
            cs_ = plt.pcolormesh(nc_med.longitude.values, nc_med.latitude.values, plot_arr[count_plot,],
                                 norm=colors.BoundaryNorm(levels,ncolors=cmap.N, clip=False),
                                 cmap=cmap,zorder=3, alpha=0.8)
            gl_ = ax_.gridlines(draw_labels=True)
            gl_.top_labels = False
            gl_.bottom_labels = False
            gl_.left_labels = False
            gl_.right_labels = False

            plt.title('Monthly '+ var + ' Anomaly : ' +titles[count_plot], fontsize=fscale*FONT_SIZE2)

            if np.mod (count_plot, ncols) == 0:
                gl_.left_labels = True

            if (nplots - count_plot -1) < ncols:
                gl_.bottom_labels = True
            ax_.coastlines(alpha=0.1, zorder=3)
            if var == 'Precip':
                ax_.add_feature(OCEAN, linewidth=0.2, zorder=3 )
                ax_.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')
            ax_.add_feature(COASTLINES, edgecolor='black', alpha=1, zorder=3)
            ax_.add_feature(BODR, linestyle='--', edgecolor='k', alpha=1, zorder=3)

            cbar = fig.colorbar(cs_, cax=cax, orientation='vertical', ticks=levels,extend=EXTEND)
            cbar.ax.tick_params(labelsize=fscale*20)
            if clabel is not None:
                cbar.set_label(clabel, fontsize=fscale*30)
            plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')

        # plot verification plot
        # Assigning values based on the conditions
        tmp_arr[np.logical_and(plot_arr[0,:] * plot_arr[1,:] < 0, plot_arr[0,:] != 0)] = 3
        tmp_arr[np.logical_and(plot_arr[0,:] < 0, plot_arr[1,:] < 0)] = 2
        tmp_arr[np.logical_and(plot_arr[0,:] > 0, plot_arr[1,:] > 0)] = 1

        if var == 'Precip':
            norms = [-1, 1.]
        else:
            norms = [-3, 3.]
        tmp_arr[np.logical_and(np.logical_and(plot_arr[0,:] >= norms[0], plot_arr[0,:] <= norms[1]),
                               np.logical_and(plot_arr[1,:] >= norms[0], plot_arr[1,:] <= norms[1]))] = 0
        plot_arr[2,:] = tmp_arr

        levels = [1, 2, 3, 4]
        style_color = plot_utils.load_table('VERIFY')
        color_arr = []
        for color in style_color:
            rgb = [float(value) / 255 for value in color]
            color_arr.append(rgb)
        cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
        cmap.set_under('white')
        cmap.set_over('white')
        count_plot = nplots-1
        cax2 = fig.add_axes(cbar_axes_horizontal)
        ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())
        cs_ = plt.pcolormesh(nc_med.longitude.values, nc_med.latitude.values, plot_arr[count_plot,],
                             norm=colors.BoundaryNorm(levels,ncolors=cmap.N, clip=False),
                             cmap=cmap,zorder=3)
        gl_ = ax_.gridlines(draw_labels=True)
        gl_.top_labels = False
        gl_.bottom_labels = False
        gl_.left_labels = False
        gl_.right_labels = False

        plt.title('Monthly '+var+ 'Anomaly Verification:' +titles[count_plot], fontsize=fscale*FONT_SIZE2)

        if np.mod (count_plot, ncols) == 0:
            gl_.left_labels = True

        if (nplots - count_plot -1) < ncols:
            gl_.bottom_labels = True
        ax_.coastlines(alpha=0.1, zorder=3)
        if var == 'Precip':
            ax_.add_feature(OCEAN, linewidth=0.2, zorder=3 )
            ax_.add_feature(cfeature.OCEAN, zorder=100, edgecolor='k')

        ax_.add_feature(COASTLINES, edgecolor='black', alpha=1, zorder=3)
        ax_.add_feature(BODR, linestyle='--', edgecolor='k', alpha=1, zorder=3)

        cbar = fig.colorbar(cs_, cax=cax2, orientation='horizontal', ticks=levels, extend=None)
        ticks = levels                     # Define the tick locations
        labels = ['Pos', 'Neg', 'Opp', ''] # Define the custom labels
        cbar.set_ticks(ticks)              # Set the tick locations
        cbar.set_ticklabels(labels)        # Assign the custom labels
        cbar.ax.tick_params(labelsize=fscale*20)

        plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')
        plt.close()
