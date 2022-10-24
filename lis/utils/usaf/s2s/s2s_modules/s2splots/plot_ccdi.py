import numpy as np
import datetime
import xarray as xr
from  scipy import stats
from scipy.interpolate import interp1d
from netCDF4 import Dataset
import sys
import os
import yaml
import argparse
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
from calendar import monthrange

vars = ['TotalPrecip_acc','Tair_f_max', 'Tair_f_min']
nbins = 41
bins = np.linspace(0, 1, nbins*1.)
HINDDIR = '/discover/nobackup/projects/usaf_lis/GHI_S2S/AFRICOM/hindcasts/'
outdir_template = '/discover/nobackup/projects/usaf_lis/smahanam/CCDI_analyses/processed/{}/{:02d}/percentile_values/'
daily_template = '{}/{:02d}/cf_{}_*{:02d}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_PA.LIS-S2S_DD.*{:02d}{:02d}_DT.0000_DF.NC'
NC = 320
NR = 320
var_name = {
    'R99d': 'TotalPrecip_acc',
    'R95d': 'TotalPrecip_acc',
    'R90d': 'TotalPrecip_acc',
    'TN10p': 'Tair_f_min',
    'TX90p': 'Tair_f_max',
    }

cutoff_name = {
    'R99d': 'PRCP_99P',
    'R95d': 'PRCP_95P',
    'R90d': 'PRCP_90P',
    'TN10p': 'TMIN_10P',
    'TX90p': 'TMAX_90P',
    }

get_boundary = {
    'EA': (22, 55, -12, 23),
    'WA': (-19, 26, -5, 25),
    'SA': (8, 52, -37, 6),
    'SA1':(24, 33, -31, -24),
    'FAME':(-20, 55, -40, 40)
    }

figwidth = 25
nrows = 2
ncols = 3
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
units = '%'
levels = np.linspace(0, 100., 21)
cmap, col_under, col_higher, extend = plt.cm.RdYlGn, 'black', '#B404AE', 'both'
style_color = load_table('L21')
color_arr = []
for color in style_color:
    rgb = [float(value) / 255 for value in color]
    color_arr.append(rgb)
norm = mpl.colors.BoundaryNorm(levels, ncolors=256)
cmap = colors.LinearSegmentedColormap.from_list('my_palette', color_arr, N=256)
cmap.set_under('white')
cmap.set_over('white')
mpl.style.use('bmh')
cbar_axes = [0.2, 0.04, 0.6, 0.03]
figure_template = '{}/{}_{}.png'
           
class getCDFs (object):
    
    def preproc(self,ds):
        ds = ds.isel(ensemble=0)
        return ds
    
    def gather_window (self,fcst_mon,date_process, model):
        # precipitation
        index_no = 0
        prcp = []
        ncdata1 = []
        # 7-day window since we count all ensemble members
        for i in range(-3,4,1):
            date = date_process + datetime.timedelta(days=i)
            daily_files = HINDDIR + daily_template.format(model, fcst_mon, model,date.month, model,date.month,date.day)
            print (daily_files)
            ncdata1.append(xr.open_mfdataset(daily_files, concat_dim='ensemble'))
            prcp.append(xr.Dataset())
            prcp[index_no]['TotalPrecip_acc']=(('ensemble','lat','lon'),ncdata1[index_no].TotalPrecip_acc.values)
            index_no += 1

        # temperation
        index_no = 0
        tair_max = []
        tair_min = []
        ncdata2 = []
        # 15-day window since we use only the first ensemble member
        for i in range(-7,8,1):
            date = date_process + datetime.timedelta(days=i)
            daily_files =  HINDDIR +'{}/{:02d}/cf_{}_*{:02d}/PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_PA.LIS-S2S_DD.*{:02d}{:02d}_DT.0000_DF.NC'.format(model, fcst_mon, model,date.month, model,date.month,date.day)
            ncdata2.append(xr.open_mfdataset(daily_files, concat_dim='ensemble',preprocess=self.preproc))
            tair_max.append(xr.Dataset())
            tair_max[index_no]['Tair_f_max']=(('ensemble','lat','lon'),ncdata2[index_no].Tair_f_max.values)
            tair_min.append(xr.Dataset())
            tair_min[index_no]['Tair_f_min']=(('ensemble','lat','lon'),ncdata2[index_no].Tair_f_min.values)
            index_no += 1
        
        return xr.concat(prcp, dim='ensemble'), xr.concat(tair_max, dim='ensemble'), xr.concat(tair_min, dim='ensemble')

    def get_cdf_func (self,ts_data):
        mx = -9999.
        mn = -9999.
        cdf = np.full((nbins), -9999.)
        if np.nanmax(ts_data, axis=None) > 0.:
            mx = np.nanmax(ts_data, axis=None)
            mn = np.nanmin(ts_data, axis=None)
            ts_norm = (ts_data - mn)/(mx - mn)
            res = stats.cumfreq(ts_norm, numbins=nbins)
            cdf = res.cumcount/len(ts_norm)
            
        return mx, mn, cdf

    def get_prcp_cutoff (self,ts_data):
        p90 = -9999.
        p95 = -9999.
        p99 = -9999.
        if np.nanmax(ts_data, axis=None) > 0.:
            mx = np.nanmax(ts_data, axis=None)
            mn = np.nanmin(ts_data, axis=None)
            ts_norm = (ts_data - mn)/(mx - mn)
            res = stats.cumfreq(ts_norm, numbins=nbins)
            cdf = res.cumcount/len(ts_norm)
            if cdf[0] > 0.9:
                cdf[0] = 0.9
            f = interp1d(cdf, bins )
            p90 = (mx - mn)*f(0.90) + mn
            p95 = (mx - mn)*f(0.95) + mn
            p99 = (mx - mn)*f(0.99) + mn

        return p90, p95, p99

    def get_temp_cutoff (self,ts_data):
        t10 = -9999.
        t90 = -9999.
        cdf = np.full((nbins), -9999.)
        if np.nanmax(ts_data, axis=None) > 0.:
            mx = np.nanmax(ts_data, axis=None)
            mn = np.nanmin(ts_data, axis=None)
            ts_norm = (ts_data - mn)/(mx - mn)
            res = stats.cumfreq(ts_norm, numbins=nbins)
            cdf = res.cumcount/len(ts_norm)
            f = interp1d(cdf, bins )
            t10 = (mx - mn)*f(0.1) + mn
            t90 = (mx - mn)*f(0.9) + mn

        return t10, t90
    
    def __init__  (self,fcst_mon, cwd, config, nmme):
        '''
        compute CDFs separatley for each month and NMME model
        '''
        for month in range (fcst_mon + 1, fcst_mon + 4):
            # Loop through the month
            start_date = datetime.datetime(2003,month,1,1,0,0)
            end_date = datetime.datetime(2003,month + 1,1,1,0,0)
            delta =  datetime.timedelta(days=1)
            
            while start_date < end_date:                
                print(start_date.strftime("%Y-%m-%d"))
                # Loop through models
                if not nmme:
                    # for nmme_model in  config['EXP']['NMME_models']:
                    nmme_models = config['EXP']['NMME_models']
                else:
                    nmme_models = [nmme]
                    
                for nmme_model in nmme_models:
                    HINDDIR = config["POST"]["hindcasts"] 
                    outdir = outdir_template.format(nmme_model,fcst_mon)
                    if not os.path.exists(outdir):
                        os.makedirs(outdir, exist_ok=True)
                    prcp, tair_max, tair_min = self.gather_window(fcst_mon,start_date, nmme_model)

                    lat_count = len(prcp.coords['lat'])
                    lon_count = len(prcp.coords['lon'])
                    pr_99  = np.full((lat_count, lon_count),-9999.)
                    pr_95  = np.full((lat_count, lon_count),-9999.)
                    pr_90  = np.full((lat_count, lon_count),-9999.)
                    tn_10  = np.full((lat_count, lon_count),-9999.)
                    tx_90  = np.full((lat_count, lon_count),-9999.)

                    for lat in range(lat_count):
                        for lon in range(lon_count):
                            pr_90[lat,lon], pr_95[lat,lon], pr_99[lat,lon] = self.get_prcp_cutoff (prcp.isel(lat=lat, lon=lon).TotalPrecip_acc.values)
                            tn_10[lat,lon], dummy = self.get_temp_cutoff (tair_min.isel(lat=lat, lon=lon).Tair_f_min.values)
                            dummy, tx_90[lat,lon] = self.get_temp_cutoff (tair_max.isel(lat=lat, lon=lon).Tair_f_max.values)
                            
#                    pr_max = np.full((lat_count, lon_count),-9999.)
#                    pr_min = np.full((lat_count, lon_count),-9999.)
#                    pr_cdf = np.full((lat_count,lon_count,nbins), -9999.)

#                    tx_max = np.full((lat_count, lon_count),-9999.)
#                    tx_min = np.full((lat_count, lon_count),-9999.)
#                    tx_cdf = np.full((lat_count,lon_count,nbins), -9999.)

#                    tn_max = np.full((lat_count, lon_count),-9999.)
#                    tn_min = np.full((lat_count, lon_count),-9999.)
#                    tn_cdf = np.full((lat_count,lon_count,nbins), -9999.)

#                    for lat in range(lat_count):
#                        for lon in range(lon_count):
#                            pr_max[lat,lon],pr_min[lat,lon],pr_cdf[lat,lon,:] = self.get_cdf_func (prcp.isel(lat=lat, lon=lon).TotalPrecip_acc.values)
#                            tx_max[lat,lon],tx_min[lat,lon],tx_cdf[lat,lon,:] = self.get_cdf_func (tair_max.isel(lat=lat, lon=lon).Tair_f_max.values)
#                            tn_max[lat,lon],tn_min[lat,lon],tn_cdf[lat,lon,:] = self.get_cdf_func (tair_min.isel(lat=lat, lon=lon).Tair_f_min.values)
#                    cdf_xr = xr.Dataset()
#                    cdf_xr['PRCP_CDF'] = (('lat','lon','bin'),pr_cdf)
#                    cdf_xr['PRCP_MaxValue'] = (('lat','lon'),pr_max)
#                    cdf_xr['PRCP_MinValue'] = (('lat','lon'),pr_min)

#                    cdf_xr['TMAX_CDF'] = (('lat','lon','bin'),tx_cdf)
#                    cdf_xr['TMAX_MaxValue'] = (('lat','lon'),tx_max)
#                    cdf_xr['TMAX_MinValue'] = (('lat','lon'),tx_min)

#                    cdf_xr['TMIN_CDF'] = (('lat','lon','bin'),tn_cdf)
#                    cdf_xr['TMIN_MaxValue'] = (('lat','lon'),tn_max)
#                    cdf_xr['TMIN_MinValue'] = (('lat','lon'),tn_min)
                    
#                    cdf_xr['BINS'] = (('bin'),bins)

                    cdf_xr = xr.Dataset()
                    cdf_xr['PRCP_90P'] = (('lat','lon'),pr_90)
                    cdf_xr['PRCP_95P'] = (('lat','lon'),pr_95)
                    cdf_xr['PRCP_99P'] = (('lat','lon'),pr_99)
                    cdf_xr['TMIN_10P'] = (('lat','lon'),tn_10)
                    cdf_xr['TMAX_90P'] = (('lat','lon'),tx_90)
                    
                    outfile = nmme_model + '_{:03d}.nc'.format(start_date.timetuple().tm_yday)
                    cdf_xr.to_netcdf(outfile)
                    command = 'nccopy -d6 ' + outfile + ' ' + outdir + '/' + outfile
                    res = os.system(command)
                    if res == 0:
                        os.remove(outfile)
                start_date += delta

class compute_CCDI (object):

    def subplot (self, fig, gs_,count_plot,ctitle,limits,plotdir,region,anom,longitude,latitude,figure):
        add_land,add_rivers,resol
        ax_ = fig.add_subplot(gs_[count_plot], projection=ccrs.PlateCarree())        
        ax_.set_extent(limits, ccrs.PlateCarree())
        cs_ = plt.contourf(longitude, latitude, anom,levels, cmap=cmap,
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
        cbar.set_label('%', fontsize=30)
        cbar.ax.tick_params(labelsize=20)
        print(ctitle)
        print(figure)
        plt.savefig(figure, dpi=150, format='png', bbox_inches='tight')
        
    def __init__  (self,fcst_year,fcst_mon, cwd, config, ccdi):
        '''
        compute CCDIs for the month
        '''
        result = []
        plotdir_template = cwd + '/plots/{:04d}{:02d}/'
        plotdir = plotdir_template.format(fcst_year,fcst_mon)+  config["EXP"]["domain"] + '/' + config["EXP"]["lsmdir"] + '/'
        first_file = True
        for month in range (fcst_mon + 1, fcst_mon + 3):
            region = 'FAME'
            this_mon = datetime.datetime(fcst_year,month,1,1,0,0)
            yyyy = str(this_mon.year)
            mm = str(this_mon.month).zfill(2)
            figure = figure_template.format(plotdir, ccdi, yyyy + mm)
            count_plot = 0
            fig = plt.figure(figsize=(figwidth,
                                      figwidth*(nrows*(get_boundary.get(region)[3]
                                                       - get_boundary.get(region)[2]))
                                      /(ncols*(get_boundary.get(region)[1]
                                               - get_boundary.get(region)[0]))))
            gs_ = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.1)
            
            result.append([])
            for nmme_model in config['EXP']['NMME_models']:
                file_path = config['SETUP']['E2ESDIR'] + \
                    's2spost/{}/{:04d}{:02d}/cf_{}_{:04d}{:02d}_all/'.format(nmme_model, fcst_year, fcst_mon, nmme_model, fcst_year, fcst_mon)
                # Loop through the month
                start_date = datetime.datetime(fcst_year,month,1,1,0,0)
                end_date = datetime.datetime(fcst_year,month + 1,1,1,0,0)
                delta =  datetime.timedelta(days=1)
                this_arr = np.zeros((NR, NC))

                while start_date < end_date:
                                       
                    file_name = file_path + 'PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.AFRICA_PA.LIS-S2S_DD.{:04d}{:02d}{:02d}_DT.0000_DF.NC'.format(nmme_model, start_date.year, start_date.month, start_date.day)
                    print(file_name)
                    dset = xr.open_dataset(file_name)
                    if first_file:
                        latitude = dset.lat
                        longitude = dset.lon
                        first_file = False

                    if (ccdi == 'R99d') or (ccdi == 'R95d') or (ccdi == 'R90d'):
                        ens_mean = np.mean(dset.TotalPrecip_acc.values, axis=0)
                    elif ccdi == 'TN10p' or ccdi == 'CSDI':
                        ens_mean = np.mean(dset.Tair_f_min.values, axis=0)
                    elif ccdi == 'TX90p' or ccdi == 'WSDI':
                        ens_mean = np.mean(dset.Tair_f_max.values, axis=0)
                    dset.close()
                        
                    # get the cut off
                    outdir = outdir_template.format(nmme_model,fcst_mon)
                    cdf_file = outdir + '/' + nmme_model + '_{:03d}.nc'.format(start_date.timetuple().tm_yday)
                    percdata = Dataset(cdf_file)
                    pval = np.array(percdata.variables[cutoff_name.get(ccdi)][:])

                    if (ccdi == 'R99d') or (ccdi == 'R95d') or (ccdi == 'R90d') or (ccdi == 'TX90p'):
                        idx = np.where(ens_mean - pval >= 0.)
                    else:
                        idx = np.where(ens_mean - pval <= 0.)
                    this_arr[idx] = this_arr[idx] + 1.

                    start_date += delta
                idx = np.where(this_arr == 0.)
                this_arr = 100. * this_arr / monthrange(this_mon.year, this_mon.month)[1]
                this_arr[idx] = -9999.
                self.subplot(fig, gs_,count_plot,nmme_model + ': '+ yyyy + mm, get_boundary.get(region),plotdir, region, this_arr, longitude, latitude, figure)   
                result[month - fcst_mon -1].append(this_arr)

                count_plot += 1
                    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--fcst_year', required=True, help='forecast start year')
    parser.add_argument('-m', '--fcst_mon', required=True, help= 'forecast end year')
    parser.add_argument('-c', '--configfile', required=True, help='config file name')
    parser.add_argument('-w', '--cwd', required=True, help='current working directory')
    parser.add_argument('-C', '--cdf', required=False, help='compute CDFs')
    parser.add_argument('-M', '--nmme', required=False, help='NMME model for CDF computation')

    args = parser.parse_args()
    configfile = args.configfile
    fcst_year = int(args.fcst_year)
    fcst_mon = int(args.fcst_mon)
    nmme     = args.nmme
    
    cwd = args.cwd
    cdf = args.cdf

    # load config file
    with open(configfile, 'r') as file:
        config = yaml.safe_load(file)
    
    if cdf == 'Y':
        getCDFs (fcst_mon, cwd, config, nmme)
        sys.exit()
    else:
        compute_CCDI(fcst_year,fcst_mon, cwd, config,'R95d')
        compute_CCDI(fcst_year,fcst_mon, cwd, config,'R99d')
        compute_CCDI(fcst_year,fcst_mon, cwd, config,'TN10p')
        compute_CCDI(fcst_year,fcst_mon, cwd, config,'TX90p')
        
        
