import os
import calendar
from datetime import datetime, date
from dateutil.relativedelta import relativedelta
import argparse
import xarray as xr
import numpy as np
import yaml
# pylint: disable=import-error
import plot_utils
# pylint: enable=import-error

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

plotdir_template = cwd + '/s2splots/{:04d}{:02d}/' + config["EXP"]["lsmdir"] + '/'
plotdir = plotdir_template.format(fcst_year, fcst_mon)
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

infile_template = '{}/{}_ANOM_init_weekly_{:02d}_{:04d}.nc'
figure_template = '{}/NMME_{}_weekly_anom_{}-{}.png'

lead_week = [0, 1, 2, 3, 4, 5]

data_dir_template = cwd + '/s2smetric/{:04d}{:02d}/DYN_ANOM/' + \
        config["EXP"]["lsmdir"] + '/'
data_dir = data_dir_template.format(fcst_year, fcst_mon)

nrows = 1
ncols = 1
clabel = 'Anomaly (K)'
cartopy_dir = config['SETUP']['supplementarydir'] + '/s2splots/share/cartopy/'

#for var_name in config_["POST"]["metric_vars"]:
for var_name in ['AirT']:
    infile = infile_template.format(data_dir, '*_' + var_name, fcst_mon, fcst_year)
    print("Reading infile {}".format(infile))
    anom = xr.open_mfdataset(infile, concat_dim='ens', combine='nested')
    median_anom = np.median(anom.anom.values, axis=0)

    BEGDATE = date(fcst_year, fcst_mon, 1)
    for lead in lead_week:
        ENDDATE = BEGDATE + relativedelta(days=6)
        titles = [var_name + ' '+  BEGDATE.strftime("%Y%m%d") + '-' + ENDDATE.strftime("%Y%m%d")]
        plot_arr = median_anom[lead, ]
        print(np.nanmin(plot_arr), np.nanmax(plot_arr))
        figure = figure_template.format(plotdir, var_name, BEGDATE.strftime("%Y%m%d"), ENDDATE.strftime("%Y%m%d"))
        print(BEGDATE.strftime("%Y%m%d"),' ' , ENDDATE.strftime("%Y%m%d"))
        plot_utils.contours (anom.longitude.values, anom.latitude.values, nrows,
                             ncols, np.expand_dims(plot_arr, axis=0), 'clim_reanaly', titles, [-90, 90, -180., 180.],
                             figure, ["white", "white"],
                             fscale=0.8, clabel=clabel, levels=np.arange(-32,32+1,2),
                             cartopy_datadir=cartopy_dir, projection="mol")
        BEGDATE += relativedelta(days=7)
