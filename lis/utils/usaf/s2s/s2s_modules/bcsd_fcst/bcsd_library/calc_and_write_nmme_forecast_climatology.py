#!/usr/bin/env python
"""
# Author: Shrad Shukla
# Usage: This is a module for the BCSD code that calculates and writes climatolgy of forecasts data.
# Output file will contain total LEAD + 1 column where LEAD can be provided as input (e.g. 6 months
# or 9 months). First column is sorted quantiles ranging from about 0 to 1 and the rest of the
# columns are just sorted climatology of forecasts for all lead times for a given forecast
# initialization month
"""
from __future__ import division
import os
import sys
import calendar
from datetime import datetime
import numpy as np
import xarray as xr
import yaml
from dateutil.relativedelta import relativedelta
# pylint: disable=import-error
from shrad_modules import read_nc_files
from bcsd_stats_functions import get_domain_info
# pylint: enable=import-error

# This function takes in a time series as input and provides sorted times series of values
# and qunatiles in return
## Using Plotting position formula to get quantiles
def create_sorted_ts(data_ts):
    """ creates sorted time series """
    clim_sort = np.sort(data_ts) ## np.sort a function to sort data into ascending order

    ## Now storing plotting position in an array
    ## formula for plotting position is (1/(n+1))
    quant_ts = np.arange(1, len(data_ts)+1)/(len(data_ts)+1)
    ## np.arange creates an array where values range from 1 to the length of time series and then
    ## simply dividing every value by length of the time series (n) + 1
    ## Now passing on sorted observed and forecast value and quantile time series
    return clim_sort, quant_ts

# Directory and file addresses
CMDARGS = str(sys.argv)
VAR = str(sys.argv[1])
INIT_FCST_MON = int(sys.argv[2])
MODEL_NAME=str(sys.argv[3])
CONFIGFILE=str(sys.argv[4])
INDIR=str(sys.argv[5])
OUTDIR = str(sys.argv[6])

with open(CONFIGFILE, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)

lat1, lat2, lon1, lon2 = get_domain_info(CONFIGFILE, extent=True)
LATS, LONS = get_domain_info(CONFIGFILE, coord=True)
CLIM_SYR = config['BCSD']['clim_start_year']
CLIM_EYR = config['BCSD']['clim_end_year']
LEAD_FINAL = config['EXP']['lead_months']
MASK_FILE = config['SETUP']['supplementarydir'] + '/bcsd_fcst/ex_raw_nmme_download.nc'
ensemble_sizes = config['EXP']['ensemble_sizes'][0]
ENS_NUM = ensemble_sizes[MODEL_NAME]

MONTH_NAME = calendar.month_abbr[INIT_FCST_MON].lower() + "01"
INFILE_TEMPLATE = '{}/{:04d}/ens{:01d}/{}.nmme.monthly.{:04d}{:02d}.nc'
OUTFILE_TEMPLATE = '{}/{}_fcst_clim.nc'

os.makedirs(OUTDIR, exist_ok=True)

### First read all forecast data
## Storing climatology of forecast for all years and ensemble members
FCST_TS = np.empty((LEAD_FINAL, ((CLIM_EYR-CLIM_SYR)+1)*ENS_NUM, len(LATS), len(LONS)))
for LEAD_NUM in range(0, LEAD_FINAL): ## Loop from lead =0 to Final Lead
    COUNT_DATA = 0
    for ens in range(ENS_NUM):
        ens1=1+ens
        for INIT_FCST_YEAR in range(CLIM_SYR, CLIM_EYR+1):
            ## Reading forecast file
            FCST_DATE = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1) + relativedelta(months=LEAD_NUM)
            FCST_YEAR, FCST_MONTH = FCST_DATE.year, FCST_DATE.month
            INFILE = INFILE_TEMPLATE.format(INDIR, INIT_FCST_YEAR, ens1, MONTH_NAME, \
                                            FCST_YEAR, FCST_MONTH)
            print(INFILE)
            FCST_TS[LEAD_NUM, COUNT_DATA, ] = read_nc_files(INFILE, VAR)
            COUNT_DATA+=1

CLIM_ARRAY = np.empty((LEAD_FINAL+1, ((CLIM_EYR-CLIM_SYR)+1)*ENS_NUM, len(LATS), len(LONS)))
## Now generating forecast climatology for all grid cells in the mask
for lat_num in range(0, len(LATS)):
	for lon_num in range(0, len(LONS)):
		## Only work with grid cells that are within the given mask
		if ((lat1<=LATS[lat_num]) and (LATS[lat_num]<=lat2) and (lon1<=LONS[lon_num]) and
                         (LONS[lon_num]<=lon2)):
			## 1st column is for sorted quantile array
			## and then rest tweleve columns have sorted climatology values for all year
			for LEAD_NUM in range(0, LEAD_FINAL):
				## Now sorting climtology time series for the given lead time
				CLIM_ARRAY[LEAD_NUM+1, :, lat_num, lon_num], CLIM_ARRAY[0, :,
                                 lat_num, lon_num] = create_sorted_ts(FCST_TS[LEAD_NUM, :,
                                                     lat_num, lon_num])

## finished storing climatology for all months
OUTFILE = OUTFILE_TEMPLATE.format(OUTDIR, VAR)
## opening outfile to write
print (CLIM_ARRAY.min(), CLIM_ARRAY.max())
# Now writing output file
CLIM_XR = xr.Dataset()
CLIM_XR['clim'] = (('DIST', 'time', 'latitude', 'longitude'), CLIM_ARRAY)
CLIM_XR.coords['DIST'] = (('DIST'), np.arange(LEAD_FINAL+1))
CLIM_XR.coords['time'] = (('time'), np.arange(((CLIM_EYR-CLIM_SYR)+1)*ENS_NUM))
CLIM_XR.coords['latitude'] = (('latitude'), LATS)
CLIM_XR.coords['longitude'] = (('longitude'), LONS)

#Now selecting only the data over the given lat lon box
SLICED_CLIM_XR = CLIM_XR.sel(longitude=slice(lon1, lon2), latitude=slice(lat1, lat2))
SLICED_CLIM_XR.to_netcdf(OUTFILE)
