#!/usr/bin/env python
# Author: Shrad Shukla
# Usage: This is a module for the BCSD code that calculates and writes climatolgy of forecasts data.
# Output file will contain total LEAD + 1 column where LEAD can be provided as input (e.g. 6 months, 9 months). First column is sorted quantiles ranging from about 0 to 1 and the rest of the columns are just sorted climatology of forecasts for all lead times for a given forecast initialization month

from __future__ import division
import numpy as np
from shrad_modules import read_nc_files
import sys
import calendar
from datetime import datetime
from dateutil.relativedelta import relativedelta
import xarray as xr
import os

# <codecell>

## This function takes in a time series as input and provides sorted times series of values and qunatiles in return
## Using Plotting position formula to get quantiles
def Create_sorted_ts(DATA_TS):
    clim_sort = np.sort(DATA_TS) ## np.sort a function to sort data into ascending order

    ## Now storing plotting position in an array
    ## formula for plotting position is (1/(n+1))
    quant_ts = np.arange(1, len(DATA_TS)+1)/(len(DATA_TS)+1)
    ## np.arange creates an array where values range from 1 to the length of time series and then simply dividing every value by length of the time series (n) + 1
    ## Now passing on sorted observed and forecast value and quantile time series
    return clim_sort, quant_ts

cmdargs = str(sys.argv)
VAR = str(sys.argv[1])
lat1, lat2, lon1, lon2 = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
INIT_FCST_MON = int(sys.argv[6])
print("INIT FCST MON:", INIT_FCST_MON)
FCST_SYR, FCST_EYR = int(sys.argv[7]), int(sys.argv[8])
LEAD_FINAL = int(sys.argv[9]);
ENS_NUM = int(sys.argv[10])
MONTH_NAME = calendar.month_abbr[INIT_FCST_MON].lower() + "01"
print("month name is:",MONTH_NAME)
CLIM_SYR, CLIM_EYR = int(sys.argv[11]), int(sys.argv[12]) ## Both FCST_SYR and CLIM_SYR, and FCST_EYR and CLIM_EYR will be the same since we are using the entire hindcast for
# Directory and file addresses
INDIR=str(sys.argv[13])
INFILE_template = '{}/{}/{:04d}/ens{:01d}/{}.cfsv2.{:04d}{:02d}.nc'
#
OUTFILE_template = '{}/{}_fcst_clim.nc'
#
MASK_FILE = str(sys.argv[14])
LATS = read_nc_files(MASK_FILE, 'lat'); LONS = read_nc_files(MASK_FILE, 'lon');

OUTDIR = str(sys.argv[15])

#print(VAR)
#print(lat1)
#print(lat2)
#print(lon1)
#print(lon2)
#print(INIT_FCST_MON)
#print(LEAD_FINAL)
#print(ENS_NUM)
#print(MONTH_NAME)
#print(INDIR)
print(MASK_FILE)
#print(OUTDIR)

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

print ("Ready to create climatology for Variable {}".format(VAR))
print ("Forecast Initialization month is {}".format(MONTH_NAME))


### First read all forecast data
FCST_TS = np.empty((LEAD_FINAL, ((CLIM_EYR-CLIM_SYR)+1)*ENS_NUM, len(LATS), len(LONS))) ## Storing climatology of forecast for all years and ensemble members
for LEAD_NUM in range(0, LEAD_FINAL): ## Loop from lead =0 to Final Lead
	COUNT_DATA = 0
	for ens in range(ENS_NUM):
		for INIT_FCST_YEAR in range(CLIM_SYR, CLIM_EYR+1):
			## Reading forecast file
			FCST_DATE = datetime(INIT_FCST_YEAR, INIT_FCST_MON, 1) + relativedelta(months=LEAD_NUM)
			FCST_YEAR, FCST_MONTH = FCST_DATE.year, FCST_DATE.month
			INFILE = INFILE_template.format(INDIR, MONTH_NAME, INIT_FCST_YEAR, ens+1, MONTH_NAME, FCST_YEAR, FCST_MONTH)
			print (INFILE)
			FCST_TS[LEAD_NUM, COUNT_DATA, ] = read_nc_files(INFILE, VAR)[0,]
			COUNT_DATA+=1
			## read_nc_files is module that I have set up to read netcdf files you just have to pass file name and the name of the variable to read

CLIM_ARRAY = np.empty((LEAD_FINAL+1, ((CLIM_EYR-CLIM_SYR)+1)*ENS_NUM, len(LATS), len(LONS)))
## Now generating forecast climatology for all grid cells in the mask
CLIM_ARRAY[0:,:,:,:] = FCST_TS.sort(axis=1)
len_ts = ((CLIM_EYR-CLIM_SYR)+1)*ENS_NUM
CLIM_ARRAY[0,:,:,:] = np.arange(1, len_ts+1)/(len_ts+1)

#for lat_num in range(0, len(LATS)):
#	for lon_num in range(0, len(LONS)):
#		## Only work with grid cells that are within the given mask
#		if ((lat1<=LATS[lat_num]) and (LATS[lat_num]<=lat2) and (lon1<=LONS[lon_num]) and (LONS[lon_num]<=lon2)):
			## 1st column is for sorted quantile array and then rest tweleve columns have sorted climatology values for all year
#			for LEAD_NUM in range(0, LEAD_FINAL):
#				## Now sorting climtology time series for the given lead time
#				CLIM_ARRAY[LEAD_NUM+1, :, lat_num, lon_num], CLIM_ARRAY[0, :, lat_num, lon_num] = Create_sorted_ts(FCST_TS[LEAD_NUM, :, lat_num, lon_num])

## finished storing climatology for all months
OUTFILE = OUTFILE_template.format(OUTDIR, VAR)
## opening outfile to write
print ("Writing {}".format(OUTFILE))
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

