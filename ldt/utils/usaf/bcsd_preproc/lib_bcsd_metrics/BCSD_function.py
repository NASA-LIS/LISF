from __future__ import division
import pandas as pd
import numpy as np
import calendar
import os.path as op
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta
from scipy.stats import percentileofscore 
from scipy.stats import scoreatpercentile, pearsonr
from math import *
import time
from BCSD_stats_functions import *
import xarray as xr
import os, errno


def CALC_BCSD(OBS_CLIM_ALL, FCST_CLIM_ALL, LEAD_FINAL, TARGET_FCST_VAL_ARR, TARGET_FCST_SYR, TARGET_FCST_EYR, FCST_SYR, ENS_NUM, MON, MONTH_NAME, count_grid, BC_VAR, TINY):

    CORRECT_FCST_COARSE = np.ones(((TARGET_FCST_EYR-TARGET_FCST_SYR)+1, LEAD_FINAL, ENS_NUM))*-999
    
    for LEAD_NUM in range(0, LEAD_FINAL): ## Loop from lead =0 to Final Lead
        TARGET_MONTH = MON + LEAD_NUM; ## This is the target forecast month
        ## Check for the cases when the target forecast month is in the next year (e.g. February 1983 forecast initialized in December 1982)
        if (TARGET_MONTH>12):
                TARGET_MONTH-=12 #subtracting 12 so 13 becomes 1 meaning the month of January and so on.
        ## Just checking if the lead and target month combination is working as expected
        if (count_grid==0): #Only printing the following for the first grid cell, no need to repeat
                print ("Initial forecast month is {} Lead is {} and Target month is {}".format(MONTH_NAME, LEAD_NUM, calendar.month_name[TARGET_MONTH]))
						
        # Retriving Observed and forecast time series for given target month
        OBS_QUANT_TS, OBS_CLIM_TS = OBS_CLIM_ALL[0, :], OBS_CLIM_ALL[TARGET_MONTH, :] ## Note that the first column is quantile time series
        FCST_QUANT_TS, FCST_CLIM_TS = FCST_CLIM_ALL[0, :], FCST_CLIM_ALL[LEAD_NUM+1, :] ## Note that the first column is quantile time series
					
        ## Now calculating mean, standard deviation and skew of both observed and forecast time series
        obs_mean, obs_sd, obs_skew = Calc_Stats(OBS_CLIM_TS, TINY)
        fcst_mean, fcst_sd, fcst_skew = Calc_Stats(FCST_CLIM_TS, TINY)
        #obs_mean, obs_sd, obs_skew = Calc_Stats(OBS_CLIM_TS.values, TINY)
        #fcst_mean, fcst_sd, fcst_skew = Calc_Stats(FCST_CLIM_TS.values, TINY)
					
        ## Ok, now getting started on the bias correction
        ## Note that bias correction is done seprately for each ensemble member of all years
					
        for fcst_yr in range(TARGET_FCST_SYR-FCST_SYR, (TARGET_FCST_EYR-FCST_SYR)+1):
                for ens_num in range (0, ENS_NUM):
                        TARGET_FCST_VAL = TARGET_FCST_VAL_ARR[fcst_yr, LEAD_NUM, ens_num]
                        ## First determine the quantile for given target forecast value
                        TARGET_FCST_QUANT = lookup(TARGET_FCST_VAL, FCST_CLIM_TS, FCST_QUANT_TS, len(FCST_CLIM_TS), BC_VAR, 'QUAN', fcst_mean, fcst_sd, fcst_skew, TINY);
                        #TARGET_FCST_QUANT = lookup(TARGET_FCST_VAL, FCST_CLIM_TS.values, FCST_QUANT_TS.values, len(FCST_CLIM_TS.values), BC_VAR, 'QUAN', fcst_mean, fcst_sd, fcst_skew, TINY);
                        ## Also note that QUAN helps the the function lookup determine if we are trying to convert a value to quantile or VICE versa
                        ## For converting a value to quantile use 'QUAN' for converting quantile to value use 'DATA'
                        ## Now using the quantile above determine the corresponding value from the observed climatology
                        BIAS_CORRECTED_VALUE = lookup(TARGET_FCST_QUANT, OBS_QUANT_TS, OBS_CLIM_TS, len(OBS_CLIM_TS), BC_VAR, 'DATA', obs_mean, obs_sd, obs_skew, TINY);
                        #BIAS_CORRECTED_VALUE = lookup(TARGET_FCST_QUANT, OBS_QUANT_TS.values, OBS_CLIM_TS.values, len(OBS_CLIM_TS.values), BC_VAR, 'DATA', obs_mean, obs_sd, obs_skew, TINY);
							
                        if (BC_VAR=='PRCP') and (BIAS_CORRECTED_VALUE<0): ## This is just a hack to check we are not getting negative value of precipitation
                                print (TARGET_FCST_VAL, TARGET_FCST_QUANT, fcst_yr, LEAD_NUM, ens_num)
							
				## Now storing the bias corrected anomaly 
                        CORRECT_FCST_COARSE[fcst_yr, LEAD_NUM, ens_num] = BIAS_CORRECTED_VALUE

    return CORRECT_FCST_COARSE

def latlon_calculations(ilat_min, ilat_max, ilon_min, ilon_max, nlats, nlons, \
                        np_OBS_CLIM_ARRAY, np_FCST_CLIM_ARRAY, \
                        LEAD_FINAL, TARGET_FCST_EYR, TARGET_FCST_SYR, FCST_SYR, ENS_NUM, MON, \
                        MONTH_NAME, BC_VAR, TINY, FCST_COARSE):

    CORRECT_FCST_COARSE = np.ones(((TARGET_FCST_EYR-TARGET_FCST_SYR)+1, LEAD_FINAL, ENS_NUM, nlats, nlons))*-999

    num_lats = ilat_max-ilat_min+1
    num_lons = ilon_max-ilon_min+1

    print("num_lats = ", num_lats, np_OBS_CLIM_ARRAY.shape)
    print("num_lons = ", num_lons, FCST_COARSE.shape)

    for ilat in range(num_lats):
        lat_num = ilat_min + ilat
        for ilon in range(num_lons):
            lon_num = ilon_min + ilon
            count_grid = ilon + ilat*num_lons

            OBS_CLIM_ALL = np_OBS_CLIM_ARRAY[:, :, ilat, ilon]
            FCST_CLIM_ALL = np_FCST_CLIM_ARRAY[:, :, ilat, ilon]
            TARGET_FCST_VAL_ARR = FCST_COARSE[:, :, :, lat_num, lon_num]

            CORRECT_FCST_COARSE[:, :, :, lat_num, lon_num] = CALC_BCSD(OBS_CLIM_ALL, FCST_CLIM_ALL, LEAD_FINAL, \
                                                                       TARGET_FCST_VAL_ARR, TARGET_FCST_SYR, \
                                                                       TARGET_FCST_EYR, FCST_SYR, ENS_NUM, MON, \
                                                                       MONTH_NAME, count_grid, BC_VAR, TINY)

    return CORRECT_FCST_COARSE
