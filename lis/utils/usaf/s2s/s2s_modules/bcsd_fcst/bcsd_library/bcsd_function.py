"""
#
"""

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------


import numpy as np
#import calendar
#import math
#import time
# pylint: disable=import-error
from bcsd_stats_functions import calc_stats, lookup
# pylint: enable=import-error
#

class VarLimits:
    '''
    This function adjusts minimum and maximum values of a variable to recorded max and minimum.
    Below limits are 6h based
    '''
    def clip_array (self, data_array, var_name=None, min_val=None, max_val=None,
                    missing=None, min_thres=None, precip=None):
        ''' Below limits are 6h based'''
        min_limit={'PRECTOT': 1.e-7,
                  'PS': 30000.,
                  'T2M': 180.,
                  'LWS': 10.,
                  'SLRSF': 0.,
                  'Q2M': 0.,
                  'WIND': 0.
        }

        max_limit={'PRECTOT': 0.04,
                  'PS': 110000.,
                  'T2M': 332.,
                  'LWS': 700.,
                  'SLRSF': 1367.,
                  'Q2M': 0.05,
                  'WIND': 70.
        }

        if min_thres is not None:
            return min_limit.get('PRECTOT')
        if min_val is None:
            min_val = min_limit.get(var_name)
        if max_val is None:
            max_val = max_limit.get(var_name)
        if missing is None:
            missing = -9999.

        if precip is None:
            clipped_array = np.where(data_array == missing, data_array,
                                     np.clip(data_array, min_val, max_val))
        else:
            # mask identifies values that are less than min_val but not equal to missing
            mask_lt_min = (data_array < min_val) & (data_array != missing)
            data_array[mask_lt_min] = 0.

            # mask identifies values that are greater than max_val but not equal to missing
            mask_gt_max = (data_array > max_val) & (data_array != missing)
            data_array[mask_gt_max] = max_val
            clipped_array = data_array

        return clipped_array

    def __init__ (self):
        self.precip_thres = self.clip_array(np.empty(1), 'PRECTOT', min_thres = True)

def get_index(ref_array, my_value):
    """
      Function for extracting the index of a Numpy array (ref_array)
      which value is closest to a given number.

      Input parameters:
        - ref_array: reference Numpy array
        - my_value:  floating point number

      Returned value:
        - An integer corresponding to the index
    """
    return np.abs(ref_array - my_value).argmin()

#pylint: disable=too-many-arguments
def calc_bcsd(obs_clim_all, fcst_clim_all, lead_final, target_fcst_val_arr, target_fcst_syr, \
              target_fcst_eyr, fcst_syr, ens_final, mon, bc_var, tiny):
#pylint: enable=too-many-arguments
    """ calculates bias correction """
    correct_fcst_coarse = np.ones(((target_fcst_eyr-target_fcst_syr)+1, lead_final, ens_final))*-9999.

    for lead_num in range(0, lead_final): ## Loop from lead =0 to Final Lead
        target_month = mon + lead_num ## This is the target forecast month
        ## Check for the cases when the target forecast month is in the next year
        ## (e.g. February 1983 forecast initialized in December 1982)
        if target_month>12:
            #subtracting 12 so 13 becomes 1 meaning the month of January and so on.
            target_month-=12

        # Reading observed time series; Note the 1st column is quantile time series
        obs_quant_ts, obs_clim_ts = obs_clim_all[0, :], obs_clim_all[target_month, :]

        # Reading forecast time series; Note the 1st column is quantile time series
        fcst_quant_ts, fcst_clim_ts = fcst_clim_all[0, :], fcst_clim_all[lead_num+1, :]

        ## Now calculating mean, standard deviation and skew of observed and forecast time series
        obs_mean, obs_sd, obs_skew = calc_stats(obs_clim_ts, tiny)
        fcst_mean, fcst_sd, fcst_skew = calc_stats(fcst_clim_ts, tiny)

        ## Ok, now getting started on the bias correction
        ## Note that bias correction is done seprately for each ensemble member of all years

        for fcst_yr in range(target_fcst_syr-fcst_syr, (target_fcst_eyr-fcst_syr)+1):
            for ens_num in range (0, ens_final):
                target_fcst_val = target_fcst_val_arr[fcst_yr, lead_num, ens_num]
                ## First determine the quantile for given target forecast value
                target_fcst_quant = lookup(target_fcst_val, fcst_clim_ts, fcst_quant_ts, \
                                           len(fcst_clim_ts), bc_var, 'QUAN', fcst_mean, fcst_sd, \
                                           fcst_skew, tiny)
                # Also note that QUAN helps the the function lookup determine if we are trying to
                # convert a value to quantile or VICE versa. For converting a value to quantile
                # use 'QUAN'. For converting quantile to value use 'DATA'

                ## Now using quantile above, determine the value from observed climatology
                bias_corrected_value = lookup(target_fcst_quant, obs_quant_ts, obs_clim_ts, \
                                              len(obs_clim_ts), bc_var, 'DATA', obs_mean, obs_sd, \
                                              obs_skew, tiny)

                ## This is just a hack to check we are not getting negative value of precipitation
                if (bc_var=='PRCP') and (bias_corrected_value<0):
                    print(target_fcst_val, target_fcst_quant, fcst_yr, lead_num, ens_num)

				        ## Now storing the bias corrected anomaly
                correct_fcst_coarse[fcst_yr, lead_num, ens_num] = bias_corrected_value

    return correct_fcst_coarse

def latlon_calculations(ilat_min, ilat_max, ilon_min, ilon_max, nlats, nlons, \
                        np_obs_clim_array, np_fcst_clim_array, \
                        lead_final, target_fcst_eyr, target_fcst_syr, fcst_syr, ens_final, mon, \
                        bc_var, tiny, fcst_coarse):
    """ lat lon calculations """

    correct_fcst_coarse = np.ones(((target_fcst_eyr-target_fcst_syr)+1, lead_final, ens_final, \
                                    nlats, nlons))*-9999.

    num_lats = ilat_max-ilat_min+1
    num_lons = ilon_max-ilon_min+1

    print("num_lats = ", num_lats, np_obs_clim_array.shape)
    print("num_lons = ", num_lons, fcst_coarse.shape)

    for ilat in range(num_lats):
        lat_num = ilat_min + ilat
        for ilon in range(num_lons):
            lon_num = ilon_min + ilon
            #count_grid = ilon + ilat*num_lons

            obs_clim_all = np_obs_clim_array[:, :, ilat, ilon]
            fcst_clim_all = np_fcst_clim_array[:, :, ilat, ilon]
            target_fcst_val_arr = fcst_coarse[:, :, :, lat_num, lon_num]

            correct_fcst_coarse[:, :, :, lat_num, lon_num] = calc_bcsd(obs_clim_all, \
                                                                       fcst_clim_all, lead_final, \
                                                                       target_fcst_val_arr, \
                                                                       target_fcst_syr, \
                                                                       target_fcst_eyr, fcst_syr, \
                                                                       ens_final, mon, \
                                                                       bc_var, tiny)

    return correct_fcst_coarse
