#!/usr/bin/env python3
"""

"""

from datetime import datetime
import sys
import math
from time import ctime as t_ctime
from time import time as t_time
import numpy as np
import yaml
#import netCDF4 as nc
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
from netCDF4 import date2num as nc4_date2num
from netCDF4 import default_fillvals as nc4_default_fillvals
# pylint: enable=no-name-in-module


def calc_bcsd(obs_clim_all, fcst_clim_all, lead_final, target_fcst_val_arr, ens_final, fcst_init_month, bc_var, tiny):
    """ calculates bias correction """
    correct_fcst_coarse = np.ones((lead_final, ens_final))*-9999.

    for lead_num in range(0, lead_final): ## Loop from lead =0 to Final Lead
        target_month = fcst_init_month + lead_num ## This is the target forecast month
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

        for ens_num in range (0, ens_final):
            target_fcst_val = target_fcst_val_arr[lead_num, ens_num]
            
            if target_fcst_val == -9999.:
                continue
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
                print(target_fcst_val, target_fcst_quant, lead_num, ens_num)

                    ## Now storing the bias corrected anomaly
            correct_fcst_coarse[lead_num, ens_num] = bias_corrected_value

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


def calc_stats(data, tiny): #,int n,float *mean,float *sd,float *skew)
    """ calculates statistics """
    ## The following function estimates mean, standard deviation and skew for a given time series
    n_len = len(data)
    if n_len <= 1:
        print ("n_len must be at least 2 in function stats")
        sys.exit(1)
    sum_val=0.0
    var = 0.0
    skew = 0.0
    #first calculate the mean
    for j in range(n_len):
        sum_val+=data[j]
    mean=sum_val/n_len

    # second calculate standard deviation and skew
    for j in range(n_len):
        sum_val=data[j]-(mean)
        skew+=pow(sum_val,3)
        var+=pow(sum_val,2)

    sd_val=math.sqrt(var/(n_len-1)) ##unbiased standard deviation

    if var >0:
        skew = (skew/n_len)/pow(var/n_len, 1.5)
        ## see section 2.2.24.1 of http://library.atgti.az/categories/economy/Zwillinger%20Daniel
        ## %20-%20CRC%20Standard%20Probability%20and%20Statistics%20Tables%20and%20Formulae.pdf
        ## for skewness formular. Python scipy uses the same formula:
        ## http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.skew.html#id1
    else:
        skew=0

        ## KRA Commented out:
        # print "No skew/kurtosis when variance={} {} {}".format(var, mean, sd)
        sd_val = tiny ### HACK for when SD == 0

    return mean, sd_val, skew

def lookup(query, vec1, vec2, dim, par, lu_type, mean, sd_val, skew, tiny):
    """ lookup function for query data """
    # Here query is a data value or qunantile which needs to be converted into quantile or
    # valueLookup value of "query" in data vector vec1[], return corresponding
    #  quantile value in vec2[], both of dimension dim */
    # "par" defines whether the data is precipitation or temperature */
    # "lu" defines whether the variable being returned by the function
    #  (from vec2[]) is a quantile or data -- function call depends on this*/
    #/* if above or below range of known data, use these estimates */
    #/* if query falls below lowest value in vector 1 */

    if query < vec1[0]:
        if lu_type == 'QUAN':
            if par == 'PRCP':
                val=get_f_from_data_weibul(mean, sd_val, skew, query, tiny)
            elif par == 'TEMP':
                val=get_f_from_data_normal(mean, sd_val, query, tiny)
        else: #(lu_type==data)
            if par == 'PRCP':
                val=get_data_from_f_weibul(mean, sd_val, skew, query, tiny)
            elif par == 'TEMP':
                val=get_data_from_f_normal(mean, sd_val, query, tiny)
        if val > vec2[0]:
            val=vec2[0]
    elif query > vec1[dim-1]: #/* if query falls above maximum value in vector 1 */
        if lu_type == 'QUAN':
            if par == 'PRCP':
                ##### HACK
                if mean == 0:
## KRA Commented out:                     print ("F HACK invoked")
                    val = tiny
                else:
                    val=get_f_from_data_evi(mean, sd_val, query, tiny)
            elif par == 'TEMP':
                val=get_f_from_data_normal(mean, sd_val, query, tiny)
        if lu_type == 'DATA':
            if par == 'PRCP':
                ##### HACK
                if mean == 0:
                    print ("data HACK invoked")
                    val=0
                else:
                    val=get_data_from_f_evi(mean, sd_val, query, tiny)
            if par=='TEMP':
                val=get_data_from_f_normal(mean, sd_val, query, tiny)

        if val<vec2[dim-1]:
            val=vec2[dim-1]
    #/* otherwise, it is within the range of known data in vector 1
    # do a linear interpolation between known points */
    else:
        for ndx in range(dim-1): ## Looping from 0 to dim-2
            if vec1[ndx] <= query <= vec1[ndx+1]:
                if (vec1[ndx+1]-vec1[ndx]) == 0:
                    a_val=1
                else:
                    a_val = (vec1[ndx+1]-query)/(vec1[ndx+1]-vec1[ndx])
                val = a_val*vec2[ndx]+(1-a_val)*vec2[ndx+1]
                break
    if par=='PRCP':
        val = max(val, 0)

    try:
        # We TRY to access the variable to see if it was created.
        return val
        # print(f" DEBUG: 'val' was created. Its value is: '{val}'")
    except UnboundLocalError:
        # If the error occurs, this block runs.
        print("DEBUG: 'val' was never created because the if-condition was False.")
        print(f"query: {query}")
        print(f"vec1: {vec1}")
        print(f"vec2: {vec2}")
        print(f"dim: {dim}") 
        print(f"par: {par}")
        print(f"lu_type: {lu_type}")
        print(f"mean: {mean}")
        print(f"sd_val: {sd_val}")
        print(f"skew: {skew}")
        print(f"tiny: {tiny}")
    
    return val

def get_f_from_data_normal(mean, sd_val, x_val, tiny):
    """ estimates quantile for given value, following normal distribution """
    ## uses approximation from Handbook of Hydrology eq. 18.2.2
    sign=1
    z_val = (x_val-mean)/sd_val

    if z_val < 0:
        sign=-1
        z_val*=sign

    f_val = 1-0.5*math.exp(-1*((83*z_val+351)*z_val+562)/(165+703/z_val))

    if f_val == 1.0:
        f_val-=tiny  # adjustment so that outlier ensembles don't.
        # print ("Invoked tiny in get_f_from_data_normal")

    if f_val==0.0:
        f_val+=tiny
        # print ("Invoked tiny in get_f_from_data_normal")

    if (f_val>=1.0) or (f_val<=0.0):
        #mostly obviated by above ifs; could be used
        print("Error in get_f_from_data_normal")
        print("Bad quantile:")
        print(f"data={x_val:0.2f} mean={mean:0.2f} sd={sd_val:0.2f}")
        print(f"quant={f_val:0.2f} z={z_val:0.2f}")
        sys.exit(1)
    else:
        if sign==-1:
            f_val = 1-f_val
        return f_val

def get_data_from_f_normal(mean, sd_val, f_val, tiny):
    """ converts quantile into value, following normal distribution """
    # uses approximation from Handbook of Hydrology eq. 18.2.3a
    if(f_val>=1.0) or (f_val<=0.0):
        print ("Error in get_data_from_f_normal")
        print ("Bad quantile:")
        print (f"mean={mean:0.2f} sd={sd_val:0.2f} quant={f_val:0.2f} tiny ={tiny:0.2f}")
        sys.exit(0)

    ## Calculating z
    z_val = (pow(f_val,0.135)-pow(1-f_val,0.135))/0.1975

    ## Now doing quality control on z

    if z_val > 5:
        print ("Pushing upper limit of applicability of equation in get_data_from_f_normal.")
        print (f"Best range is 0<z<5. z = {z_val:0.2f}")

    x_val = z_val*sd_val+mean
    return x_val

def get_f_from_data_evi(mean, sd_val, x_val, tiny):
    """ estimates quantile for given value, following Gumbel EVI type-1 distribution """
    #Gumbel (EV Type I distribution) for maxima
    b_val = 3.14159/(sd_val*math.sqrt(6))
    a_val = mean-0.5772/b_val
    f_val = math.exp(-1*math.exp(-b_val*(x_val-a_val)))

    ## For outliers in quantile
    if f_val==1.0:
        f_val-=tiny
        #print ("Invoked tiny in get_f_from_data_evi")
    if f_val==0.0:
        f_val+=tiny
        #print ("Invoked tiny in get_f_from_data_evi")

    if (f_val>=1.0) or (f_val<=0.0):
        print ("Error in get_f_from_data_evi")
        print ("Bad quantile:")
        print (f"data={x_val:0.2f} mean={mean:0.2f} sd={sd_val:0.2f} quant={f_val:0.2f}")
        sys.exit(1)
    else:
        return f_val

def get_data_from_f_evi(mean, sd_val, f_val, tiny):
    """ converts quantile into value, following Gumbel EVI type-1 distribution """
    #Gumbel (EV Type I distribution) for maxima
    if(f_val>=1.0) or (f_val<=0.0):
        print ("Error in get_data_from_f_evi")
        print ("Bad quantile:")
        print (f"mean={mean:0.2f} sd={sd_val:0.2f} quant={f_val:0.2f} tiny={tiny:0.2f}")
        sys.exit(0)
    b_val = 3.14159/(sd_val*math.sqrt(6))
    a_val = mean-0.5772/b_val
    x_val = a_val-(1/b_val)*(math.log(-math.log(f_val)))
    return x_val

def get_f_from_data_weibul(mean, sd_val, skew, x_val, tiny):
    """ estimates quantile for given value, following Weibul distribution """
    #Weibull (EV Type III distribution) for minima, bounded at zero
    #approximation for a (alpha) is eq. 11-32 from Kite (1977)

    #if(skew>8.214) or (skew<-1):
    #    print ("Outside limit for table in get_data_from_f_weibul")
    #    print ("best range is -1<skew<8.2 skew= {:0.2f}".format(skew))

    alpha, a_alpha, b_alpha = weibul_params(skew)
    b_val = a_alpha*sd_val+mean
    bound = b_val-b_alpha*sd_val
    #lower bound minimum of zero, but allow larger minima if data say so
    bound = max(bound, 0)
    if bound > x_val:
        bound=0
    c_val = (x_val-bound)/(b_val-bound)
    c_val = max(c_val, 0)
        #print ('Invoked c_val=0 in get_f_from_data_weibul')
    f_val=1-math.exp(-1*pow(c_val,alpha))

    if f_val==1.0:
        # For outliers
        f_val-=tiny
        #print ("Invoked tiny in get_f_from_data_weibul")

    if f_val==0.0:
        f_val+=tiny
        #print ("Invoked tiny in get_f_from_data_weibul")

    if(f_val>=1.0) or (f_val<=0.0):
        print ("Error in get_f_from_data_weibul")
        print ("Bad quantile:")
        print (f"data={x_val:0.2f} mean={mean:0.2f} sd={sd_val:0.2f} quant={f_val:0.2f}")
        print (f"Other info: A={a_alpha:0.3f} b={b_val:0.3f} B={b_alpha:0.3f} bound={bound:0.3f}")
        sys.exit(1)
    else:
        return f_val

def get_data_from_f_weibul(mean, sd_val, skew, f_val, tiny):
    """ converts quantile into value, following Weibul distribution """
    #/* Weibull (EV Type III distribution) for minima, bounded at zero */
    #approximation for a (alpha) is eq. 11-32 from Kite (1977) */
    if(f_val>=1.0) or (f_val<=0.0):
        print("Error entering get_f_from_data_weibul")
        print("Bad quantile:")
        print(f"mean={mean:0.2f} sd={sd_val:0.2f} skew={skew:0.2f}")
        print(f"quant={f_val:0.2f} tiny={tiny:0.2f}")
        sys.exit(1)

    #if(skew>8.2) or (skew<-1):
    #    print ("Outside limit for table in get_data_from_f_weibul")
    #    print ("best range is -1<skew<8.2 skew= {:0.2f}".format(skew))

    alpha, a_alpha, b_alpha = weibul_params(skew)

    b_val = a_alpha*sd_val+mean
    bound = b_val-b_alpha*sd_val

    #/*  lower bound minimum of zero, but allow larger minima if data say so */
    bound = max(bound, 0)

    x_out = ((pow(-math.log(1-f_val),1/alpha))*(b_val-bound)) + bound
    return x_out

def weibul_params(skew):
    """ returns alpha, Aalpha, and Balpha for the Weibull distrubution """
    #table taken from Statistical Methods in Hydrology, Haan, shown below
    sk_val = [-1.000, -0.971, -0.917, -0.867, -0.638, -0.254, 0.069, 0.359, 0.631, 0.896,  1.160,  \
          1.430,  1.708,  2.000, 2.309, 2.640, 2.996,  3.382,  3.802,  4.262,  4.767,  5.323,  \
          5.938, 6.619, 7.374,  8.214]
    inva = [0.020, 0.030, 0.040, 0.050, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, \
            0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900, 2.000, \
            2.100, 2.200]
    a_vec = [0.446,  0.444,  0.442,  0.439,  0.425,  0.389,  0.346,  0.297, 0.246,  0.193,  0.142, \
             0.092,  0.044,  0.000, -0.040, -0.077, -0.109, -0.136, -0.160, -0.180, -0.196, \
             -0.208, -0.217, -0.224, -0.227, -0.229]
    b_vec = [40.005, 26.987, 20.481, 16.576, 8.737, 4.755, 3.370, 2.634, 2.159,  1.815,  1.549, \
             1.334, 1.154, 1.000, 0.867, 0.752, 0.652,  0.563,  0.486,  0.418, 0.359, 0.308, \
             0.263, 0.224, 0.190, 0.161]

    if skew>sk_val[-1]:
        skew=sk_val[-1]
    elif skew<sk_val[0]:
        skew=sk_val[0]

    for ndx in range(len(sk_val)-1): ## looping from 0 to 24
        #if(skew <= sk_val[ndx+1]) and (skew >= sk_val[ndx]):
        if sk_val[ndx] <= skew <= sk_val[ndx+1]:
            if (sk_val[ndx+1]-sk_val[ndx]) == 0:
                a_val=1
            else:
                a_val = (sk_val[ndx+1]-skew)/(sk_val[ndx+1]-sk_val[ndx])
            alpha = 1/(a_val*inva[ndx]+(1-a_val)*inva[ndx+1])
            a_alpha = a_val*a_vec[ndx]+(1-a_val)*a_vec[ndx+1]
            b_alpha = a_val*b_vec[ndx]+(1-a_val)*b_vec[ndx+1]
            break
    return alpha, a_alpha, b_alpha
