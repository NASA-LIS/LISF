"""
# Author: Shrad Shukla
# coding: utf-8
#Author: Shrad Shukla
"""
import sys
import math
import numpy as np
import xarray as xr

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
                  'LWGAB': 10.,
                  'SWGDN': 0.,
                  'QV2M': 0.,
                  'WIND': 0.
        }

        max_limit={'PRECTOT': 0.04,
                  'PS': 110000.,
                  'T2M': 332.,
                  'LWGAB': 700.,
                  'SWGDN': 1367.,
                  'QV2M': 0.05,
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
        val = min(val, vec2[0])

    elif query > vec1[dim-1]: #/* if query falls above maximum value in vector 1 */
        if lu_type == 'QUAN':
            if par == 'PRCP':
                if mean == 0:
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

        val = max(val, vec2[dim - 1])

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
        try:
            val = max(val, 0)
        except NameError:
            val = tiny
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

def apply_regridding_with_mask(data, regridder, source_land_mask,
                               target_land_mask=None, method_type='conservative'):
    """
    Apply land mask and regrid the data.

    Parameters:
    -----------
    data : xarray.DataArray or xarray.Dataset
    regridder : xesmf.Regridder
    source_land_mask : xarray.Dataset (source grid land mask)
    target_land_mask : xarray.Dataset optional (target grid land mask)
    method_type : str ('bilinear' or 'conservative')
    
    Returns:
    --------
    xarray.DataArray or xarray.Dataset
    """
    any_land = source_land_mask.LANDMASK > 0

    if isinstance(data, xr.DataArray):
        masked_data = data.where(any_land)
        result = regridder(masked_data)

    elif isinstance(data, xr.Dataset):
        masked_data = data.copy(deep=True)
        for var_name in masked_data.data_vars:
            masked_data[var_name] = masked_data[var_name].where(any_land)
        result = regridder(masked_data)

    else:
        raise TypeError(f"Expected xarray.DataArray or xarray.Dataset, got {type(data)}")

    # Apply target land mask for bilinear interpolation to remove ocean extrapolation
    if method_type == 'bilinear' and target_land_mask is not None:
        target_mask = target_land_mask.LANDMASK > 0

        if isinstance(result, xr.DataArray):
            result = result.where(target_mask)
        elif isinstance(result, xr.Dataset):
            for var_name in result.data_vars:
                result[var_name] = result[var_name].where(target_mask)

    return result
