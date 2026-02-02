"""
# Author: Shrad Shukla
# coding: utf-8
#Author: Shrad Shukla
"""
import sys
import math
import numpy as np
import xarray as xr
from ghis2s.shared.utils import load_ncdata

class VarLimits:
    '''
    This function adjusts minimum and maximum values of a variable to recorded max and minimum.
    Below limits are 6h based
    '''
    def __init__ (self, enable_rounding=False, missing_value=-9999.):
        self.precip_thres = self.clip_array(np.empty(1), 'PRECTOT', min_thres = True)
        self.enable_rounding = enable_rounding
        self.missing_value = missing_value

    def clip_array (self, data_array, var_name=None, min_val=None, max_val=None,
                    missing=None, min_thres=None, precip=None, round_decimals=None):
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

        # Optimal rounding for each variable (balances precision vs file size)
        round_precision = {
            'PS': 1,        # ~10 Pa = 0.1 hPa
            'T2M': 3,       # ~0.001 K
            'LWGAB': 2,     # ~0.01 W/m²
            'SWGDN': 2,     # ~0.01 W/m²
            'QV2M': 8,      # ~1e-8 kg/kg
        }

        if min_thres is not None:
            return min_limit.get('PRECTOT')
        if min_val is None:
            min_val = min_limit.get(var_name)
        if max_val is None:
            max_val = max_limit.get(var_name)
        if missing is None:
            missing = self.missing_value

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

        # Apply rounding if specified
        if self.enable_rounding and var_name in round_precision:
            round_decimals = round_precision[var_name]
        
        if round_decimals is not None:
            # Only round non-missing values
            mask_valid = clipped_array != missing
            clipped_array = np.where(mask_valid,
                                    np.round(clipped_array, round_decimals),
                                    clipped_array)
            
        return clipped_array

    def clip_forcing_variables(self, ds_out, var_configs):
        """
        Efficiently clip forcing variables using xarray operations.
        """
        # Clip each variable
        for var_name, kwargs in var_configs.items():
            if var_name in ds_out:
                limit_name = kwargs.get('var_name', var_name)
                clipped = xr.apply_ufunc(
                    self.clip_array,
                    ds_out[var_name],
                    kwargs={'var_name': limit_name, **kwargs},
                    dask='parallelized',
                    output_dtypes=[ds_out[var_name].dtype]
                )
                
                # Compute immediately and replace
                ds_out[var_name] = clipped

        return ds_out

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

def add_fcorr_vars(ds):
    """
    Add LHOUR, DOY and LAT to input xarray dataset.
    """
    import pandas as pd
    def get_utc_hour(time_array):
        """Extract UTC hour from datetime64 array"""
        hours = np.zeros(len(time_array), dtype=np.float32)
        for i, t in enumerate(time_array):
            t_pd = pd.Timestamp(t)
            hours[i] = t_pd.hour + t_pd.minute / 60.0 + t_pd.second / 3600.0
        return hours

    def get_day_of_year(time_array):
        """Extract day of year from datetime64 array"""
        doy = np.zeros(len(time_array), dtype=np.int16)
        for i, t in enumerate(time_array):
            t_pd = pd.Timestamp(t)
            doy[i] = t_pd.dayofyear
        return doy
    
    # Calculate UTC hours for all time steps
    utc_hours = xr.DataArray(
        get_utc_hour(ds.time.values),
        coords={'time': ds.time},
        dims=['time']
    )
    
    # Calculate solar offset from longitude
    solar_offset = ds.lon / 15.0
    
    # utc_hours: (time,) + solar_offset: (lon,) -> (time, lon)
    # Then broadcast to (time, lat, lon)
    lhour = utc_hours + solar_offset
    
    # Wrap to 0-24 range
    lhour = xr.where(lhour < 0, lhour + 24, lhour)
    lhour = xr.where(lhour >= 24, lhour - 24, lhour)
    lhour = lhour.expand_dims({'lat': ds.lat})
    lhour = lhour.transpose('time', 'lat', 'lon')
    
    lhour.attrs = {
        'long_name': 'Local Solar Hour',
        'units': 'hours',
        'description': 'Local solar hour (0-24) based on longitude, independent of political time zones',
        'valid_range': [0, 24]
    }
    ds['LHOUR'] = lhour
    
    # Add Day of Year (DOY)
    doy = xr.DataArray(
        get_day_of_year(ds.time.values),
        coords={'time': ds.time},
        dims=['time']
    )
    doy = doy.expand_dims({'lat': ds.lat, 'lon': ds.lon})
    doy = doy.transpose('time', 'lat', 'lon')

    doy.attrs = {
        'long_name': 'Day of Year',
        'units': 'day',
        'description': 'Day of year (1-366)',
        'valid_range': [1, 366]
    }
    ds['DOY'] = doy

    # Add 2D latitude array for vectorized operations
    lat_2d = xr.DataArray(
        np.broadcast_to(ds.lat.values[:, np.newaxis], (len(ds.lat), len(ds.lon))),
        coords={'lat': ds.lat, 'lon': ds.lon},
        dims=['lat', 'lon']
    )
    
    lat_2d.attrs = {
        'long_name': 'Latitude (2D)',
        'units': 'degrees_north',
        'description': '2D latitude array for vectorized calculations'
    }
    ds['LAT_2D'] = lat_2d
    
    return ds

def apply_fcorr(mask, slope, aspect, elevdiff, lat,
                force_tmp, force_hum, force_lwd, force_prs,
                swd, lhour, doy, force_corr=None):
    '''
    Apply forcing corrections (lapse rate and aspect).
    
    Parameters:
    -----------
    mask : float - land mask value [point]
    slope : float - terrain slope (degrees) [point]
    aspect : float - terrain aspect (degrees) [point]
    elevdiff : float - elevation difference (m) [point]
    lat : float - latitude (degrees) [point]
    force_tmp : array (time,) - air temperature (K) [n_times]
    force_hum : array (time,) - specific humidity (kg/kg) [n_times]
    force_lwd : array (time,) - longwave radiation (W/m2) [n_times]
    force_prs : array (time,) - surface pressure (Pa) [n_times]
    swd : array (time,) - shortwave radiation (W/m2)[n_times] 
    lhour : array (time,) - local solar hour (0-24) [n_times]
    doy : array (time,) - day of year (1-366) [n_times]
    force_corr : dict - correction flags {'lapsrate': bool, 'aspect': bool}
    
    Returns:
    --------
    corrected_force : array (n_outvar, time) - corrected forcing variables
    '''
    # constants
    bb = 2016
    tdry = 287.
    lis_g = 9.80616
    lis_tfrz = 273.16
    lapserate = -0.0065
    deg2rad = np.pi / 180.0
    
    def correct_swddirect(aslope, aspect, saz, solzen, swddirect):
        delaz = abs(aspect - saz) * deg2rad
        sdircorr = np.sin(solzen) * np.sin(aslope) * np.cos(delaz)
        sdircorr = np.where(solzen <= 85 * deg2rad, sdircorr / np.cos(solzen), sdircorr / np.cos(85))
        swddirect = swddirect * (np.cos(aslope) + sdircorr)
        swddirect = np.where(swddirect < 0, 0, swddirect)
        return swddirect
    
    n_times = len(force_tmp)
    n_outvar = 0
    var_no = 0

    if force_corr['lapsrate']:
        n_outvar += 4
    if force_corr['aspect']:
        n_outvar += 1

    corrected_force = np.ones((n_outvar, n_times), dtype=np.float32)*-9999.
    if mask > 0:
        # (1) lapse rate correction
        if force_corr['lapsrate']:
            # Temperature
            tcforce = force_tmp + (lapserate*elevdiff)
            
            # Pressure
            tbar = (force_tmp + tcforce)/2
            pcforce = force_prs / (np.exp((lis_g * elevdiff) / (tdry * tbar)))
            
            # Humidity
            force_hum = np.where(force_hum == 0, 1e-08, force_hum)
            ee = (force_hum * force_prs) / 0.622               
            esat = 611.2 * np.exp(\
                                  (17.67 * (force_tmp - lis_tfrz)) / ((force_tmp - lis_tfrz) + 243.5)
                                  )
            qsat = (0.622 * esat) / (force_prs - (0.378 * esat))
            rh = (force_hum / qsat) * 100.0
            fesat = 611.2 * np.exp(\
                                   (17.67 * (tcforce - lis_tfrz))/((tcforce - lis_tfrz) + 243.5)
                                   )
            fqsat = (0.622 * fesat) / (pcforce - (0.378 * fesat))
            hcforce = (rh * fqsat) / 100.0
            
            # Longwave Radiation
            fe = (hcforce * pcforce) / 0.622
            mee = ee / 100.0
            mfe = fe / 100.0
        
            # Correct for negative vapor pressure at very low temperatures at high latitudes
            mee = np.where(mee < 0, 1e-08, mee)
            mfe = np.where(mfe < 0, 1e-08, mfe)
            emiss  = 1.08 * (1 - np.exp(-mee**(force_tmp / bb)))
            femiss = 1.08 * (1 - np.exp(-mfe**(tcforce / bb)))
            ratio = (femiss * (tcforce**4))/(emiss * (force_tmp**4))
            lcforce = force_lwd * ratio
        
            # lapse rate corrected forcings
            corrected_force[var_no + 0,:] = np.round(tcforce,2)
            corrected_force[var_no + 1,:] = np.round(hcforce,6)
            corrected_force[var_no + 2,:] = np.round(lcforce,1)
            corrected_force[var_no + 3,:] = np.round(pcforce,0)
            var_no = 4

        # (1) Aspect correction
        if force_corr['aspect']:
            lat_r = lat * deg2rad
            # Generate cosz and decl
            gamma = 2 * np.pi * (doy - 1) / 365.
            decl = (
                0.006918
                - 0.399912 * np.cos(gamma)
                + 0.070257 * np.sin(gamma)
                - 0.006758 * np.cos(2 * gamma)
                + 0.000907 * np.sin(2 * gamma)
                - 0.002697 * np.cos(3 * gamma)
                + 0.00148 * np.sin(3 * gamma)
            )

            # Solar hour angle from local solar hour
            # lhour is already the local solar hour (0-24)
            # Solar noon is at lhour=12, so hour angle is:
            hour_angle_deg = (lhour - 12.0) * 15.0  # degrees from solar noon
            omega = hour_angle_deg * deg2rad
            
            # Cosine of solar zenith angle
            cosz = np.sin(decl) * np.sin(lat_r) + np.cos(decl) * np.cos(lat_r) * np.cos(omega)
            cosz = np.where(cosz < 0, 0, cosz)
            cosz = np.where(cosz > 1, 1, cosz)
            
            # Cloud fraction estimation
            sunang = np.where(cosz < 0.01764, 0.01764, cosz)
            cloud  = (1160.0 * sunang - swd) / (963.0 * sunang)
            cloud  = np.where(cloud < 0, 0, cloud)
            cloud  = np.where(cloud > 1, 1, cloud)

            difrat = np.where(abs(sunang - 0.0223) > 0.0001, 0.0604 / (sunang - 0.0223) + 0.0683, 1)
            difrat = np.where(difrat < 0, 0, difrat)
            difrat = np.where(difrat > 1, 1, difrat)
            difrat = difrat + (1.0 - difrat) * cloud
            
            # Visible/NIR split
            vnrat = (580.0 - cloud * 464.0) / ((580.0 - cloud * 499.0) + (580.0 - cloud * 464.0))
            swddirect  = swd * ((1.0 - difrat) * vnrat + (1.0 - difrat)*(1.0 - vnrat))
            swddiffuse = swd * (difrat * vnrat + difrat * (1.0 - vnrat))

            # Using the simplified approach since we have local solar hour
            cosphi = cosz
            phi = np.arccos(cosphi)

            # Avoid division by zero
            costheta = np.where(np.cos(lat_r) * np.sin(phi) == 0,
                                np.where(np.sin(lat_r) * cosphi - np.sin(decl) > 0, 1, -1),
                                (np.sin(lat_r) * cosphi - np.sin(decl)) / (np.cos(lat_r) * np.sin(phi))
                                )

            # Avoid floating point errors
            costheta = np.where(np.abs(costheta) > 1,
                                np.where(costheta > 0, 1, -1),
                                costheta)

            # Calculate solar azimuth based on hemisphere and time of day
            # omega is the hour angle: negative before solar noon, positive after
            saz = np.where(lat_r >= 0,  # Northern hemisphere
                           np.where(omega < 0,  # Morning (before solar noon)
                                    180.0 - np.arccos(costheta) / deg2rad,
                                    180.0 + np.arccos(costheta) / deg2rad),
                           np.where(omega < 0,  # Southern hemisphere, morning
                                    np.arccos(costheta) / deg2rad,
                                    360.0 - np.arccos(costheta) / deg2rad)
                           )

            slope_rad = slope * deg2rad
            aslope = np.clip(slope_rad, 0, 1.57)
            solzen = np.arccos(cosz)
                        
            swddirect = np.where((swd > 0) & (aslope > 0) & (aslope < 90 * deg2rad),
                                 correct_swddirect(aslope, aspect, saz, solzen, swddirect),
                                 swddirect)

            corrected_force[var_no,:] = np.round(swddirect,1) + np.round(swddiffuse,1)
    
    return corrected_force

    
    





