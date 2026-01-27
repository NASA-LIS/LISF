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

BB = 2016
TDRY = 287.
LIS_CONST_G = 9.80616
LIS_CONST_TKFRZ = 273.16
LAPSE_RATE = -0.0065
DEG2RAD = math.pi / 180.0

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

def read_geosv3_elevation(static_file, logger):
    ''' reads GEOS5 surface geopotential (PHIS) and returns the height at the surface '''
    phyis = load_ncdata(static_file, logger, var_name='PHIS').isel(time = 0)
    return phyis/LIS_CONST_G

def apply_lapse_rate_correction_to_data():
    # HydroSFS Names
    tmp_name = get_hydrosfs_name('T')
    hum_name = get_hydrosfs_name('Q')
    lwd_name = get_hydrosfs_name('LWdown')
    prs_name = get_hydrosfs_name('PS')

    for itime in range(ds_in.time.size):

        # Convert to numpy values
        force_tmp = ds_in[tmp_name].values[itime,:,:]
        force_hum = ds_in[hum_name].values[itime,:,:]
        force_lwd = ds_in[lwd_name].values[itime,:,:]
        force_prs = ds_in[prs_name].values[itime,:,:]

        # Temperature
        tcforce = force_tmp + (lapse*elevdiff)

        # Pressure
        tbar = (force_tmp + tcforce)/2
        pcforce = force_prs / (np.exp((LIS_CONST_G * elevdiff) / (rdry * tbar)))

        # Humidity
        force_hum = np.where(force_hum == 0, 1e-08, force_hum)
        ee = (force_hum * force_prs) / 0.622               
        esat = 611.2 * np.exp(
            (17.67 * (force_tmp - LIS_CONST_TKFRZ)) / ((force_tmp - LIS_CONST_TKFRZ) + 243.5)
        )
        qsat = (0.622 * esat) / (force_prs - (0.378 * esat))
        rh = (force_hum / qsat) * 100.0
        fesat = 611.2 * np.exp(
            (17.67 * (tcforce - LIS_CONST_TKFRZ))/((tcforce - LIS_CONST_TKFRZ) + 243.5)
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
        
        # Output
        ds_out[tmp_name].values[itime,:,:] = tcforce
        ds_out[hum_name].values[itime,:,:] = hcforce
        ds_out[lwd_name].values[itime,:,:] = lcforce
        ds_out[prs_name].values[itime,:,:] = pcforce

    return ds_out

def apply_slope_aspect_correction_to_data(ds_in : xr.Dataset):
    
    ds_out = ds_in.copy(deep=True)
    
    # HydroSFS Names
    swd_name = get_hydrosfs_name('SWdown')
    
    # Constants
    deg2rad = math.pi / 180.0

    # Read in Sl
    file_directory_merit = f"{dir_supplementary}/lis_darun"
    file_path_merit = f"{file_directory_merit}/{file_name_landmask}"
    ds_merit = xr.open_dataset(file_path_merit)

    lon_d  = ds_merit.lon.values
    lat_d  = ds_merit.lat.values
    lat_r  = lat_d*deg2rad
    slope  = ds_merit['SLOPE'].values / deg2rad
    aspect = ds_merit['ASPECT'].values / deg2rad

    for index_time in range(ds_in.time.size):
        # Setup variables
        swd = ds_in[swd_name].isel(time = index_time).values

        # We include the full calculation here, even though minute & second are 0
        timestamp = pd.Timestamp(ds_in.time.isel(time = index_time).values)
        gmt = timestamp.hour + timestamp.minute / 60 + timestamp.second / 3600
        hr = timestamp.hour
        yr = timestamp.year
        mn = timestamp.minute
        doy = timestamp.day_of_year

        # Calculate time zone (1-24)
        zone = np.vectorize(get_timezone)(lon_d)

        # Calculate local hour (0-23) 0 = midnight, 23 = 11:00 p.m.
        change = zone - 13
        lhour = gmt + change
        lhour = np.where(lhour < 0, lhour + 24, lhour)
        lhour = np.where(lhour > 23, lhour - 24, lhour)

        # Generate cosz and decl
        gamma = 2 * math.pi * (doy - 1) / 365.
        decl = (
            0.006918
            - 0.399912 * math.cos(gamma)
            + 0.070257 * math.sin(gamma)
            - 0.006758 * math.cos(2 * gamma)
            + 0.000907 * math.sin(2 * gamma)
            - 0.002697 * math.cos(3 * gamma)
            + 0.00148 * math.sin(3 * gamma)
        )
        
        # Equation of Time
        et = 229.18 * (
            0.000075
            + 0.001868 * math.cos(gamma)
            - 0.032077 * math.sin(gamma)
            - 0.014615 * math.cos(2 * gamma)
            - 0.04089 * math.sin(2 * gamma)
        )

        #
        ls = ((zone - 1) * 15) - 180.
        lcorr = 4.*(ls - lon_d) * (-1)

        #
        latime = lhour + lcorr / 60. + et / 60.
        latime = np.where(latime < 0, latime + 24, latime)
        latime = np.where(latime > 24, latime - 24, latime)

        #
        omegad = (latime - 12) * (-15)
        omega = omegad * deg2rad

        #
        cosz = math.sin(decl) * np.sin(lat_r) + math.cos(decl) * np.cos(lat_r) * np.cos(omega)
        cosz = np.where(cosz < 0, 0, cosz)
        cosz = np.where(cosz > 1, 1, cosz)

        #
        sunang = np.where(cosz < 0.01764, 0.01764, cosz)
        cloud  = (1160.0 * sunang - swd) / (963.0 * sunang)
        cloud  = np.where(cloud < 0, 0, cloud)
        cloud  = np.where(cloud > 1, 1, cloud)

        #
        difrat = np.where(abs(sunang - 0.0223) > 0.0001, 0.0604 / (sunang - 0.0223) + 0.0683, 1)
        difrat = np.where(difrat < 0, 0, difrat)
        difrat = np.where(difrat > 1, 1, difrat)
        difrat = difrat + (1.0 - difrat) * cloud

        #
        vnrat = (580.0 - cloud * 464.0) / ((580.0 - cloud * 499.0) + (580.0 - cloud * 464.0))
        swddirect  = swd * ((1.0 - difrat) * vnrat + (1.0 - difrat)*(1.0 - vnrat))
        swddiffuse = swd * (difrat * vnrat + difrat * (1.0 - vnrat))

        #
        thour = 0 if (hr - 24) <= 0 else hr
        mody = yr - math.floor(yr * .25)*4
        if abs(mody) > 0:
            fyear = (2.0*math.pi / 365.0) * (doy - 1.0 + (thour - 12.0) / 24.0)
        else:
            fyear = (2.0*math.pi / 366.0) * (doy - 1.0 + (thour - 12.0) / 24.0)

        # Calculate the equation of time in minutes
        eqtime = 229.18 * (
            7.5e-5
            + 1.868e-3 * math.cos(fyear)
            - 3.2077e-2 * math.sin(fyear)
            - 1.4615e-2 * math.cos(2 * fyear)
            - 4.0849e-2 * math.sin(2 * fyear)
        )

        # Calculate the true solar time 
        time_offset = eqtime + 4.0 * lon_d + 60.0 * abs(lhour - hr)
        tst = thour * 60.0 + mn + time_offset

        # Solar hour angle
        ha = (tst * 0.25 - 180.0) * deg2rad

        cosphi = np.sin(lat_r) * math.sin(decl) + np.cos(lat_r) * math.cos(decl) * np.cos(ha)
        phi = np.acos(cosphi)

        #
        costheta = np.where(np.cos(lat_r) * np.sin(phi) == 0,
                            np.where(np.sin(lat_r)*cosphi - math.sin(decl) > 0, 1, -1),
                            np.sin(lat_r) * cosphi - math.sin(decl)) / (np.cos(lat_r) * np.sin(phi))
        
        # Avoid floating point errors
        costheta = np.where(abs(costheta) > 1,
                            np.where(costheta > 0, 1, -1),
                            costheta)

        #
        saz = np.where(lat_r >= 0,
                       np.where(ha < 0,
                                180.0 - np.acos(costheta) / deg2rad,
                                180.0 + np.acos(costheta) / deg2rad),
                       np.where(ha < 0,
                                np.acos(costheta) / deg2rad,
                                360 - np.acos(costheta) / deg2rad))


        #
        aslope = np.vectorize(min)(1.57, slope * deg2rad)
        aslope = np.vectorize(max)(0, aslope)
        solzen = np.acos(cosz)

        swddirect = np.where((swd > 0) & (aslope > 0) & (aslope < 90 * deg2rad),
                       correct_swddirect(aslope, aspect, saz, solzen, swddirect),
                       swddirect)

        #
        ds_out[swd_name].values[index_time, :, :] = swddirect + swddiffuse

    return ds_out

    
    





