

from __future__ import division
import numpy as np
from shrad_modules import read_nc_files
import calendar
import os.path as op
import sys
from datetime import datetime
from scipy.stats import percentileofscore 
from scipy.stats import scoreatpercentile, pearsonr
from math import *
import time
import yaml

# In[11]:

#This function write 4 dimensional arrays
def get_domain_info (configfile, extent=None, coord=None):
    # load config file
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)
    sys.path.append(config['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.shared import utils

    if extent is not None:
        return utils.get_domain_info(configfile, extent=True)
    if coord is not None:
        return utils.get_domain_info(configfile, coord=True)
    
def write_4d_netcdf(infile, var, varname, DESCRIPTION, SOURCE, VAR_UNITS, SIG_DIGIT, lons, lats, ENS_NUM, LEAD_NUM, SDATE, dates):
    from datetime import datetime, timedelta
    import netCDF4 as nc
    rootgrp = nc.Dataset(infile, 'w', format='NETCDF4')
    longitude = rootgrp.createDimension('longitude', len(lons))
    latitude = rootgrp.createDimension('latitude', len(lats))
    lead = rootgrp.createDimension('Lead', LEAD_NUM)
    ens = rootgrp.createDimension('Ens', ENS_NUM)
    time = rootgrp.createDimension('time', None)
    leads = rootgrp.createVariable('Lead','d',('Lead',))
    enss = rootgrp.createVariable('Ens','d',('Ens',))
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    times = rootgrp.createVariable('time','f8',('time',))
    # two dimensions unlimited.
    varname = rootgrp.createVariable(varname,'f4',('time', 'Lead', 'Ens', 'latitude','longitude'), fill_value=nc.default_fillvals['f4'], zlib=True, complevel=6, shuffle=True)
    import time
    rootgrp.description = DESCRIPTION
    rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.source = SOURCE
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    enss.units = 'unitless'
    varname.units = VAR_UNITS
    STRING_DATE = datetime.strftime(SDATE, "%Y-%m-%d")
    times.units = 'days since ' + STRING_DATE
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    leads[:]=np.arange(0.5, LEAD_NUM+0.5)
    enss[:]=np.arange(0, ENS_NUM)
    varname[:,:,:,:,:] = var
    times[:] = nc.date2num(dates,units=times.units,calendar=times.calendar)
    rootgrp.close()


## The following function estimates mean, standard deviation and skew for a given time series
def Calc_Stats(DATA, TINY): #,int n,float *mean,float *sd,float *skew)
    n = len(DATA)
    if (n<=1):
        print ("n must be at least 2 in function stats")
        sys.exit(1);
    s=0.0; ep= 0.0; var = 0.0; skew = 0.0
    #first calculate the mean
    for j in range(n):
        s+=DATA[j];
    mean=s/n;
  
    # second calculate standard deviation and skew
    for j in range(n):
        s=DATA[j]-(mean);
        skew+=pow(s,3);
        var+=pow(s,2)
    sd=sqrt(var/(n-1)); ##unbiased standard deviation 
    
    if(var):
        skew = (skew/n)/pow(var/n, 1.5)
        ## see section 2.2.24.1 of http://library.atgti.az/categories/economy/Zwillinger%20Daniel%20-%20CRC%20Standard%20Probability%20and%20Statistics%20Tables%20and%20Formulae.pdf for skewness formular
        ## Python scipy uses the same formula http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.skew.html#id1
    else:
        skew=0;

## KRA Commented out:        print "No skew/kurtosis when variance={} {} {}".format(var, mean, sd);
        sd = TINY ### HACK for when SD == 0
		#sys.exit(1);
    
    return mean, sd, skew
def lookup(query, vec1, vec2, dim, par, lu, mean, sd, skew, TINY):
    # Here query is a data value or qunantile which needs to be converted into quantile or valueLookup value of "query" in data vector vec1[], return corresponding
    #  quantile value in vec2[], both of dimension dim */
    # "par" defines whether the data is precipitation or temperature */
    # "lu" defines whether the variable being returned by the function
    #  (from vec2[]) is a quantile or data -- function call depends on this*/
    #/* if above or below range of known data, use these estimates */
    #/* if query falls below lowest value in vector 1 */
    
    if(query < vec1[0]):
        if(lu=='QUAN'):
            if(par=='PRCP'):
                val=get_F_from_data_weibul(mean, sd, skew, query, TINY);
            elif(par=='TEMP'):
                val=get_F_from_data_normal(mean, sd, query, TINY);
        else: #(lu==DATA)
            if(par=='PRCP'):
                val=get_data_from_F_weibul(mean, sd, skew, query, TINY);
            elif(par=='TEMP'):
                val=get_data_from_F_normal(mean, sd, query, TINY);
        if(val>vec2[0]):
            val=vec2[0];
    elif (query > vec1[dim-1]): #/* if query falls above maximum value in vector 1 */
        if(lu=='QUAN'):
            if(par=='PRCP'):
                ##### HACK
                if (mean==0):
## KRA Commented out:                     print ("F HACK invoked")
                    val = TINY
                else:
                    val=get_F_from_data_EVI(mean, sd, query, TINY);
            elif(par=='TEMP'):
                val=get_F_from_data_normal(mean, sd, query, TINY);
        if(lu=='DATA'):
            if(par=='PRCP'):
                ##### HACK
                if (mean==0):
                    print ("DATA HACK invoked")
                    val=0
                else:
                    val=get_data_from_F_EVI(mean, sd, query, TINY);
            if(par=='TEMP'):
                val=get_data_from_F_normal(mean, sd, query, TINY);

        if(val<vec2[dim-1]):
            val=vec2[dim-1];
    #/* otherwise, it is within the range of known data in vector 1 
    # do a linear interpolation between known points */
    else:  
        for ndx in range(dim-1): ## Looping from 0 to dim-2
            if(query <= vec1[ndx+1]) and (query >= vec1[ndx]):
                if((vec1[ndx+1]-vec1[ndx]) == 0):
                    A=1;
                else: 
                    A = (vec1[ndx+1]-query)/(vec1[ndx+1]-vec1[ndx]);
                val = A*vec2[ndx]+(1-A)*vec2[ndx+1];
                break;
    if (par=='PRCP'):
        if (val<0):
            val=0
 
    return val;
### The function for bias-correction

#The following function estimates quantile for given data, following normal distribution

def get_F_from_data_normal(mean, sd, x, TINY):
    
    ## uses approximation from Handbook of Hydrology eq. 18.2.2
    sign=1;
    z = (x-mean)/sd;
    
    if(z<0):
        sign=-1; 
        z*=sign
        
    F = 1-0.5*exp(-1*((83*z+351)*z+562)/(165+703/z))
  
    if(F==1.0):
        F-=TINY  # adjustment so that outlier ensembles don't.
        # print ("Invoked TINY in get_F_from_data_normal")

    if(F==0.0):
        F+=TINY
        # print ("Invoked TINY in get_F_from_data_normal")
    
    if(F>=1.0) or (F<=0.0):
        #mostly obviated by above ifs; could be used
        print ("Error in get_F_from_data_normal")
        print ("Bad quantile: data={:0.2f} mean={:0.2f} sd={:0.2f} quant={:0.2f} z={:0.2f}".format(x,mean,sd,F,z))
        sys.exit(1)
    else:
        if(sign==-1):
            F = 1-F
        return F;

## This function converts quantile into value following normal distribution
def get_data_from_F_normal(mean, sd, F, TINY):
    # uses approximation from Handbook of Hydrology eq. 18.2.3a 
    if(F>=1.0) or (F<=0.0):
        print ("Error in get_data_from_F_normal")
        print ("Bad quantile: mean={:0.2f} sd={:0.2f} quant={:0.2f}".format(mean,sd,F))
        sys.exit(0)
    
    ## Calculating z
    z = (pow(F,0.135)-pow(1-F,0.135))/0.1975;
    
    ## Now doing quality control on z
    
    if(z>5):
        print ("Pushing upper limit of applicability of equation in get_data_from_F_normal. Best range is 0<z<5. z= {:0.2f}".format(z))
  
    x = z*sd+mean;
    return x;
## The following function provides quantile for a given value following Gumbel EVI type-1 distribution
def get_F_from_data_EVI(mean, sd, x, TINY):
    #Gumbel (EV Type I distribution) for maxima
    b = 3.14159/(sd*sqrt(6));
    a = mean-0.5772/b;
    F = exp(-1*exp(-b*(x-a)));
    
    ## For outliers in quantile
    if(F==1.0):
        F-=TINY; 
        #print ("Invoked TINY in get_F_from_data_EVI")
    if(F==0.0):
        F+=TINY;     
        #print ("Invoked TINY in get_F_from_data_EVI")
    
    if(F>=1.0) or (F<=0.0):
        print ("Error in get_F_from_data_EVI")
        print ("Bad quantile: data={:0.2f} mean={0.2f} sd={0.2f} quant={0.2f}".format(x,mean,sd,F));
        sys.exit(1);
    else:
        return F;
## The following function provides value for a given quantiles following Gumbel EVI type-1 distribution
def get_data_from_F_EVI(mean, sd, F, TINY):
    #Gumbel (EV Type I distribution) for maxima
    if(F>=1.0) or (F<=0.0):
        print ("Error in get_data_from_F_EVI")
        print ("Bad quantile: mean={0.2f} sd={0.2f} quant={0.2f}".format(mean,sd,F))
        sys.exit(0);
    b = 3.14159/(sd*sqrt(6))
    a = mean-0.5772/b
    x = a-(1/b)*(log(-log(F)));
    return x;
## The following function provide quantile for a given value following Weibul distribution
def get_F_from_data_weibul(mean, sd, skew, x, TINY):
    #Weibull (EV Type III distribution) for minima, bounded at zero 
    #approximation for a (alpha) is eq. 11-32 from Kite (1977)
    
    #if(skew>8.214) or (skew<-1): 
    #    print ("Outside limit for table in get_data_from_F_weibul")
    #    print ("best range is -1<skew<8.2 skew= {:0.2f}".format(skew));
    
    a, A, B = weibul_params(skew);
    b = A*sd+mean;
    bound = b-B*sd;
    #lower bound minimum of zero, but allow larger minima if data say so
    if(bound<0):
        bound=0;
    if(bound > x): 
        bound=0;
    c = (x-bound)/(b-bound)
    if (c<0):
        c=0
        #print ('Invoked C=0 in get_F_from_data_weibul')
    F=1-exp(-1*pow(c,a));
    
    if(F==1.0):
        # For outliers
        F-=TINY;
        #print ("Invoked TINY in get_F_from_data_weibul")
    
    if(F==0.0):
        F+=TINY;
        #print ("Invoked TINY in get_F_from_data_weibul")
    
    if(F>=1.0) or (F<=0.0):
        print ("Error in get_F_from_data_weibul")
        print ("Bad quantile: data={:0.2f} mean={:0.2f} sd={:0.2f} quant={:0.2f}".format(x,mean,sd,F));
        print ("Other info: A={:0.3f} b={:0.3f} B={:0.3f} bound={:0.3f}".format(A,b,B,bound));
        sys.exit(1);
    else: 
        return F;
# The following function provides data for a given quantile following weibull distributions
def get_data_from_F_weibul(mean, sd, skew, F, TINY):
    #/* Weibull (EV Type III distribution) for minima, bounded at zero */
    #approximation for a (alpha) is eq. 11-32 from Kite (1977) */
    if(F>=1.0) or (F<=0.0):
        print ("Error entering get_F_from_data_weibul")
        print ("Bad quantile: mean={:0.2f} sd={:0.2f} skew={:0.2f} quant={:0.2f}".format(mean,sd,skew,F));
        sys.exit(1);
    
    #if(skew>8.2) or (skew<-1):
    #    print ("Outside limit for table in get_data_from_F_weibul")
    #    print ("best range is -1<skew<8.2 skew= {:0.2f}".format(skew));
    
    a, A, B = weibul_params(skew);
    
    b = A*sd+mean;
    bound = b-B*sd;
    
    #/*  lower bound minimum of zero, but allow larger minima if data say so */
    if(bound<0): 
        bound=0;
    
    x = ((pow(-log(1-F),1/a))*(b-bound)) + bound;
    return x;
def weibul_params(skew):
    #returns alpha, Aalpha, and Balpha for the Weibull distrubution
    #table taken from Statistical Methods in Hydrology, Haan, shown below 
    sk = [-1.000, -0.971, -0.917, -0.867, -0.638, -0.254, 0.069, 0.359, 0.631, 0.896,  1.160,  1.430,  1.708,  2.000, 2.309, 2.640, 2.996,  3.382,  3.802,  4.262,  4.767,  5.323, 5.938, 6.619, 7.374,  8.214];
    inva = [0.020, 0.030, 0.040, 0.050, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900, 2.000, 2.100, 2.200];
    Avec = [0.446,  0.444,  0.442,  0.439,  0.425,  0.389,  0.346,  0.297, 0.246,  0.193,  0.142,  0.092,  0.044,  0.000, -0.040, -0.077, 
           -0.109, -0.136, -0.160, -0.180, -0.196, -0.208, -0.217, -0.224, -0.227, -0.229];
    Bvec = [40.005, 26.987, 20.481, 16.576, 8.737, 4.755, 3.370, 2.634, 2.159,  1.815,  1.549,  1.334, 1.154, 1.000, 0.867, 0.752, 
           0.652,  0.563,  0.486,  0.418, 0.359, 0.308, 0.263, 0.224, 0.190, 0.161];
    
    if(skew>sk[-1]):
        skew=sk[-1];
    elif (skew<sk[0]): 
        skew=sk[0];
  
    for ndx in range(len(sk)-1): ## looping from 0 to 24
        if(skew <= sk[ndx+1]) and (skew >= sk[ndx]):
            if((sk[ndx+1]-sk[ndx]) == 0):
                A=1;
            else:
                A = (sk[ndx+1]-skew)/(sk[ndx+1]-sk[ndx]);
            a = 1/(A*inva[ndx]+(1-A)*inva[ndx+1]);
            Aa = A*Avec[ndx]+(1-A)*Avec[ndx+1];
            Ba = A*Bvec[ndx]+(1-A)*Bvec[ndx+1];
            break
    return a, Aa, Ba
