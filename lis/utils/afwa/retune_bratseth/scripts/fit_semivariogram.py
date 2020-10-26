#!/usr/bin/env python

# Standard library
import configparser
import math
import os
import sys

# Other libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#------------------------------------------------------------------------------
def usage():
    print("Usage: %s config.cfg" %(sys.argv[0]))

#------------------------------------------------------------------------------
def fit_func_gaussian(x, a, b, c):
    if a < 0:
        return -9999
    if b < 0:
        return -9999
    if c < 30:
        return -9999
    # Here a is sigma2_o, b is sigma2_b, and c is L_b
    return a  + b*(1. - np.exp(-1*x*x/c/c))

#------------------------------------------------------------------------------
def fit_func_soar(x, a, b, c):
    if a < 0:
        return -9999
    if b < 0:
        return -9999
    if c < 0:
        return -9999
    # Here a is sigma2_o, b is sigma2_b, and c is L_b
    return a  + b*(1. - ((1. + x/c)*np.exp(-1*x/c)))

#------------------------------------------------------------------------------
def fit_func_invexp(x, a, b, c):
    if a < 0:
        return -9999
    if b < 0:
        return -9999
    if c < 0:
        return -9999
    # Here a is sigma2_o, b is sigma2_b, and c is L_b
    return a  + b*(1. - np.exp(-1*x/c))

#------------------------------------------------------------------------------
fit_func_dict = {
    "Gaussian" : fit_func_gaussian,
    "InvExp"   : fit_func_invexp,
    "SOAR"     : fit_func_soar,
}

#------------------------------------------------------------------------------
def readdata(filename, maxdist):
    distvector = []
    variovector = []
    countvector = []
    lines = open(filename,"r").readlines()
    samplesize = 0
    for line in lines:
        if "#" in line:
            continue
        dist = float(line.split()[1])
        vario = float(line.split()[3])
        count = int(line.split()[5])

        if dist == 0:
            continue
        if dist > maxdist:
            continue

        samplesize += count

        distvector.append(dist)
        variovector.append(vario)
        countvector.append(count)

    # Convert to numpy arrays
    distvector = np.array(distvector)
    variovector = np.array(variovector)
    return distvector, variovector, samplesize

#------------------------------------------------------------------------------
# Check command line
if len(sys.argv) != 2:
    print("[ERR] Bad command line arguments!")
    usage()
    sys.exit(1)

# Read config file
cfgfile = sys.argv[1]
if not os.path.exists(cfgfile):
    print("[ERR] Config file %s does not exist!" %(cfgfile))
    sys.exit(1)
config = configparser.ConfigParser()
config.read(cfgfile)

# Handle the Input section
try:
    vario_filename = config.get('Input', 'vario_filename')
    max_distance   = config.getfloat('Input', 'max_distance')
except:
    print("[ERR] Problem reading from config file!")
    raise

if not os.path.exists(vario_filename):
    print("[ERR] %s does not exist!" %(vario_filename))
    sys.exit(1)
if not max_distance > 0:
    print("[ERR] Maximum distance must be positive!")
    sys.exit(1)

# Handle the Fit section
try:
    function_type  = config.get('Fit', 'function_type')
except:
    print("[ERR] Problem reading from config file!")
    raise

function_types = fit_func_dict.keys()
function_types = list(function_types)
function_types.sort()
if function_type not in function_types:
    print('[ERR] function type %s is not supported!' %(function_type))
    print("Currently only the following functions can be fit:")
    for type in function_types:
        print("  %s" %(type))
    sys.exit(1)

# Handle the Plot section
try:
    title  = config.get('Plot', 'title')
    xlabel = config.get('Plot', 'xlabel')
    ylabel = config.get('Plot', 'ylabel')
    oblabel = config.get('Plot', 'oblabel')
    bglabel = config.get('Plot', 'bglabel')
except:
    print("[ERR] Problem reading from config file!")
    raise

# Read the datafile
distvector, variovector, samplesize = \
    readdata(vario_filename, max_distance)

# Fit function
# Comparing gages to background field.  Three parameters must be fit.
sigma2_gage_guess = np.amin(variovector)
sigma2_back_guess = np.amax(variovector) - sigma2_gage_guess
L_back_guess = distvector[2]
fit_func = fit_func_dict[function_type]
sigma2_gage_min = 0.1*sigma2_gage_guess
sigma2_gage_max = np.amax(variovector)
sigma2_back_min = 0.1*sigma2_gage_guess
sigma2_back_max = np.amax(variovector)
L_back_min = L_back_guess
L_back_max = distvector[-1]
popt, pconv = \
             curve_fit(fit_func, distvector, variovector,
                       p0=(sigma2_gage_guess, sigma2_back_guess, L_back_guess),
                       bounds=([sigma2_gage_min, sigma2_back_min, L_back_min],
                               [sigma2_gage_max, sigma2_back_max, L_back_max]))

sigma2_gage = popt[0]
sigma2_back = popt[1]
L_back      = popt[2]

# Plot the semivariogram
distvector_tmp = np.array([0])
distvector     = np.concatenate((distvector_tmp, distvector))
print(distvector)

variovector_tmp = np.array([np.nan])
variovector = np.concatenate((variovector_tmp, variovector))
print(variovector)

plt.plot(distvector, variovector, "b+",
         distvector, fit_func(distvector, *popt), "r")

# Annotate
fulltitle  = "%s\n"%(title)
fulltitle += "Based on %s comparisons of innovations\n" %(samplesize)
plt.title(fulltitle)
plt.xlabel(r"%s" %(xlabel))
plt.ylabel(r"%s"%(ylabel))
plt.legend(["Data", "%s" %(function_type)],
           loc='lower right')

params = r"$\sigma_{%s}^2 = %f, \sigma_{%s}^2=%f, L_{%s} = %f$" \
         %(oblabel, sigma2_gage, bglabel, sigma2_back, bglabel, L_back)
plt.figtext(0.2, 0.9, params)
plt.grid(True)
plt.show()
