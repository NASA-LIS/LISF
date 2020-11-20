#!/usr/bin/env python3
"""
SCRIPT: plot_semivariogram.py

Script for plotting empirical and fitted semivariograms based on data from
procOBA_NWP or procOBA_Sat, plus fit_semivariogram.py.  For interactive
use.

REQUIREMENTS:
* Python 3
* Matplotlib
* Numpy

REVISION HISTORY:
20 Nov 2020: Eric Kemp. Initial specification.
"""

# Standard library
import configparser
import os
import sys

# Other libraries
import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------------
def usage():
    """Print usage statement to standard out."""
    print("Usage: %s CONFIGFILE PARAMFILE" %(sys.argv[0]))
    print("  CONFIG is config file for this script")
    print("  PARAMFILE contains best-fit parameters from fit_semivariogram.py")

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
    """Reads semivariogram data from file, and returns in lists."""
    dist_vector = []
    vario_vector = []
    count_vector = []
    lines = open(filename,"r").readlines()
    sample_size = 0
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

        sample_size += count

        dist_vector.append(dist)
        vario_vector.append(vario)
        count_vector.append(count)

    # Convert to numpy arrays
    dist_vector = np.array(dist_vector)
    vario_vector = np.array(vario_vector)
    return dist_vector, vario_vector, sample_size

#------------------------------------------------------------------------------
def read_param_file(paramfile):
    """Reads fitted parameters for semivariogram for plotting."""
    lines = open(paramfile, "r").readlines()
    sigma2_gage = None
    sigma2_back = None
    L_back = None
    for line in lines:
        key, value = line.split(":")
        value = float(value)
        if key == "SIGMA2_obs":
            sigma2_gage = value
        elif key == "SIGMA2_back":
            sigma2_back = value
        elif key == "L_back":
            L_back = value
    return sigma2_gage, sigma2_back, L_back

#------------------------------------------------------------------------------
# Check command line
if len(sys.argv) != 3:
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
if max_distance <= 0:
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
    for function_type in function_types:
        print("  %s" %(function_type))
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

# Get the param file
paramfile = sys.argv[2]
if not os.path.exists(paramfile):
    print("[ERR] Paramfile %s does not exist!" %(paramfile))
    sys.exit(1)

# Read the datafile
distvector, variovector, samplesize = \
    readdata(vario_filename, max_distance)

# Read the paramfile
sigma2_gage, sigma2_back, L_back = read_param_file(paramfile)
popt = [sigma2_gage, sigma2_back, L_back]

# Plot the semivariogram
distvector_tmp = np.array([0])
distvector     = np.concatenate((distvector_tmp, distvector))

variovector_tmp = np.array([np.nan])
variovector = np.concatenate((variovector_tmp, variovector))

fit_func = fit_func_dict[function_type]
plt.plot(distvector, variovector, "b+",
         distvector, fit_func(distvector, *popt), "r")

# Annotate
fulltitle  = "%s\n"%(title)
fulltitle += "Based on %s comparisons of innovations\n" %(samplesize)
plt.title(fulltitle)
plt.xlabel(r"%s" %(xlabel))
plt.ylabel(r"%s"%(ylabel))
plt.legend(["Data", "%s Best Fit" %(function_type)],
           loc='lower right')

params = r"$\sigma_{%s}^2 = %f, \sigma_{%s}^2=%f, L_{%s} = %f$" \
         %(oblabel, sigma2_gage, bglabel, sigma2_back, bglabel, L_back)
plt.figtext(0.2, 0.9, params)
plt.grid(True)
plt.show()
