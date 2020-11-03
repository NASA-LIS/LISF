#!/usr/bin/env python3
"""
SCRIPT: fit_semivariogram.py

Script for creating a fit to an empirical semivariogram.  Input data are
from the procOBA_NWP or procOBA_Sat programs.  Blacklist file is from
create_blacklist.py.

REVISION HISTORY:
26 Oct 2020: Eric Kemp. Initial Specification.
03 Nov 2020: Eric Kemp. Removed plotting.  Fitted parameters are now saved
               to file.
"""

# Standard library
import configparser
import os
import sys

# Other libraries
import numpy as np
from scipy.optimize import curve_fit

#------------------------------------------------------------------------------
def usage():
    """Print usage statement to standard out"""
    print("Usage: %s CONFIGFILE PARAMFILE" %(sys.argv[0]))
    print("  CONFIGFILE is config file for this script")
    print("  PARAMFILE is output file storing best-fit parameters")

#------------------------------------------------------------------------------
# NOTE: Pylint complains about the single-character variable names not
# conforming to snake_case convention.  For sanity, we disable this test
# here.
# pylint: disable=invalid-name
def fit_func_gaussian(x, a, b, c):
    """Fits a Gaussian function to the semivariogram."""
    if a < 0:
        return -9999
    if b < 0:
        return -9999
    if c < 30:
        return -9999
    # Here a is sigma2_o, b is sigma2_b, and c is L_b
    return a  + b*(1. - np.exp(-1*x*x/c/c))
# pylint: enable=invalid-name

#------------------------------------------------------------------------------
# NOTE: Pylint complains about the single-character variable names not
# conforming to snake_case convention.  For sanity, we disable this test
# here.
# pylint: disable=invalid-name
def fit_func_soar(x, a, b, c):
    """Fits a second-order auto-regressive function to the semivariogram."""
    if a < 0:
        return -9999
    if b < 0:
        return -9999
    if c < 0:
        return -9999
    # Here a is sigma2_o, b is sigma2_b, and c is L_b
    return a  + b*(1. - ((1. + x/c)*np.exp(-1*x/c)))
# pylint: enable=invalid-name

#------------------------------------------------------------------------------
# NOTE: Pylint complains about the single-character variable names not
# conforming to snake_case convention.  For sanity, we disable this test
# here.
# pylint: disable=invalid-name
def fit_func_invexp(x, a, b, c):
    """Fits an inverse exponential function to the semivariogram."""
    if a < 0:
        return -9999
    if b < 0:
        return -9999
    if c < 0:
        return -9999
    # Here a is sigma2_o, b is sigma2_b, and c is L_b
    return a  + b*(1. - np.exp(-1*x/c))
# pylint: enable=invalid-name

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


# Get the output files
paramfile = sys.argv[2]

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

# NOTE: Pylint gives a false-positive warning about unbalanced tuple unpacking
# for the values returned from curve_fit. For sanity, we disable the test here.
# pylint: disable=unbalanced-tuple-unpacking
popt, pconv = \
             curve_fit(fit_func, distvector, variovector,
                       p0=(sigma2_gage_guess, sigma2_back_guess, L_back_guess),
                       bounds=([sigma2_gage_min, sigma2_back_min, L_back_min],
                               [sigma2_gage_max, sigma2_back_max, L_back_max]))
# pylint: enable=unbalanced-tuple-unpacking

sigma2_gage = popt[0]
sigma2_back = popt[1]
L_back      = popt[2]

# Write this to the param file
fd = open(paramfile, "w")
fd.write("SIGMA2_obs: %s\n" %(sigma2_gage))
fd.write("SIGMA2_back: %s\n" %(sigma2_back))
fd.write("L_back: %s\n" %(L_back))
fd.close()

