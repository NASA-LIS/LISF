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
import semivar

#------------------------------------------------------------------------------
def usage():
    """Print usage statement to standard out."""
    print("Usage: %s CONFIGFILE PARAMFILE" %(sys.argv[0]))
    print("  CONFIG is config file for this script")
    print("  PARAMFILE contains best-fit parameters from fit_semivariogram.py")

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

vario_filename, max_distance = semivar.read_input_section_cfg(config)
function_type = semivar.read_fit_section_cfg(config)
title, xlabel, ylabel, oblabel, bglabel = \
    semivar.read_plot_section_cfg(config)

# Get the param file
paramfile = sys.argv[2]
if not os.path.exists(paramfile):
    print("[ERR] Paramfile %s does not exist!" %(paramfile))
    sys.exit(1)

# Read the datafile
distvector, variovector, samplesize = \
    semivar.readdata(vario_filename, max_distance)

# Read the paramfile
sigma2_gage, sigma2_back, L_back = read_param_file(paramfile)
popt = [sigma2_gage, sigma2_back, L_back]

# Plot the semivariogram
distvector_tmp = np.array([0])
distvector     = np.concatenate((distvector_tmp, distvector))

variovector_tmp = np.array([np.nan])
variovector = np.concatenate((variovector_tmp, variovector))

fit_func = semivar.fit_func_dict[function_type]
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
