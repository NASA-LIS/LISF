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
import semivar

#------------------------------------------------------------------------------
def usage():
    """Print usage statement to standard out"""
    print("Usage: %s CONFIGFILE PARAMFILE" %(sys.argv[0]))
    print("  CONFIGFILE is config file for this script")
    print("  PARAMFILE is output file storing best-fit parameters")

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

# Get the output files
paramfile = sys.argv[2]

# Read the datafile
distvector, variovector, samplesize = \
    semivar.readdata(vario_filename, max_distance)

# Fit function
# Comparing gages to background field.  Three parameters must be fit.
sigma2_gage_guess = np.amin(variovector)
sigma2_back_guess = np.amax(variovector) - sigma2_gage_guess
L_back_guess = distvector[2]
fit_func = semivar.fit_func_dict[function_type]
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

