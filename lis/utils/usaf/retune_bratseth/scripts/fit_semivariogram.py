#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
SCRIPT: fit_semivariogram.py

Script for creating a fit to an empirical semivariogram.  Input data are
from the procOBA_NWP or procOBA_Sat programs.  Blacklist file is from
create_blacklist.py.

REVISION HISTORY:
26 Oct 2020: Eric Kemp. Initial Specification.
03 Nov 2020: Eric Kemp. Removed plotting.  Fitted parameters are now saved
               to file.
13 Dec 2021: Eric Kemp. Changed first guess for covariance settings.
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
    print(f"Usage: {sys.argv[0]} CONFIGFILE PARAMFILE")
    print("  CONFIGFILE is config file for this script")
    print("  PARAMFILE is output file storing best-fit parameters")

#------------------------------------------------------------------------------
# Check command line
if len(sys.argv) != 3:
    print("[ERR] Bad command line arguments!")
    usage()
    sys.exit(1)

# Read config file
CFGFILE = sys.argv[1]
if not os.path.exists(CFGFILE):
    print(f"[ERR] Config file {CFGFILE} does not exist!")
    sys.exit(1)
config = configparser.ConfigParser()
config.read(CFGFILE)

vario_filename, max_distance = semivar.read_input_section_cfg(config)
function_type = semivar.read_fit_section_cfg(config)

# Get the output files
paramfile = sys.argv[2]

# Read the datafile
distvector, variovector, samplesize = \
    semivar.readdata(vario_filename, max_distance)

# Fit function
# Comparing gages to background field.  Three parameters must be fit.
sigma2_gage_guess = (np.amin(variovector) + np.amax(variovector)) / 2.
sigma2_back_guess = sigma2_gage_guess
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
with open(paramfile, "w", encoding="ascii") as fd:
    fd.write(f"SIGMA2_obs: {sigma2_gage}\n")
    fd.write(f"SIGMA2_back: {sigma2_back}\n")
    fd.write(f"L_back: {L_back}\n")
