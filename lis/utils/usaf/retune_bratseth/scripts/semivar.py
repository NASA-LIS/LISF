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
SCRIPT: semivar.py

Contains common functions used by both fit_semivariogram.py and
plot_semivariogram.py.

REVISION HISTORY:
20 Nov 2020: Eric Kemp. Initial specification.
"""

# Standard library
import os
import sys

# Other libraries
import numpy as np

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
    with open(filename,"r",encoding="ascii") as file:
        lines = file.readlines()
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
def read_input_section_cfg(config):
    """Reads the Input section of the config file."""
    try:
        vario_filename = config.get('Input', 'vario_filename')
        max_distance   = config.getfloat('Input', 'max_distance')
    except:
        print("[ERR] Problem reading from config file!")
        raise
    if not os.path.exists(vario_filename):
        print(f"[ERR] {vario_filename} does not exist!")
        sys.exit(1)
    if max_distance <= 0:
        print("[ERR] Maximum distance must be positive!")
        sys.exit(1)
    return vario_filename, max_distance

#------------------------------------------------------------------------------
def read_fit_section_cfg(config):
    """Reads the Fit section of the config file."""
    try:
        function_type  = config.get('Fit', 'function_type')
    except:
        print("[ERR] Problem reading from config file!")
        raise
    function_types = fit_func_dict.keys()
    function_types = list(function_types)
    function_types.sort()
    if function_type not in function_types:
        print(f'[ERR] function type {function_type} is not supported!')
        print("Currently only the following functions can be fit:")
        for func in function_types:
            print(f"  {func}")
            sys.exit(1)
    return function_type

#------------------------------------------------------------------------------
def read_plot_section_cfg(config):
    """Reads the Plot section of the config file."""
    try:
        title  = config.get('Plot', 'title')
        xlabel = config.get('Plot', 'xlabel')
        ylabel = config.get('Plot', 'ylabel')
        oblabel = config.get('Plot', 'oblabel')
        bglabel = config.get('Plot', 'bglabel')
    except:
        print("[ERR] Problem reading from config file!")
        raise
    return title, xlabel, ylabel, oblabel, bglabel
