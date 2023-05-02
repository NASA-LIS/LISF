#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
# -----------------------------------------------------------------------------
#
# SCRIPT: read_surf.py
#
# PURPOSE:  Reads SURF files (or lookalike files based on LIS), and prints
# out metadata for inspection.
#
# FUTURE:  Add plotting capability with matplotlib.
#
# REQUIREMENTS as of 31 March 2020:
# * Python 3.6
# * UKMO MULE Python library 2020.01.1
# * NumPy Python library (for array objects)
#
# REVISIONS
# 31 Mar 2020: Eric Kemp (SSAI), initial version
# 06 Dec 2022: Eric Kemp (SSAI), changed to improve pylint score.
#
# -----------------------------------------------------------------------------
"""

# Standard modules
import sys

# Non-standard modules
#import matplotlib.pyplot as plt
import mule
#import numpy

#if len(sys.argv) in [2, 3]:
if len(sys.argv) in [2]:
    FILENAME = sys.argv[1]
else:
    print(f"Usage: {sys.argv[0]} filename")
    sys.exit(1)

# plot = False
# if len(sys.argv) == 3:
#     if sys.argv[2] == "--plot":
#         plot = True

ANCIL = mule.AncilFile.from_file(FILENAME)

if ANCIL.fixed_length_header is not None:
    print("***fixed_length_header***")
    for header in ANCIL.fixed_length_header.HEADER_MAPPING:
        key, i = header[0], header[1]
        print(i, key, ' ', ANCIL.fixed_length_header.raw[i])
        # Useful to store the number of level dependent constants for later
        if key == "level_dependent_constants_dim2":
            length_ldc = int(ANCIL.fixed_length_header.raw[i])

# NOTE:  Pylint cannot see most of the components of the AncilFile class,
# and issues error messages and penalizes accordingly.  So, we disable the
# no-member check for sanity.
# pylint: disable=no-member
if ANCIL.integer_constants is not None:
    print("***integer_constants***")
    for header in ANCIL.integer_constants.HEADER_MAPPING:
        key, i = header[0], header[1]
        print(i, key, ' ', ANCIL.integer_constants.raw[i])

if ANCIL.real_constants is not None:
    print("***real_constants***")
    for header in ANCIL.real_constants.HEADER_MAPPING:
        key, i = header[0], header[1]
        print(i, key, ' ', ANCIL.real_constants.raw[i])

if ANCIL.level_dependent_constants is not None:
    print("***level_dependent_constants***")
    # NOTE...MULE 2020.01.1 removed level dependent constants from
    # ancil files, but this is required by the UM RECON preprocessor.
    # So, we use the FieldsFile equivalent here.
    # Since older SURF files may have a subset of the FieldsFile level
    # dependent constants, we have additional logic to check the number
    # of constants in the file.
    for header in mule.ff.FF_LevelDependentConstants.HEADER_MAPPING:
        key, i = header[0], header[1]
        IDX = header[1][1]
        if IDX <= length_ldc:
            print(IDX, key, ' ', ANCIL.level_dependent_constants.raw[i])

if ANCIL.row_dependent_constants is not None:
    print("***row_dependent_constants***")
    for header in ANCIL.level_dependent_constants.HEADER_MAPPING:
        key, i = header[0], header[1]
        print(i, key, ' ', ANCIL.row_dependent_constants.raw[i])

if ANCIL.column_dependent_constants is not None:
    print("***row_dependent_constants***")
    for header in ANCIL.column_dependent_constants.HEADER_MAPPING:
        key, i = header[0], header[1]
        print(i, key, ' ', ANCIL.column_dependent_constants.raw[i])

print('*** ', len(ANCIL.fields), ' fields in file')
for jj, field in enumerate(ANCIL.fields):
    print(jj+1, field)
    for header in field.HEADER_MAPPING:
        key, i = header[0], header[1]
        print(i, key, ' ', field.raw[i])

#------------------------------------------------------------------------------
# FUTURE:  Add plotting support for matplotlib
# if plot:
#     for field in ANCIL.fields:
#         # if field.lbfc not in [93,0,122,23,1510,329,330]:
#         #    continue
#         # print field
#         # for i in range(0,len(field.raw)):
#         #    print i, field.raw[i]
#         data = field.get_data()
#         data = numpy.ma.masked_where(data == -9999, data)
#         data = numpy.ma.masked_where(data < 0, data)

#         plt.pcolormesh(data)
#         plt.axis("tight")
#         if field.lbfc == 16:
#             plt.title("water_temp")
#             plt.clim(270, 310)
#         if field.lbfc == 37:
#             plt.title("aice")
#             plt.clim(0, 1)
#         if field.lbfc == 687:
#             plt.title("hi")
#             plt.clim(0, 2)
#         if field.lbfc == 93:
#             plt.title("SWE")
#             plt.clim(0, 250)
#         if field.lbfc == 0:
#             plt.title("SnowDepth")
#             plt.clim(0, 3)
#         if field.lbfc == 122:
#             plt.title("SoilMoist %s" % (field.blev))
#             plt.clim(0, 1)
#         if field.lbfc == 23:
#             plt.title("SoilTemp %s" % (field.blev))
#             plt.clim(210, 330)
#         if field.lbfc == 1510:
#             plt.title("AvgSurfT")
#             plt.clim(210, 330)
#         if field.lbfc == 329:
#             plt.title("JULES_SM_WILT")
#             plt.clim(0, 1)
#         if field.lbfc == 330:
#             plt.title("JULES_SM_CRIT")
#             plt.clim(0, 1)
#         plt.colorbar()
#         plt.show()
