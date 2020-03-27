#!/usr/bin/env python
# -----------------------------------------------------------------------------
#
# SCRIPT: read_surf.py
#
# PURPOSE:  Reads SURF files (or lookalike files based on LIS), and prints
# out metadata for inspection.
#
# FUTURE:  Add plotting capability with matplotlib.
#
# REQUIREMENTS as of 27 March 2020:
# * Python 3.6
# * UKMO MULE Python library 2020.01.1
# * NumPy Python library (for array objects)
#
# -----------------------------------------------------------------------------
# Standard modules
import os
import sys

# Non-standard modules
#import matplotlib.pyplot as plt
import mule
import numpy

#if len(sys.argv) in [2, 3]:
if len(sys.argv) in [2]:
    filename = sys.argv[1]
else:
    print("Usage: %s filename" %(sys.argv[0]))
    sys.exit(1)

# plot = False
# if len(sys.argv) == 3:
#     if sys.argv[2] == "-plot":
#         plot = True

ancil = mule.AncilFile.from_file(filename)

if ancil.fixed_length_header != None:
    print("***fixed_length_header***")
    for header in mule._UM_FIXED_LENGTH_HEADER:
        key, i = header[0], header[1]
        print(i, key, ' ', ancil.fixed_length_header.raw[i])
        # Useful to store the number of level dependent constants for later
        if key == "level_dependent_constants_dim2":
            length_ldc = int(ancil.fixed_length_header.raw[i])

if ancil.integer_constants != None:
    print("***integer_constants***")
    for header in mule.ancil._ANCIL_INTEGER_CONSTANTS:
        key, i = header[0], header[1]
        print(i, key, ' ', ancil.integer_constants.raw[i])

if ancil.real_constants != None:
    print("***real_constants***")
    for header in mule.ancil._ANCIL_REAL_CONSTANTS:
        key, i = header[0], header[1]
        print(i, key, ' ', ancil.real_constants.raw[i])

if ancil.level_dependent_constants != None:
    print("***level_dependent_constants***")
    # NOTE...MULE 2020.01.1 removed level dependent constants from
    # ancil files, but this is required by the UM RECON preprocessor.
    # So, we use the FieldsFile equivalent here.
    # Since older SURF files may have a subset of the FieldsFile level
    # dependent constants, we have additional logic to check the number
    # of constants in the file.
    for header in mule.ff._FF_LEVEL_DEPENDENT_CONSTANTS:
        key, i = header[0], header[1]
        idx = header[1][1]
        if idx <= length_ldc:
            print(idx, key, ' ', ancil.level_dependent_constants.raw[i])

if ancil.row_dependent_constants != None:
    print("***row_dependent_constants***")
    for header in mule.ancil._ANCIL_ROW_DEPENDENT_CONSTANTS:
        key, i = header[0], header[1]
        print(i, key, ' ', ancil.row_dependent_constants.raw[i])

if ancil.column_dependent_constants != None:
    print("***row_dependent_constants***")
    for header in mule.ancil._ANCIL_COLUMN_DEPENDENT_CONSTANTS:
        key, i = header[0], header[1]
        print(i, key, ' ', ancil.column_dependent_constants.raw[i])

print('*** ', len(ancil.fields), ' fields in file')
counter = 0
for field in ancil.fields:
    counter += 1
    print(counter, field)
    for header in mule._LOOKUP_HEADER_3:
        key, i = header[0], header[1]
        print(i, key, ' ', field.raw[i])

# if plot:
#     for field in ancil.fields:
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
