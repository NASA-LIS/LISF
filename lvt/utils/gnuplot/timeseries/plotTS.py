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

# currently assumes eps files - need to modify for ESMPy to read config
# options

import sys

print("Generating the gnuplot script for time series plots ...")

infile = open('ts.config', 'r')
for line in infile:
    try:
        (option, out_img_format) = line.split("Output format for images:")
    except ValueError as e:
        print(e)
    try:
        (option, lvt_config) = line.split("LVT config file:")
    except ValueError as e:
        print(e)
    try:
        (option, lvt_ts_file) = line.split("LVT time series locations file:")
    except ValueError as e:
        print(e)
    try:
        (option, outfilename) = line.split("gnuplot output filename:")
    except ValueError as e:
        print(e)

    try:
        (option, plt_summ_stats) = line.split(
            "Plot summary stats in the image:")
    except ValueError as e:
        print(e)
    try:
        (option, n_plt_vars) = line.split("Number of variables to plot:")
    except ValueError as e:
        print(e)
    try:
        (option, d1_label) = line.split("data series 1 label:")
    except ValueError as e:
        print(e)
    try:
        (option, d1_style) = line.split("data series 1 style:")
    except ValueError as e:
        print(e)
    try:
        (option, d2_label) = line.split("data series 2 label:")
    except ValueError as e:
        print(e)
    try:
        (option, d2_style) = line.split("data series 2 style:")
    except ValueError as e:
        print(e)
    try:
        (option, n_x_tics) = line.split("Number of xtics:")
    except ValueError as e:
        print(e)
    try:
        (option, vnames) = line.split("Variable names:")
    except ValueError as e:
        print(e)
    try:
        (option, vunits) = line.split("Variable units:")
    except ValueError as e:
        print(e)


# open lvt.config file and read entries.
lvtfile = open(lvt_config.strip(), 'r')
for line in lvtfile:
    try:
        (option, statsodir) = line.split("Stats output directory:")
    except ValueError as e:
        print(e)

outfile = open(outfilename.strip(), 'w')
#outfile = open('ts.plt','w')

if out_img_format.strip() == "gif":
    cline = 'set terminal gif large font "Times-Roman" 24\n'
    outfile.writelines(cline)
elif out_img_format.strip() == "eps":
    cline = 'set terminal postscript enhanced eps color "Times-Roman" 24\n'
    outfile.writelines(cline)
else:
    sys.exit(-1)

# 'set output "'//trim(stnname(i))//'_'//&
#                      trim(varname(j))//'_ts.eps"'
cline = 'set size 2,1\n'
outfile.writelines(cline)

cline = 'set datafile missing "-0.999900E+04"\n'
outfile.writelines(cline)

# cline= 'set title "'//trim(stnname(i))//'-'//&
# trim(varname(j))//'('//trim(vunits(j))//')"'

cline = 'set xdata time\n'
outfile.writelines(cline)
cline = 'set timefmt "%Y %m %d %H %M"\n'
outfile.writelines(cline)
cline = 'set format x "%Y/%m"\n'
outfile.writelines(cline)

plt_summ_stats = int(plt_summ_stats)
if plt_summ_stats == 1:
    rmsefile = statsodir.strip() + '/RMSE_SUMMARY_STATS.dat'
    print(rmsefile)
outfile.close()

print('Sucessfully generated the gnuplot script '+outfilename)
