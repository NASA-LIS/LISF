#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.3
#
# Copyright (c) 2020 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

import datetime
import os
import sys

template = "template/lvt.config.template.jules50"

startdt = datetime.datetime(2007, 12, 1, 0)
enddt = datetime.datetime(2007, 12, 2, 0)

output = "netcdf"
#output = "grib2"

# Most variables are processed independently, and are listed below.
var_attributes = {
    "Evap_tavg":
    "Evap        1  1  kg/m2s -  1  1 Evap        1  1  kg/m2s -  1  1",
    "LWdown_f_tavg":
        "LWdown_f    1  1  W/m2   -  1  1 LWdown_f    1  1  W/m2   -  1  1 ",
    "SoilMoist_tavg":
        "SoilMoist   1  4  m3/m3  -  1  4 SoilMoist   1  4  m3/m3  -  1  4",
    "SoilTemp_tavg":
        "SoilTemp    1  4  K      -  1  4 SoilTemp    1  4  K      -  1  4",
    "SWdown_f_tavg":
        "SWdown_f    1  1  W/m2   -  1  1 SWdown_f    1  1  W/m2   -  1  1",
    "Tair_f_max":
        "Tair_f_max  1  1  K      -  1  1 Tair_f_max  1  1  K      -  1  1",
    "Tair_f_tavg":
        "Tair_f      1  1  K      -  1  1 Tair_f      1  1  K      -  1  1",
    "TotalPrecip_acc":
        "TotalPrecip 1  1  kg/m2  -  3  1 TotalPrecip 1  1  kg/m2  -  3  1",
    "Wind_f_tavg":
        "Wind_f      1  1  m/s    -  1  1 Wind_f      1  1  m/s    -  1  1",
}

# RHMin must be processed with Tair_f_min, so these are listed together
var_attributes_special = {
    "Tair_f_min":
    "Tair_f_min  1  1  K      -  1  1 Tair_f_min  1  1  K      -  1  1",
    "RHMin_inst":
        "RHMin       1  1  %      -  0  1 RHMin       1  1  %      -  0  1",
}

# Smooth variables that are perturbed, derived from perturbed variables,
# or are LSM outputs that are affected by perturbed variables via physics.
smooth_vars = ["Evap_tavg", "LWdown_f_tavg", "SoilMoist_tavg",
               "SoilTemp_tavg", "SWdown_f_tavg",
               "Tair_f_max", "Tair_f_tavg",
               "TotalPrecip_acc", "Tair_f_min", "RHMin_inst"]

lines = open(template, 'r').readlines()

vars = list(var_attributes.keys())
vars.append("RHMin_inst")  # RHMin will be handled specially below
vars.sort()
for var in vars:
    newlines = []
    for line in lines:
        if "LVT output format:" in line:
            line = "LVT output format: %s\n" % (output)
        elif "Process HYCOM data:" in line:
            line = "Process HYCOM data: 0\n"
        elif "Apply noise reduction filter:" in line:
            if var in smooth_vars:
                line = "Apply noise reduction filter: 1\n"
            else:
                line = "Apply noise reduction filter: 0\n"
        elif "Starting year:" in line:
            line = "Starting year: %s\n" % (startdt.year)
        elif "Starting month:" in line:
            line = "Starting month: %s\n" % (startdt.month)
        elif "Starting day:" in line:
            line = "Starting day: %s\n" % (startdt.day)
        elif "Starting hour:" in line:
            line = "Starting hour: %s\n" % (startdt.hour)
        elif "Ending year:" in line:
            line = "Ending year: %s\n" % (enddt.year)
        elif "Ending month:" in line:
            line = "Ending month: %s\n" % (enddt.month)
        elif "Ending day:" in line:
            line = "Ending day: %s\n" % (enddt.day)
        elif "Ending hour:" in line:
            line = "Ending hour: %s\n" % (enddt.hour)
        elif "LVT clock timestep:" in line:
            line = 'LVT clock timestep: "24hr"\n'
        elif "LVT diagnostic file:" in line:
            line = "LVT diagnostic file: logs/lvtlog.%s.24hr" % (var)
        elif "LVT datastream attributes table::" in line:
            line = "LVT datastream attributes table::\n"
            # Special handling for RHMin_inst, which must be processed with
            # Tair_f_min
            if var == "RHMin_inst":
                keys = sorted(list(var_attributes_special.keys()))
                for key in keys:
                    line += "%s\n" % (var_attributes_special[key])
            # The general case
            else:
                line += "%s\n" % (var_attributes[var])
        elif "Metrics computation frequency:" in line:
            line = 'Metrics computation frequency: "24hr"\n'
        elif "Metrics output directory:" in line:
            line = "Metrics output directory: OUTPUT/STATS.%s.24hr\n" % (var)
        elif "Metrics output frequency:" in line:
            line = 'Metrics output frequency: "24hr"\n'
        elif "LIS output attributes file:" in line:
            line = "LIS output attributes file:"
            line += " ./tables/MODEL_OUTPUT_LIST.TBL.lvt_557post.%s.24hr\n" % (var)

        newlines.append(line)

    newfile = "configs/lvt.config.%s.24hr" % (var)
    print("Writing %s" % (newfile))
    f = open(newfile, "w")
    for line in newlines:
        f.write(line)
    f.close()
