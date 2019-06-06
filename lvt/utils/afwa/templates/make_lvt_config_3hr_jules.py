#!/usr/bin/env python

import datetime
import os
import sys

template = "lvt.config.template"

#startdt = datetime.datetime(2008,5,15,0)
#enddt = datetime.datetime(2008,5,16,0)
#startdt = datetime.datetime(2018, 2,28, 0)
#enddt =   datetime.datetime(2018, 3, 1, 0)
#startdt = datetime.datetime(2011, 2, 28, 21)
#enddt =   datetime.datetime(2011, 3,  1, 21)
#startdt = datetime.datetime(2011, 8, 31, 21)
#enddt =   datetime.datetime(2011, 9,  1, 21)
startdt = datetime.datetime(2018,  7, 29,  0)
enddt =   datetime.datetime(2018,  7, 31,  0)

#output = "netcdf"
output = "grib2"

# Most variables are processed independently, and are listed below.
var_attributes = {
    "AvgSurfT_inst" : \
        "AvgSurfT    1  1  K      -  0  1 AvgSurfT    1  1  K      -  0  1",
    "AvgSurfT_tavg" : \
        "AvgSurfT    1  1  K      -  1  1 AvgSurfT    1  1  K      -  1  1",
    "CanopInt_inst" : \
        "CanopInt    1  1  kg/m2  -  0  1 CanopInt    1  1  kg/m2  -  0  1",
    "Elevation_inst" : \
        "Elevation   1  1  m      -  0  1 Elevation   1  1  m      -  0  1",
    "Evap_tavg" : \
        "Evap        1  1  kg/m2s -  1  1 Evap        1  1  kg/m2s -  1  1",
    "Greenness_inst" : \
        "Greenness   1  1  -      -  0  1 Greenness   1  1  -      -  0  1",
    "Landmask_inst" : \
        "Landmask    1  1  -      -  0  1 Landmask    1  1  -      -  0  1",
    "Landcover_inst" : \
        "Landcover   1  1  -      -  0  1 Landcover   1  1  -      -  0  1",
    "LWdown_f_inst" : \
        "LWdown_f    1  1  W/m2   -  0  1 LWdown_f    1  1  W/m2   -  0  1 ",
    "LWdown_f_tavg" : \
        "LWdown_f    1  1  W/m2   -  1  1 LWdown_f    1  1  W/m2   -  1  1 ",
    "Psurf_f_inst" : \
        "Psurf_f     1  1  Pa     -  0  1 Psurf_f     1  1  Pa     -  0  1",
    "Psurf_f_tavg" : \
        "Psurf_f     1  1  Pa     -  1  1 Psurf_f     1  1  Pa     -  1  1",
    "Qair_f_inst" : \
        "Qair_f      1  1  kg/kg  -  0  1 Qair_f      1  1  kg/kg  -  0  1",
    "Qair_f_tavg" : \
        "Qair_f      1  1  kg/kg  -  1  1 Qair_f      1  1  kg/kg  -  1  1",
    "Qh_tavg" : \
        "Qh          1  1  W/m2   -  1  1 Qh          1  1  W/m2   -  1  1",
    "Qle_tavg" : \
        "Qle         1  1  W/m2   -  1  1 Qle         1  1  W/m2   -  1  1",
    "Qs_acc" : \
        "Qs          1  1  kg/m2  -  3  1 Qs          1  1  kg/m2  -  3  1",
    "Qsb_acc" : \
        "Qsb         1  1  kg/m2  -  3  1 Qsb         1  1  kg/m2  -  3  1",
    "RelSMC_inst" : \
        "RelSMC      1  4  m3/m3  -  0  4 RelSMC      1  4  m3/m3  -  0  4",
    "SnowDepth_inst" : \
        "SnowDepth   1  1  m      -  0  1 SnowDepth   1  1  m      -  0  1",
    "SoilMoist_inst" : \
        "SoilMoist   1  4  m3/m3  -  0  4 SoilMoist   1  4  m3/m3  -  0  4",
    "SoilMoist_tavg" : \
        "SoilMoist   1  4  m3/m3  -  1  4 SoilMoist   1  4  m3/m3  -  1  4",
    "SoilTemp_inst" : \
        "SoilTemp    1  4  K      -  0  4 SoilTemp    1  4  K      -  0  4",
    "SoilTemp_tavg" : \
        "SoilTemp    1  4  K      -  1  4 SoilTemp    1  4  K      -  1  4",
    "Soiltype_inst" : \
        "Soiltype    1  1  -      -  0  1 Soiltype    1  1  -      -  0  1",
    "SWdown_f_inst" : \
        "SWdown_f    1  1  W/m2   -  0  1 SWdown_f    1  1  W/m2   -  0  1",
    "SWdown_f_tavg" : \
        "SWdown_f    1  1  W/m2   -  1  1 SWdown_f    1  1  W/m2   -  1  1",
    "SWE_inst" : \
        "SWE         1  1  kg/m2  -  0  1 SWE         1  1  kg/m2  -  0  1",
    "Tair_f_inst" : \
        "Tair_f      1  1  K      -  0  1 Tair_f      1  1  K      -  0  1",
    "Tair_f_max" : \
        "Tair_f_max  1  1  K      -  1  1 Tair_f_max  1  1  K      -  1  1",
    "Tair_f_tavg" : \
        "Tair_f      1  1  K      -  1  1 Tair_f      1  1  K      -  1  1",
    "TotalPrecip_acc" : \
        "TotalPrecip 1  1  kg/m2  -  3  1 TotalPrecip 1  1  kg/m2  -  3  1",
    "Wind_f_inst" : \
        "Wind_f      1  1  m/s    -  0  1 Wind_f      1  1  m/s    -  0  1",
    "Wind_f_tavg" : \
        "Wind_f      1  1  m/s    -  1  1 Wind_f      1  1  m/s    -  1  1",
}

# RHMin must be processed with Tair_f_min, so these are listed together
var_attributes_special = {
    "Tair_f_min" : \
        "Tair_f_min  1  1  K      -  1  1 Tair_f_min  1  1  K      -  1  1",
    "RHMin_inst" : \
        "RHMin       1  1  %      -  0  1 RHMin       1  1  %      -  0  1",
}

lines = open(template,'r').readlines()

vars = var_attributes.keys()
vars.append("RHMin_inst") # RHMin will be handled specially below
vars.sort()
firstVar = True
for var in vars:
    newlines = []
    for line in lines:
        if "LVT output format:" in line:
            line = "LVT output format: %s\n" %(output)
        elif "Process HYCOM data:" in line:
            if firstVar:
                line = "Process HYCOM data: 1\n"
            else:
                line = "Process HYCOM data: 0\n"
        elif "Starting year:" in line:
            line = "Starting year: %s\n" %(startdt.year)
        elif "Starting month:" in line:
            line = "Starting month: %s\n" %(startdt.month)
        elif "Starting day:" in line:
            line = "Starting day: %s\n" %(startdt.day)
        elif "Starting hour:" in line:
            line = "Starting hour: %s\n" %(startdt.hour)
        elif "Ending year:" in line:
            line = "Ending year: %s\n" %(enddt.year)
        elif "Ending month:" in line:
            line = "Ending month: %s\n" %(enddt.month)
        elif "Ending day:" in line:
            line = "Ending day: %s\n" %(enddt.day)
        elif "Ending hour:" in line:
            line = "Ending hour: %s\n" %(enddt.hour)
        elif "LVT diagnostic file:" in line:
            line = "LVT diagnostic file: lvtlog.%s.3hr" %(var)
        elif "LVT datastream attributes table::" in line:
            line = "LVT datastream attributes table::\n"
            # Special handling for RHMin_inst, which must be processed with
            # Tair_f_min
            if var == "RHMin_inst": 
                keys = var_attributes_special.keys()
                keys.sort()
                for key in keys:
                    line += "%s\n" %(var_attributes_special[key])
            # The general case
            else:
                line += "%s\n" %(var_attributes[var])
        elif "Metrics output directory:" in line:
            line = "Metrics output directory: STATS.%s.3hr\n" %(var)
        elif "LIS output attributes file:" in line:
            line = "LIS output attributes file:"
            line += " ./MODEL_OUTPUT_LIST.TBL.lvt_557post.%s.3hr\n" %(var)

        newlines.append(line)

    firstVar = False
    newfile = "lvt.config.%s.3hr" %(var)
    print "Writing %s" %(newfile)
    f = open(newfile,"w")
    for line in newlines:
        f.write(line)
    f.close()
