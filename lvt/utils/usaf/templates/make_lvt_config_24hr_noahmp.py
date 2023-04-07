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
Sample script to customize lvt.config files for noahmp401 postprocessing for
557WW.
"""

import datetime
import os


_TEMPLATE = "templates/lvt.config.template.noahmp401"

_STARTDT = datetime.datetime(2022, 8, 1, 12)
_ENDDT = datetime.datetime(2022, 8, 2, 12)

#_OUTPUT = "netcdf"
_OUTPUT = "grib2"

# Most variables are processed independently, and are listed below.
_VAR_ATTRIBUTES = {
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
_VAR_ATTRIBUTES_SPECIAL = {
    "Tair_f_min":
    "Tair_f_min  1  1  K      -  1  1 Tair_f_min  1  1  K      -  1  1",
    "RHMin_inst":
        "RHMin       1  1  %      -  0  1 RHMin       1  1  %      -  0  1",
}

# Smooth variables that are perturbed, derived from perturbed variables,
# or are LSM outputs that are affected by perturbed variables via physics.
_SMOOTH_VARS = ["Evap_tavg", "LWdown_f_tavg",
               "SoilMoist_tavg",
               "SoilTemp_tavg", "SWdown_f_tavg",
               "Tair_f_max", "Tair_f_tavg",
               "TotalPrecip_acc", "Tair_f_min", "RHMin_inst"]

def _main():
    """Main Driver"""

    with open(_TEMPLATE, 'r', encoding="ascii") as file:
        lines = file.readlines()

    varlist = list(_VAR_ATTRIBUTES.keys())
    varlist.append("RHMin_inst")  # RHMin will be handled specially below
    varlist.sort()
    for var in varlist:
        newlines = []
        for line in lines:
            if "LVT output format:" in line:
                line = f"LVT output format: {_OUTPUT}\n"
            elif "Process HYCOM data:" in line:
                line = "Process HYCOM data: 0\n"
            elif "Apply noise reduction filter:" in line:
                if var in _SMOOTH_VARS:
                    line = "Apply noise reduction filter: 1\n"
                else:
                    line = "Apply noise reduction filter: 0\n"
            elif "Starting year:" in line:
                line = f"Starting year: {_STARTDT.year}\n"
            elif "Starting month:" in line:
                line = f"Starting month: {_STARTDT.month}\n"
            elif "Starting day:" in line:
                line = f"Starting day: {_STARTDT.day}\n"
            elif "Starting hour:" in line:
                line = f"Starting hour: {_STARTDT.hour}\n"
            elif "Ending year:" in line:
                line = f"Ending year: {_ENDDT.year}\n"
            elif "Ending month:" in line:
                line = f"Ending month: {_ENDDT.month}\n"
            elif "Ending day:" in line:
                line = f"Ending day: {_ENDDT.day}\n"
            elif "Ending hour:" in line:
                line = f"Ending hour: {_ENDDT.hour}\n"
            elif "LVT clock timestep:" in line:
                line = 'LVT clock timestep: "24hr"\n'
            elif "LVT diagnostic file:" in line:
                line = f"LVT diagnostic file: logs/lvtlog.{var}.24hr"
            elif "LVT datastream attributes table::" in line:
                line = "LVT datastream attributes table::\n"
                # Special handling for RHMin_inst, which must be processed with
                # Tair_f_min
                if var == "RHMin_inst":
                    keys = sorted(list(_VAR_ATTRIBUTES_SPECIAL.keys()))
                    for key in keys:
                        line += f"{_VAR_ATTRIBUTES_SPECIAL[key]}\n"
                else:
                    line += f"{_VAR_ATTRIBUTES[var]}\n"
            elif "Metrics attributes file:" in line:
                line = 'Metrics attributes file: "templates/METRICS.TBL"\n'
            elif "Metrics computation frequency:" in line:
                line = 'Metrics computation frequency: "24hr"\n'
            elif "Metrics output directory:" in line:
                line = f"Metrics output directory: OUTPUT/STATS.{var}.24hr\n"
            elif "Metrics output frequency:" in line:
                line = 'Metrics output frequency: "24hr"\n'
            elif "LIS output attributes file:" in line:
                line = "LIS output attributes file:"
                line += f" ./templates/MODEL_OUTPUT_LIST.TBL.lvt_557post.{var}.24hr\n"

            newlines.append(line)

        if not os.path.exists("configs"):
            os.mkdir("configs")
        newfile = f"configs/lvt.config.{var}.24hr"
        print(f"Writing {newfile}")
        with open(newfile, "w", encoding="ascii") as file:
            for line in newlines:
                file.write(line)

if __name__ == "__main__":
    _main()
