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
Sample script to customize lvt.config files for noah39 postprocessing for
557WW.  Intended for use for producing STR 2d ensemble gridspace files.
"""

import datetime
import os

_TEMPLATE = "templates/lvt.config.template.noah39.str"

_STARTDT = datetime.datetime(2025, 1, 20, 12)
_ENDDT = datetime.datetime(2025, 1, 21, 0)
#_STARTDT = datetime.datetime(2025, 1, 21, 0)
#_ENDDT = datetime.datetime(2025, 1, 21, 12)
#_STARTDT = datetime.datetime(2025, 1, 21, 6)
#_ENDDT = datetime.datetime(2025, 1, 21, 18)
#_STARTDT = datetime.datetime(2025, 1, 21, 12)
#_ENDDT = datetime.datetime(2025, 1, 22, 0)
#_STARTDT = datetime.datetime(2025, 1, 21, 18)
#_ENDDT = datetime.datetime(2025, 1, 22, 6)
#_STARTDT = datetime.datetime(2025, 1, 22, 0)
#_ENDDT = datetime.datetime(2025, 1, 22, 12)

_OUTPUT = "netcdf" # For 557 WW operations, for 2d ensemble gridspace

# Most variables are processed independently, and are listed below.
_VAR_ATTRIBUTES = {
    "Qs_tavg":
        "Qs          1  1  kg/m2s -  1  1 Qs          1  1  kg/m2s  - 1  1",
    "Qsb_tavg":
        "Qsb         1  1  kg/m2s -  1  1 Qsb         1  1  kg/m2s -  1  1",
    "RelSMC_inst":
        "RelSMC      1  4  -      -  0  4 RelSMC      1  4  -      -  0  4",
    "SmLiqFrac_inst":
        "SmLiqFrac   1  4  m3/m3  -  0  4 SmLiqFrac   1  4  m3/m3  -  0  4",
    "Snowcover_inst":
        "Snowcover   1  1  %      -  0  1 Snowcover   1  1  %      -  0  1",
    "SnowDepth_inst":
        "SnowDepth   1  1  m      -  0  1 SnowDepth   1  1  m      -  0  1",
    "SoilMoist_inst":
        "SoilMoist   1  4  m3/m3  -  0  4 SoilMoist   1  4  m3/m3  -  0  4",
    "SoilMoist_tavg":
        "SoilMoist   1  4  m3/m3  -  1  4 SoilMoist   1  4  m3/m3  -  1  4",
    "SoilTemp_inst":
        "SoilTemp    1  4  K      -  0  4 SoilTemp    1  4  K      -  0  4",
    "SoilTemp_tavg":
        "SoilTemp    1  4  K      -  1  4 SoilTemp    1  4  K      -  1  4",
    "SWE_inst":
        "SWE         1  1  kg/m2  -  0  1 SWE         1  1  kg/m2  -  0  1",
    "TotalPrecip_acc":
        "TotalPrecip 1  1  kg/m2  -  3  1 TotalPrecip 1  1  kg/m2  -  3  1",
}

# Smooth variables that are perturbed, derived from perturbed variables,
# or are LSM outputs that are affected by perturbed variables via physics.
_SMOOTH_VARS = ["Qs_tavg", "Qsb_tavg", "RelSMC_inst",
                "SmLiqFrac_inst", "SnowDepth_inst",
                "Snowcover_inst", "SoilMoist_inst",
                "SoilMoist_tavg", "SoilTemp_inst",
                "SoilTemp_tavg",  "SWE_inst"]

def _main():
    """Main Driver"""

    with open(_TEMPLATE, 'r', encoding="ascii") as file:
        lines = file.readlines()

    varlist = list(_VAR_ATTRIBUTES.keys())
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
            elif "LVT diagnostic file:" in line:
                line = f"LVT diagnostic file: logs/lvtlog.{var}.3hr.str"
            elif "LVT datastream attributes table::" in line:
                line = "LVT datastream attributes table::\n"
                line += f"{_VAR_ATTRIBUTES[var]}\n"
            elif "Metrics output directory:" in line:
                line = f"Metrics output directory: OUTPUT/STATS.{var}.3hr.str\n"
            elif "LIS output attributes file:" in line:
                line = "LIS output attributes file:"
                line += f" ./tables/MODEL_OUTPUT_LIST.TBL.lvt_557post.{var}.3hr\n"

            newlines.append(line)

        first_var = False
        if not os.path.exists("configs"):
            os.mkdir("configs")
        newfile = f"configs/lvt.config.{var}.3hr.str"
        print(f"Writing {newfile}")
        with open(newfile, "w", encoding="ascii") as file:
            for line in newlines:
                file.write(line)

if __name__ == "__main__":
    _main()
