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
Sample script to customize lvt.config files for jules50 postprocessing for
557WW.
"""

import datetime
import os

_TEMPLATE = "templates/lvt.config.template.jules50"

_STARTDT = datetime.datetime(2022, 8, 2,  6)
_ENDDT = datetime.datetime(2022, 8, 2,  12)

#_OUTPUT = "netcdf" # For 557 Ops (Fields file)
_OUTPUT = "grib2"

# Most variables are processed independently, and are listed below.
_VAR_ATTRIBUTES = {
    "AvgSurfT_inst":
    "AvgSurfT    1  1  K      -  0  1 AvgSurfT    1  1  K      -  0  1",
    "AvgSurfT_tavg":
        "AvgSurfT    1  1  K      -  1  1 AvgSurfT    1  1  K      -  1  1",
    "Albedo_tavg":
        "Albedo      1  1  %      -  1  1 Albedo      1  1  %      -  1  1",
    "CanopInt_inst":
        "CanopInt    1  1  kg/m2  -  0  1 CanopInt    1  1  kg/m2  -  0  1",
    "Elevation_inst":
        "Elevation   1  1  m      -  0  1 Elevation   1  1  m      -  0  1",
    "Evap_tavg":
        "Evap        1  1  kg/m2s -  1  1 Evap        1  1  kg/m2s -  1  1",
    "Landmask_inst":
        "Landmask    1  1  -      -  0  1 Landmask    1  1  -      -  0  1",
    "Landcover_inst":
        "Landcover   1  1  -      -  0  1 Landcover   1  1  -      -  0  1",
    "LWdown_f_inst":
        "LWdown_f    1  1  W/m2   -  0  1 LWdown_f    1  1  W/m2   -  0  1 ",
    "LWdown_f_tavg":
        "LWdown_f    1  1  W/m2   -  1  1 LWdown_f    1  1  W/m2   -  1  1 ",
    "Psurf_f_inst":
        "Psurf_f     1  1  Pa     -  0  1 Psurf_f     1  1  Pa     -  0  1",
    "Psurf_f_tavg":
        "Psurf_f     1  1  Pa     -  1  1 Psurf_f     1  1  Pa     -  1  1",
    "Qair_f_inst":
        "Qair_f      1  1  kg/kg  -  0  1 Qair_f      1  1  kg/kg  -  0  1",
    "Qair_f_tavg":
        "Qair_f      1  1  kg/kg  -  1  1 Qair_f      1  1  kg/kg  -  1  1",
    "Qh_tavg":
        "Qh          1  1  W/m2   -  1  1 Qh          1  1  W/m2   -  1  1",
    "Qle_tavg":
        "Qle         1  1  W/m2   -  1  1 Qle         1  1  W/m2   -  1  1",
    "Qs_acc":
        "Qs          1  1  kg/m2  -  3  1 Qs          1  1  kg/m2  -  3  1",
    "Qsb_acc":
        "Qsb         1  1  kg/m2  -  3  1 Qsb         1  1  kg/m2  -  3  1",
    "RelSMC_inst":
        "RelSMC      1  4  -      -  0  4 RelSMC      1  4  -      -  0  4",
    "SmLiqFrac_inst":
        "SmLiqFrac   1  4  m3/m3  -  0  4 SmLiqFrac   1  4  m3/m3  -  0  4",
    "SnowDepth_inst":
        "SnowDepth   1  1  m      -  0  1 SnowDepth   1  1  m      -  0  1",
    "Snowcover_inst":
        "Snowcover   1  1  %      -  0  1 Snowcover   1  1  %      -  0  1",
    "SoilMoist_inst":
        "SoilMoist   1  4  m3/m3  -  0  4 SoilMoist   1  4  m3/m3  -  0  4",
    "SoilMoist_tavg":
        "SoilMoist   1  4  m3/m3  -  1  4 SoilMoist   1  4  m3/m3  -  1  4",
    "SoilTemp_inst":
        "SoilTemp    1  4  K      -  0  4 SoilTemp    1  4  K      -  0  4",
    "SoilTemp_tavg":
        "SoilTemp    1  4  K      -  1  4 SoilTemp    1  4  K      -  1  4",
    "SWdown_f_inst":
        "SWdown_f    1  1  W/m2   -  0  1 SWdown_f    1  1  W/m2   -  0  1",
    "SWdown_f_tavg":
        "SWdown_f    1  1  W/m2   -  1  1 SWdown_f    1  1  W/m2   -  1  1",
    "SWE_inst":
        "SWE         1  1  kg/m2  -  0  1 SWE         1  1  kg/m2  -  0  1",
    "Tair_f_inst":
        "Tair_f      1  1  K      -  0  1 Tair_f      1  1  K      -  0  1",
    "Tair_f_max":
        "Tair_f_max  1  1  K      -  1  1 Tair_f_max  1  1  K      -  1  1",
    "Tair_f_tavg":
        "Tair_f      1  1  K      -  1  1 Tair_f      1  1  K      -  1  1",
    "TotalPrecip_acc":
        "TotalPrecip 1  1  kg/m2  -  3  1 TotalPrecip 1  1  kg/m2  -  3  1",
    "Wind_f_inst":
        "Wind_f      1  1  m/s    -  0  1 Wind_f      1  1  m/s    -  0  1",
    "Wind_f_tavg":
        "Wind_f      1  1  m/s    -  1  1 Wind_f      1  1  m/s    -  1  1",

    # "ActSnowNL_inst":
    #     "ActSnowNL   1  1    -    -  0  1 ActSnowNL   1  1   -     -  0  1",
    # "GrndSnow_inst":
    #     "GrndSnow    1  1  kg/m2  -  0  1 GrndSnow    1  1  kg/m2  -  0  1",
    # "LayerSnowDensity_inst":
    #     "LayerSnowDensity 1 1 kg/m3 - 0 3 LayerSnowDensity 1 1 kg/m3 - 0 3",
    # "LayerSnowDepth_inst":
    #     "LayerSnowDepth   1 1 m   -  0  3 LayerSnowDepth  1 1 m    -  0  3",
    # "LayerSnowGrain_inst":
    #     "LayerSnowGrain 1 1 microns - 0 3 LayerSnowGrain 1 1 microns - 0 3",
    # "SnowDensity_inst":
    #     "SnowDensity 1  1  kg/m3  -  0  1 SnowDensity 1  1 kg/m3   -  0  1",
    # "SnowGrain_inst":
    #     "SnowGrain   1  1 microns -  0  1 SnowGrain   1  1 microns -  0  1",
    # "SnowIce_inst":
    #     "SnowIce     1  1  kg/m2  -  0  3 SnowIce     1  1 kg/m2   -  0  3",
    # "SnowLiq_inst":
    #     "SnowLiq     1  1  kg/m2  -  0  3 SnowLiq     1  1 kg/m2   -  0  3",
    # "SnowTProf_inst":
    #     "SnowTProf   1  1    K    -  0  3 SnowTProf   1  1   K     -  0  3",
    # "SurftSnow_inst":
    #     "SurftSnow   1  1  kg/m2  -  0  1 SurftSnow   1  1  kg/m2  -  0  1",
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
_SMOOTH_VARS = ["AvgSurfT_inst", "AvgSurfT_tavg",
               "Albedo_tavg", "CanopInt_inst",
               "Evap_tavg", "LWdown_f_inst",
               "LWdown_f_tavg", "Qh_tavg", "Qle_tavg",
               "Qs_acc", "Qsb_acc", "RelSMC_inst",
               "SmLiqFrac_inst", "SnowDepth_inst",
               "Snowcover_inst", "SoilMoist_inst",
               "SoilMoist_tavg", "SoilTemp_inst",
               "SoilTemp_tavg", "SWdown_f_inst",
               "SWdown_f_tavg", "SWE_inst",
               "Tair_f_inst", "Tair_f_max",
               "Tair_f_tavg", "TotalPrecip_acc",
               "Tair_f_min", "RHMin_inst"]

#               "GrndSnow_inst", "LayerSnowDensity_inst",
#               "LayerSnowDepth_inst", "LayerSnowGrain_inst",
#               "SnowDensity_inst", "SnowGrain_inst",
#               "SnowIce_inst", "SnowLiq_inst",
#               "SnowTProf_inst", "SurftSnow_inst"]

def _main():
    """Main driver"""

    with open(_TEMPLATE, 'r', encoding="ascii") as file:
        lines = file.readlines()

    varlist = list(_VAR_ATTRIBUTES.keys())
    varlist.append("RHMin_inst")  # RHMin will be handled specially below
    varlist.sort()
    first_var = True
    for var in varlist:
        newlines = []
        for line in lines:
            if "LVT output format:" in line:
                line = f"LVT output format: {_OUTPUT}\n"
            elif "Process HYCOM data:" in line:
                if first_var:
                    line = "Process HYCOM data: 1\n"
                else:
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
                line = f"Ending year: {_STARTDT.year}\n"
            elif "Ending month:" in line:
                line = f"Ending month: {_ENDDT.month}\n"
            elif "Ending day:" in line:
                line = f"Ending day: {_ENDDT.day}\n"
            elif "Ending hour:" in line:
                line = f"Ending hour: {_ENDDT.hour}\n"
            elif "LVT diagnostic file:" in line:
                line = f"LVT diagnostic file: logs/lvtlog.{var}.3hr"
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
            elif "Metrics output directory:" in line:
                line = f"Metrics output directory: OUTPUT/STATS.{var}.3hr\n"
            elif "LIS output attributes file:" in line:
                line = "LIS output attributes file:"
                line += f" ./templates/MODEL_OUTPUT_LIST.TBL.lvt_557post.{var}.3hr\n"

            newlines.append(line)

        first_var = False
        if not os.path.exists("configs"):
            os.mkdir("configs")
        newfile = f"configs/lvt.config.{var}.3hr"
        print(f"Writing {newfile}")
        with open(newfile, "w", encoding="ascii") as file:
            for line in newlines:
                file.write(line)

if __name__ == "__main__":
    _main()
