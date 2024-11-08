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
#------------------------------------------------------------------------------
#
# SCRIPT: make_lvt_config_sm_anomaly_00Z.py
#
# PURPOSE: Customizes LVT config file to update 00Z soil moisture climatology,
# and output the climatology for the current month, plus anomalies.
#
# REVISION HISTORY:
# 28 Jun 2021: Eric Kemp (SSAI), first version.
# 08 Dec 2022: Eric Kemp (SSAI), refactored to improve pylint score.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import sys

_NEWFILE = "configs/lvt.config.sm_anomaly"

_VAR_ATTRIBUTES = \
    "SoilMoist   1  4  m3/m3  -  1  4 SoilMoist   1  4  m3/m3  -  1  4"

def _usage():
    """Print command line usage."""
    print(f"[INFO] Usage: {sys.argv[0]} restart_dir template yyyymmdd")
    print("[INFO]   where:")
    print("[INFO]    restart_dir is the directory with the LVT restart file.")
    print("[INFO]    template is input lvt.config template file to customize.")
    print("[INFO]    yyyymmdd is the year/month/day to process.")

def _read_cmd_args():
    """Read command line arguments."""
    # Check if argument count is correct.
    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if restart_dir exists.
    restart_dir = sys.argv[1]
    if not os.path.exists(restart_dir):
        print(f"[ERR] {restart_dir} does not exist!")
        sys.exit(1)

    # Check if lvt.config template exists.
    template = sys.argv[2]
    if not os.path.exists(template):
        print(f"[ERR] {template} does not exist!")
        sys.exit(1)

    # Parse the valid date
    yyyymmdd = sys.argv[3]
    date = datetime.date(year=int(yyyymmdd[0:4]),
                         month=int(yyyymmdd[4:6]),
                         day=int(yyyymmdd[6:8]))

    return restart_dir, template, date

def _main():
    """Main driver"""
    restart_dir, template, rundate = _read_cmd_args()

    rundt = datetime.datetime(year=rundate.year,
                              month=rundate.month,
                              day=rundate.day,
                              hour=0)
    startdt = rundt - datetime.timedelta(days=1)

    with open(template, 'r', encoding="ascii") as file:
        lines = file.readlines()
    newlines = []
    for line in lines:
        if "LVT running mode:" in line:
            line = 'LVT running mode: "Data intercomparison"\n'
        elif "LVT output format:" in line:
            line = "LVT output format: netcdf\n"
        elif "LVT output methodology:" in line:
            line = 'LVT output methodology: "2d gridspace"\n'
        elif "Analysis data sources:" in line:
            line = 'Analysis data sources: "LIS output" "none"\n'
        elif "Process HYCOM data:" in line:
            line = "Process HYCOM data 0\n"
        elif "Apply noise reduction filter: " in line:
            line = "Apply noise reduction filter: 0\n"
        elif "Start mode:" in line:
            line = "Start mode: restart\n"
        elif "LVT output restart files:" in line:
            line = "LVT output restart files: 1\n"
        elif "LVT restart output interval:" in line:
            line = "LVT restart output interval: 24hr\n"
        elif "LVT restart filename:" in line:
            line = "LVT restart filename: "
            line += f"{restart_dir}/"
            line += f"LVT.{startdt.year:04}{startdt.month:02}"
            line += f"{startdt.day:02}{startdt.hour:02}00.rst\n"
        elif "Starting year:" in line:
            line = f"Starting year: {startdt.year}\n"
        elif "Starting month:" in line:
            line = f"Starting month: {startdt.month}\n"
        elif "Starting day:" in line:
            line = f"Starting day: {startdt.day}\n"
        elif "Starting hour:" in line:
            line = f"Starting hour: {startdt.hour}\n"
        elif "Starting minute:" in line:
            line = "Starting minute: 00\n"
        elif "Starting second:" in line:
            line = "Starting second: 00\n"
        elif "Ending year:" in line:
            line = f"Ending year: {rundt.year}\n"
        elif "Ending month:" in line:
            line = f"Ending month: {rundt.month}\n"
        elif "Ending day:" in line:
            line = f"Ending day: {rundt.day}\n"
        elif "Ending hour:" in line:
            line = f"Ending hour: {rundt.hour}\n"
        elif "Ending minute:" in line:
            line = "Ending minute: 00\n"
        elif "Ending second:" in line:
            line = "Ending second: 00\n"
        elif "LVT clock timestep:" in line:
            line = "LVT clock timestep: 24hr\n"
        elif "LVT diagnostic file:" in line:
            line = "LVT diagnostic file: logs/lvtlog.sm_anomaly\n"
        elif "LVT datastream attributes table::" in line:
            line = "LVT datastream attributes table::\n"
            line += f"{_VAR_ATTRIBUTES}\n"
        elif "Metrics attributes file:" in line:
            line = 'Metrics attributes file: tables/METRICS.TBL.anomaly\n'
        elif "Metrics computation frequency:" in line:
            line = "Metrics computation frequency: 24hr\n"
        elif "Metrics output directory:" in line:
            line = "Metrics output directory: OUTPUT/STATS.sm_anomaly\n"
        elif "Metrics output frequency:" in line:
            line = "Metrics output frequency: 24hr\n"
        elif "LIS output interval:" in line:
            line = "LIS output interval: 3hr\n"
        elif "LIS output attributes file:" in line:
            line = "LIS output attributes file: ./tables/"
            line += "MODEL_OUTPUT_LIST.TBL.lvt_557post.SoilMoist_inst.24hr\n"
        elif "LIS model timestep:" in line:
            line = "LIS model timestep: 15mn\n"

        newlines.append(line)

    if not os.path.exists("configs"):
        os.mkdir("configs")
    print(f"[INFO] Writing {_NEWFILE}")
    with open(_NEWFILE, "w", encoding="ascii") as file:
        for line in newlines:
            file.write(line)

if __name__ == "__main__":
    _main()
