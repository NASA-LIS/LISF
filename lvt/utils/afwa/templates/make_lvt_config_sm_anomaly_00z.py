#!/usr/bin/env python3

"""
#------------------------------------------------------------------------------
#
# SCRIPT: make_lvt_config_sm_anomaly.py
#
# PURPOSE: Customizes LVT config file to update 00Z soil moisture climatology,
# and output the climatology for the current month, plus anomalies.
#
# REVISION HISTORY:
# 28 Jun 2021: Eric Kemp (SSAI), first version.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import sys

NEWFILE = "configs/lvt.config.sm_anomaly"

VAR_ATTRIBUTES = \
    "SoilMoist   1  4  m3/m3  -  1  4 SoilMoist   1  4  m3/m3  -  1  4"

def _usage():
    """Print command line usage."""
    print("[INFO] Usage: %s restart_dir template yyyymmdd" %(sys.argv[0]))
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
    _restart_dir = sys.argv[1]
    if not os.path.exists(_restart_dir):
        print("[ERR] %s does not exist!" %(_restart_dir))
        sys.exit(1)

    # Check if lvt.config template exists.
    _template = sys.argv[2]
    if not os.path.exists(_template):
        print("[ERR] %s does not exist!" %(_template))
        sys.exit(1)

    # Parse the valid date
    yyyymmdd = sys.argv[3]
    _date = datetime.date(year=int(yyyymmdd[0:4]),
                              month=int(yyyymmdd[4:6]),
                              day=int(yyyymmdd[6:8]))

    return _restart_dir, _template, _date

if __name__ == "__main__":
    restart_dir, template, rundate = _read_cmd_args()

    rundt = datetime.datetime(year=rundate.year,
                              month=rundate.month,
                              day=rundate.day,
                              hour=0)
    startdt = rundt - datetime.timedelta(days=1)

    lines = open(template, 'r').readlines()
    newlines = []
    for LINE in lines:
        if "LVT running mode:" in LINE:
            LINE = 'LVT running mode: "Data intercomparison"\n'
        elif "LVT output format:" in LINE:
            LINE = "LVT output format: netcdf\n"
        elif "LVT output methodology:" in LINE:
            LINE = 'LVT output methodology: "2d gridspace"\n'
        elif "Analysis data sources:" in LINE:
            LINE = 'Analysis data sources: "LIS output" "none"\n'
        elif "Process HYCOM data:" in LINE:
            LINE = "Process HYCOM data 0\n"
        elif "Apply noise reduction filter: " in LINE:
            LINE = "Apply noise reduction filter: 0\n"
        elif "Start mode:" in LINE:
            LINE = "Start mode: restart\n"
        elif "LVT output restart files:" in LINE:
            LINE = "LVT output restart files: 1\n"
        elif "LVT restart output interval:" in LINE:
            LINE = "LVT restart output interval: 24hr\n"
        elif "LVT restart filename:" in LINE:
            LINE = "LVT restart filename: "
            LINE += "%s/LVT.%4.4d%2.2d%2.2d%2.2d00.rst\n" \
                %(restart_dir, startdt.year, startdt.month,
                  startdt.day, startdt.hour)
        elif "Starting year:" in LINE:
            LINE = "Starting year: %s\n" %(startdt.year)
        elif "Starting month:" in LINE:
            LINE = "Starting month: %s\n" %(startdt.month)
        elif "Starting day:" in LINE:
            LINE = "Starting day: %s\n" %(startdt.day)
        elif "Starting hour:" in LINE:
            LINE = "Starting hour: %s\n" %(startdt.hour)
        elif "Starting minute:" in LINE:
            LINE = "Starting minute: 00\n"
        elif "Starting second:" in LINE:
            LINE = "Starting second: 00\n"
        elif "Ending year:" in LINE:
            LINE = "Ending year: %s\n" %(rundt.year)
        elif "Ending month:" in LINE:
            LINE = "Ending month: %s\n" %(rundt.month)
        elif "Ending day:" in LINE:
            LINE = "Ending day: %s\n" %(rundt.day)
        elif "Ending hour:" in LINE:
            LINE = "Ending hour: %s\n" %(rundt.hour)
        elif "Ending minute:" in LINE:
            LINE = "Ending minute: 00\n"
        elif "Ending second:" in LINE:
            LINE = "Ending second: 00\n"
        elif "LVT clock timestep:" in LINE:
            LINE = "LVT clock timestep: 24hr\n"
        elif "LVT diagnostic file:" in LINE:
            LINE = "LVT diagnostic file: logs/lvtlog.sm_anomaly\n"
        elif "LVT datastream attributes table::" in LINE:
            LINE = "LVT datastream attributes table::\n"
            LINE += "%s\n" %(VAR_ATTRIBUTES)
        elif "Metrics attributes file:" in LINE:
            LINE = 'Metrics attributes file: tables/METRICS.TBL.anomaly\n'
        elif "Metrics computation frequency:" in LINE:
            LINE = "Metrics computation frequency: 24hr\n"
        elif "Metrics output directory:" in LINE:
            LINE = "Metrics output directory: OUTPUT/STATS.sm_anomaly\n"
        elif "Metrics output frequency:" in LINE:
            LINE = "Metrics output frequency: 24hr\n"
        elif "LIS output interval:" in LINE:
            LINE = "LIS output interval: 3hr\n"
        elif "LIS output attributes file:" in LINE:
            LINE = "LIS output attributes file: ./tables/"
            LINE += "MODEL_OUTPUT_LIST.TBL.lvt_557post.SoilMoist_inst.24hr\n"
        elif "LIS model timestep:" in LINE:
            LINE = "LIS model timestep: 15mn\n"

        newlines.append(LINE)


    print("[INFO] Writing %s" %(NEWFILE))
    f = open(NEWFILE, "w")
    for line in newlines:
        f.write(line)
    f.close()
