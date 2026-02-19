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
#--------------------------------------------------------------------------
#
# SCRIPT: update_lis_config.py
#
# PURPOSE:  Updates the LIS MR run configuration (lis.config) to use
# either GALWEM GE or GALWEM GD as forcing.
#
# REQUIREMENTS as of 08 July 2025:
# * Python 3.11 or higher
#
# REVISION HISTORY:
# 08 July 2025:  Yeosang Yoon, first version.
# 10 July 2025:  Eric Kemp, clean up to pacify pylint.
#
#--------------------------------------------------------------------------
"""

import argparse
from datetime import datetime, timedelta
import sys

def update_config(date, cycle, met_type, lsm):
    """Customize the lis.config file"""

    # Validate met forcing
    if met_type not in ['GALWEM-GE', 'GALWEM']:
        txt = "Invalid met forcing. Only 'GALWEM-GE' or 'GALWEM' " + \
            "are supported."
        raise ValueError(txt)

    # Validate LSM
    valid_lsms = ['noah39', 'noahmp401']
    if lsm not in valid_lsms:
        txt = "Invalid LSM. Only 'noah39' or 'noahmp401' " + \
            " are supported."
        raise ValueError(txt)

    # Parse dates
    start_dt = datetime.strptime(date + cycle, "%Y%m%d%H")
    forecast_days = 16 if met_type == "GALWEM-GE" else 10
    end_dt = start_dt + timedelta(days=forecast_days)
    date_str = start_dt.strftime("%Y%m%d%H")

    # File paths
    input_file = f"./input/template/lis.config.mr.{lsm}.rapid.template"
    output_file = \
       f"lis.config.mr.{lsm}.rapid.{met_type.lower().replace('-', '_')}.76"
    restart_prefix = "LIS_RST_NOAH39" if lsm == "noah39" \
        else "LIS_RST_NOAHMP401"

    with open(input_file, "r", encoding="ascii") as f:
        lines = f.readlines()

    new_lines = []
    forcings_inserted = False

    for _, line in enumerate(lines):
        # Start/End time edits (shared)
        if 'Starting year:' in line:
            line = 'Starting year:                            ' + \
                f'{start_dt.year}\n'
        elif 'Starting month:' in line:
            line = 'Starting month:                           ' + \
                f'{start_dt.month:02d}\n'
        elif 'Starting day:' in line:
            line = 'Starting day:                             ' + \
                f'{start_dt.day:02d}\n'
        elif 'Starting hour:' in line:
            line = 'Starting hour:                            ' + \
                f'{start_dt.hour:02d}\n'
        elif 'Ending year:' in line:
            line = 'Ending year:                              ' + \
                f'{end_dt.year}\n'
        elif 'Ending month:' in line:
            line = 'Ending month:                             ' + \
                f'{end_dt.month:02d}\n'
        elif 'Ending day:' in line:
            line = 'Ending day:                               ' + \
                f'{end_dt.day:02d}\n'
        elif 'Ending hour:' in line:
            line = 'Ending hour:                              ' + \
                f'{end_dt.hour:02d}\n'
        # Common output settings
        elif 'Output directory:' in line:
            line = 'Output directory:                       ' + \
                f'./output/{date_str}/' + \
                f'{met_type.lower().replace("-", "_")}/rapid/{lsm}\n'
        elif 'Diagnostic output file:' in line:
            line = 'Diagnostic output file:                 ' + \
                f'./output/{date_str}/' + \
                f'{met_type.lower().replace("-", "_")}' + \
                f'/rapid/{lsm}/log/lislog\n'
        # Noah restart file
        elif ('Noah.3.9 restart file:' in line and \
              lsm == 'noah39') or \
             ('Noah-MP.4.0.1 restart file:' in line and \
              lsm == 'noahmp401'):
            if met_type == "GALWEM":
                rst_file = f"./input/rstfile/{date_str}/{lsm}/" + \
                    f"{restart_prefix}_{date_str}00_EN01.d01.nc"
            else:  # GALWEM-GE
                rst_file = f"./input/rstfile/{date_str}/{lsm}/" + \
                    f"{restart_prefix}_{date_str}00_EN10.d01.nc"
            # Set correct prefix
            label = 'Noah.3.9 restart file:' if lsm == 'noah39' \
                else 'Noah-MP.4.0.1 restart file:'
            line = f'{label:<45}{rst_file}\n'

        # RAPID restart file
        elif 'RAPID routing model restart file:' in line:
            if met_type == "GALWEM":
                rst = f"./input/rstfile/{date_str}/{lsm}/" + \
                    f"LIS_RST_RAPID_router_{date_str}00.d01.nc"
            else:
                rst = f"./input/rstfile/{date_str}/{lsm}/" + \
                    f"LIS_RST_RAPID_router_{date_str}00_EN10.d01.nc"
            line = f'RAPID routing model restart file:       {rst}\n'

        # RAPID ensemble mode
        elif 'RAPID run in ensemble mode:' in line:
            if met_type == "GALWEM":
                line = 'RAPID run in ensemble mode:             ' + \
                    '1 # 0=open loop; 1=ensemble mean; 2=ensemble\n'
            else:
                line = 'RAPID run in ensemble mode:             ' + \
                    '2 # 0=open loop; 1=ensemble mean; 2=ensemble\n'
        # --- Met forcing-specific sections ---
        elif met_type == "GALWEM-GE":
            if 'Met forcing sources:' in line:
                line = 'Met forcing sources:                    ' + \
                    '"GALWEM-GE forecast"\n'
            elif 'Output methodology:' in line:
                line = 'Output methodology:                     ' + \
                    '"2d ensemble gridspace"\n'
            elif 'Number of ensembles per tile:' in line:
                line = 'Number of ensembles per tile:           10\n'
            elif '#--------------------------------FORCINGS---------' \
                 in line and not forcings_inserted:
                new_lines.append(line)
                new_lines.extend([
                    "\n# GALWEM-GE forecast\n",
                    "GALWEM-GE forecast forcing directory:          " + \
                    "./input/GALWEM_GE\n",
                    "GALWEM-GE forecast run mode:                   " + \
                    "forecast\n",
                    "GALWEM-GE forecast number of ensemble members: " + \
                    "10\n\n",
                    "Apply GALWEM-GE precipitation bias correction: " + \
                    "1 #enter 1 - use; or 0\n",
                    "GALWEM-GE model CDF directory:                 " + \
                    "./input/cdf\n\n"
                ])
                forcings_inserted = True
                continue

        elif met_type == "GALWEM":
            if 'Met forcing sources:' in line:
                line = 'Met forcing sources:                    ' + \
                    '"GALWEM forecast"\n'
            elif 'Output methodology:' in line:
                line = 'Output methodology:                     ' + \
                    '"2d gridspace"\n'
            elif 'Number of ensembles per tile:' in line:
                line = 'Number of ensembles per tile:           1\n'
            elif '#--------------------------------FORCINGS' \
                 in line and not forcings_inserted:
                new_lines.append(line)
                new_lines.extend([
                    "\n# GALWEM forecast\n",
                    "GALWEM forecast forcing directory:         " + \
                    "./input/GALWEM_GD/\n",
                    "GALWEM forecast resolution:                " + \
                    "17            # 17(=17km) or 25(=25 deg)\n",
                    "GALWEM forecast run mode:                  " + \
                    "forecast\n\n",
                ])
                forcings_inserted = True
                continue

        new_lines.append(line)

    with open(output_file, "w", encoding="ascii") as f:
        f.writelines(new_lines)

    print(f"[SUCCESS] Generated config saved to '{output_file}'")
    return output_file

if __name__ == "__main__":
    DESCRIPTION = "Generate LIS config from template."
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--date", required=True, \
                        help="Forecast date (YYYYMMDD)")
    parser.add_argument("--cycle", required=True, \
                        help="Cycle hour (e.g., 00, 12)")
    parser.add_argument("--met", required=True, \
                        help="Met forcing: 'GALWEM' or 'GALWEM-GE'")
    parser.add_argument("--lsm", required=True, \
                        choices=["noah39", "noahmp401"], \
                        help="Land surface model")

    args = parser.parse_args()

    try:
        update_config(args.date, args.cycle, args.met, args.lsm)
    except (ValueError, UnicodeError, OSError) as e:
        print(f"[ERR] {e}")
        sys.exit(1)
