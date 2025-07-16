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
# SCRIPT: fault_tolerance_galwem_ge.py
#
# PURPOSE:  Checks for GALWEM-GE GRIB2 files, and based on results,
# updates the LIS MR run configuration (lis.config) to use either
# GALWEM-GE or GALWEM-GD as forcing.
#
# REQUIREMENTS as of 10 July 2025:
# * Python 3.11 or higher
# * xarray Python library
# * cfgrib Python library
#
# REVISION HISTORY:
# 24 June 2025:  Yeosang Yoon, first version.
# 26 June 2025:  Eric Kemp, updates to satisfy pylint.
# 08 July 2025:  Yeosang Yoon, update the codes to generate lis.config
#                              insted of submitting to Slurm.
# 10 July 2025:  Eric Kemp, minor cleanup.
#
#--------------------------------------------------------------------------
"""

import argparse
from datetime import datetime
import os
import sys

import xarray as xr
from update_lis_config import update_config

def generate_forecast_hours():
    """Generate a list of forecast hours."""
    return [f"{i:03d}" for i in range(0, 195, 3)] + \
        [f"{i:03d}" for i in range(198, 385, 6)]

def generate_members():
    """Generate a list of member directory names."""
    return ["member000"] + \
        [f"member{str(i).zfill(3)}" for i in range(28, 37)]

def build_filename(base_prefix, grib_suffix, memb_id, fh):
    """Construct the expected GRIB2 filename for a given member and
forecast hour."""
    return f"{base_prefix}-MEMB{memb_id}{grib_suffix}.{fh}_DF.GR2"

def is_grib2_valid(filepath, deep_check=False):
    """Test if opened file is a GRIB file"""
    try:
        # Step 1: quick size filter
        if os.path.getsize(filepath) < 10 * 1024:
            return False
        # Step 2: magic header check
        with open(filepath, 'rb') as f:
            if f.read(4) != b'GRIB':
                return False
        # Step 3: optional deep read
        if deep_check:
            # We just need to check one GRIB2 message
            filters = {
                'typeOfLevel' : 'surface',
                'stepType' : 'accum',
            }
            ds = xr.open_dataset(filepath, engine='cfgrib', \
                                 filter_by_keys=filters, \
                                 decode_timedelta=False)
            del ds
        return True
    except (OSError, ValueError, FileNotFoundError) as e:
        print(f"An unexpected error occurred:  {e}")
        return False

def check_member_files(member_path, memb_id, forecast_hours, \
                       base_prefix, \
                       grib_suffix, deep_check=False):
    """Checks files for a given GALWEM GE member"""
    missing_or_corrupted = []

    for fh in forecast_hours:
        filename = build_filename(base_prefix, grib_suffix, memb_id, fh)
        filepath = os.path.join(member_path, filename)

        if not os.path.isfile(filepath):
            missing_or_corrupted.append((filename, "missing"))
        elif not is_grib2_valid(filepath, deep_check=deep_check):
            missing_or_corrupted.append((filename, "corrupted"))

    return missing_or_corrupted

def check_files(base_path, date_str, cycle_str, deep_check=False):
    """Main function to check all members and their forecast files."""
    forecast_hours = generate_forecast_hours()
    members = generate_members()

    base_prefix = "PS.557WW_SC.U_DI.C_GP.GALWEM-GE"
    grib_suffix = f"_GR.C20KM_AR.GLOBAL_DD.{date_str}_CY.{cycle_str}_FH"

    missing_members = []
    bad_files_by_member = {}

    for member in members:
        member_path = os.path.join(base_path, member)
        memb_id = member[-3:]

        if not os.path.isdir(member_path):
            missing_members.append(member)
            continue

        bad_files = check_member_files(member_path, memb_id, \
                                       forecast_hours, base_prefix, \
                                       grib_suffix, deep_check=deep_check)
        if bad_files:
            bad_files_by_member[member] = bad_files

    return missing_members, bad_files_by_member

def print_results(missing_members, bad_files):
    """Prints summary of file checks to screen"""
    print("==== File Check Result ====")

    if missing_members:
        print(f"[WARN] Missing member folders: {', '.join(missing_members)}")
    else:
        print("[INFO] All member folders are present.")

    if bad_files:
        print("\n[WARN] Missing or corrupted files by member:")
        for member, files in bad_files.items():
            print(f"  - {member}: {len(files)} problematic files")
            for f, reason in files:
                print(f"      {f} ({reason})")
    else:
        print("[INFO] All files are present and valid.")

def write_log_file(log_path, date_str, cycle_str, missing_members, \
                   bad_files, config_path):
    """Writes information on available ensemble files to a log file"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_filename = f"check_result_{date_str}_CY{cycle_str}.log"
    full_path = os.path.join(log_path, log_filename)

    with open(full_path, "w", encoding="ascii") as log:
        txt = f"[{timestamp}] File Check Result for {date_str}" + \
              f" Cycle {cycle_str}\n"
        log.write(txt)
        log.write("="*60 + "\n")

        if missing_members:
            txt = f"Missing member folders ({len(missing_members)}):" + \
                f"  {', '.join(missing_members)}\n"
            log.write(txt)
        else:
            log.write("All member folders are present.\n")

        if bad_files:
            log.write("\nMissing or corrupted files by member:\n")
            for member, files in bad_files.items():
                txt = f"  - {member}: {len(files)} problematic files\n"
                log.write(txt)
                for f, reason in files:
                    log.write(f"      {f} ({reason})\n")
        else:
            log.write("All GRIB2 files are valid.\n")

        log.write(f"\nGenerated config saved to '{config_path}'\n")
    print(f"[INFO] Log saved to: {full_path}")

def main():
    """Main driver"""
    txt = "Check GALWEM-GE GRIB2 files and generate LIS config"
    parser = argparse.ArgumentParser(description=txt)
    txt = "Root directory containing <date>/<cycle>/memberXXX folders"
    parser.add_argument("--root", type=str, required=True, \
                        help=txt)
    parser.add_argument("--date", type=str, required=True, \
                        help="Forecast start date (YYYYMMDD)")
    parser.add_argument("--cycle", type=str, required=True, \
                        help="Cycle (e.g., 00, 12)")
    txt = "Directory to save log files (default: current)"
    parser.add_argument("--log-dir", type=str, default=".", \
                        help=txt)
    parser.add_argument("--lsm", type=str, required=True, \
                        choices=["noah39", "noahmp401"], \
                        help="Land surface model")
    parser.add_argument("--deep-check", action="store_true", \
                        help="Enable deep GRIB2 read check using cfgrib")

    args = parser.parse_args()

    missing_members, bad_files = check_files(args.root, args.date, \
                                             args.cycle, \
                                             deep_check=args.deep_check)
    print_results(missing_members, bad_files)

    if not missing_members and not bad_files:
        met_type = "GALWEM-GE"
    else:
        met_type = "GALWEM"

    try:
        config_path = update_config(args.date, args.cycle, met_type, \
                                    args.lsm)
    except (ValueError, UnicodeError, OSError) as e:
        print(f"[ERR] Failed to generate LIS config: {e}")
        sys.exit(1)

    write_log_file(args.log_dir, args.date, args.cycle, missing_members, \
                   bad_files, config_path)

if __name__ == "__main__":
    main()
