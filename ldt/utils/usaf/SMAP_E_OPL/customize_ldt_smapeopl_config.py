#!/usr/bin/env python3

"""
SCRIPT: customize_ldt_smapeopl_config.py

Creates a customized ldt.config file for generating SMAP_E_OPL soil
moisture retrieval with LDT.

User must provide template ldt.config file, valid date/time, and
filelist suffix number via command line arguments.

To see usage statement, run script on command line w/o arguments.

REVISION HISTORY:
25 Jan 2024:  Eric Kemp.  Initial specification.
26 Jan 2024:  Eric Kemp.  Added filelist suffix number.
"""

import datetime
import os
import sys

def _usage():
    """Prints usage statement for script"""
    stmt = f"USAGE: {sys.argv[0]} YYYYMMDDHH ldt_config_tmpl_name"
    stmt += " filelist_suffix_number"
    print(stmt)

def _proc_valid_datetime():
    """Process command line argument for valid time"""
    if len(sys.argv[1]) != 10:
        print("ERR, invalid date/time for SMAP_E_OPL!")
        print(f"ERR, received {sys.argv[1]}")
        _usage()
        sys.exit(1)

    yyyymmddhh = sys.argv[1]
    try:
        yyyy = int(yyyymmddhh[0:4])
        mm = int(yyyymmddhh[4:6])
        dd = int(yyyymmddhh[6:8])
        hh = int(yyyymmddhh[8:])
        validdt = \
            datetime.datetime(year=yyyy,
                              month=mm,
                              day=dd,
                              hour=hh)
    except (TypeError, ValueError):
        print("ERR, invalid date/time for SMAP_E_OPL!")
        print(f"ERR, received {sys.argv[1]}")
        _usage()
        sys.exit(1)

    return validdt

def _proc_ldt_config_tmpl_name():
    """Process command line for ldt config template file name"""
    ldt_config_tmpl = sys.argv[2]
    if not os.path.exists(ldt_config_tmpl):
        print(f"ERR, cannot file template file {ldt_config_tmpl}")
        _usage()
        sys.exit(1)
    return ldt_config_tmpl

def _proc_filelist_suffix_number():
    """Process command line for filelist suffix number"""
    suffix_number = sys.argv[3]
    try:
        int_value = int(suffix_number)
    except ValueError:
        print("ERR, did not receive valid filelist suffix number!")
        _usage()
        sys.exit(1)
    if int_value < 0:
        print("ERR, filelist suffix must be nonnegative!")
        _usage()
        sys.exit(1)
    return int_value

def _proc_cmd_line():
    """Process command line and return valid datetime for SMAP_E_OPL"""
    if len(sys.argv) != 4:
        print("ERR, invalid number of command line arguments!")
        _usage()
        sys.exit(1)
    validdt = _proc_valid_datetime()
    ldt_config_tmpl = _proc_ldt_config_tmpl_name()
    filelist_suffix_number = _proc_filelist_suffix_number()
    return validdt, ldt_config_tmpl, filelist_suffix_number

def _create_new_ldt_config(validdt, ldt_config_tmpl, \
                           filelist_suffix_number):
    """Create new ldt.config customized to valid date"""
    with open(ldt_config_tmpl, "r", encoding='ascii') as file:
        lines = file.readlines()
    yyyymmddhh = f"{validdt.year:04d}{validdt.month:02d}" + \
        f"{validdt.day:02d}{validdt.hour:02d}"
    newfile = f"ldt.config.smapeopl.{yyyymmddhh}"
    with open(newfile, "w", encoding='ascii') as file:
        for line in lines:
            if "SMAP_E_OPL valid date (YYYYMMDDHH):" in line:
                newline = "SMAP_E_OPL valid date (YYYYMMDDHH): "
                newline += f"{yyyymmddhh}\n"
                file.write(newline)
                continue
            if "LDT diagnostic file:" in line:
                newline = "LDT diagnostic file: "
                newline += f"ldtlog.{yyyymmddhh}\n"
                file.write(newline)
                continue
            if "SMAP_E_OPL filelist suffix number:" in line:
                newline = "SMAP_E_OPL filelist suffix number: "
                newline += f"{filelist_suffix_number}\n"
                file.write(newline)
                continue
            # Pass through all other lines
            file.write(line)
            continue

def _main():
    """Main driver"""
    validdt, ldt_config_tmpl, filelist_suffix_number = _proc_cmd_line()
    _create_new_ldt_config(validdt, ldt_config_tmpl, \
                           filelist_suffix_number)

if __name__ == "__main__":
    _main()
