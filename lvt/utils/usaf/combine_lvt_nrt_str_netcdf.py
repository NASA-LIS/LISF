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
# SCRIPT: combine_lvt_nrt_netcdf.py
#
# PURPOSE:  Runs the ncks utility (part of the netCDF Operator utilities)
# to consolidate LVT output netCDF files into a master file for netCDF
# streamflow users.  The bulk of the consolidation (that is, variable,
# dimension, and attribute copying) is done by ncks,
# while this script keeps track of the input and output files and the
# required variables.
#
# REQUIREMENTS:
# * Python 3.13
# * NetCDF Operator (NCO) utilities
#
# REVISION HISTORY:
# 18 Aug 2026:  Eric Kemp (SSAI), based on earlier run_ncks.py script.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import subprocess
import sys

#------------------------------------------------------------------------------

# Path to NCO ncks program
_NCKS_PATH = "/usr/local/other/nco/5.1.7/bin/ncks" # On Discover

# Supported LIS LSMs
_LIS_LSMS = ["NOAH", "NOAHMP"]

# Supported LIS ROUTING models
_LIS_ROUTING = ["RAPID"]

# The LVT invocations for Noah LSM output.  Each invocation handles a subset
# of the total variable list due to memory limitations.
_LVT_NOAH_INVOCATIONS_3HR = ['Qs_tavg', 'Qsb_tavg',
                             'RelSMC_inst', 'SmLiqFrac_inst',
                             'SnowDepth_inst', 'Snowcover_inst',
                             'SoilMoist_inst', 'SoilMoist_tavg',
                             'SoilTemp_inst', 'SoilTemp_tavg',
                             'SWE_inst', 'TotalPrecip_acc']


# The LVT invocations for NoahMP LSM output.  Each invocation handles a subset
# of the total variable list due to memory limitations.
_LVT_NOAHMP_INVOCATIONS_3HR = ['Qs_tavg', 'Qsb_tavg',
                               'RelSMC_inst', 'SmLiqFrac_inst',
                               'SnowDepth_inst', 'Snowcover_inst',
                               'SoilMoist_inst', 'SoilMoist_tavg',
                               'SoilTemp_inst', 'SoilTemp_tavg',
                               'SWE_inst', 'TotalPrecip_acc']

# The combined invocation dictionary for all supported LSMs.
_INVOCATIONS = {
    "NOAH_3HR": _LVT_NOAH_INVOCATIONS_3HR,
    "NOAHMP_3HR": _LVT_NOAHMP_INVOCATIONS_3HR,
}

# The Noah variables handled by each LVT invocation.
_LIS_NOAH_VARIABLES_3HR = {}
for var in _LVT_NOAH_INVOCATIONS_3HR:
    _LIS_NOAH_VARIABLES_3HR[var] = [var]

# The NoahMP variables handled by each LVT invocation.
_LIS_NOAHMP_VARIABLES_3HR = {}
for var in _LVT_NOAHMP_INVOCATIONS_3HR:
    _LIS_NOAHMP_VARIABLES_3HR[var] = [var]

# Combined breakdown of all variables handled by LSM and LVT invocation.
_LIS_VARIABLES = {
    "NOAH_3HR": _LIS_NOAH_VARIABLES_3HR,
    "NOAHMP_3HR": _LIS_NOAHMP_VARIABLES_3HR,
}

# These variables are processed by all invocations, and are generated
# automatically.
_OTHER_VARIABLES = ["latitude", "longitude", "time"]

#------------------------------------------------------------------------------
def _usage():
    """Print command line usage"""
    print(f"Usage: {sys.argv[0]} yyyymmddhh fhr lsm routing")
    print("   where:")
    print("           yyyymmddhh is start year/month/day/hour in UTC")
    print("           fhr forecast hour")
    print("           lsm is name of land surface model used by LIS")
    print("           routing is name of routing model used by LIS")

#------------------------------------------------------------------------------
def _read_cmd_args():
    """Read command line arguments"""
    # Check if argument count is correct
    if len(sys.argv) not in [5]:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Convert yyyymmddhh argument to a datetime object
    yyyymmddhh = sys.argv[1]
    try:
        year = int(yyyymmddhh[0:4])
        month = int(yyyymmddhh[4:6])
        day = int(yyyymmddhh[6:8])
        hour = int(yyyymmddhh[8:10])
        startdt = datetime.datetime(year, month, day, hour)
    except ValueError:
        print("[ERR] Cannot process yyyymmddhh argument!")
        _usage()
        sys.exit(1)

    try:
        fhr = int(sys.argv[2])
    except ValueError:
        print("[ERR] Cannot process fhr argument!")
        _usage()
        sys.exit(1)

    # Get lsm name
    lsm = None
    if sys.argv[3] in _LIS_LSMS:
        lsm = sys.argv[3]
    if lsm is None:
        print("[ERR] Invalid lsm selection!")
        print(f" lsm value is {sys.argv[3]}")
        text = " Supported lsms:"
        for lsm in _LIS_LSMS:
            text += f" {lsm}"
        print(text)
        sys.exit(1)

    # Get routing name
    routing = None
    if sys.argv[4] in _LIS_ROUTING:
        routing = sys.argv[4]
    if routing is None:
        print("[ERR] Invalid routing selection!")
        print(f" routing value is {sys.argv[4]}")
        text = " Supported routing models:"
        for routing in _LIS_ROUTING:
            text += f" {routing}"
        print(text)
        sys.exit(1)

    # See if ncks exists and is executable by current user (the script)
    # This used to be specified on the command line, but is now hardwired
    # to better comply with Air Force security requirements.
    ncks = _NCKS_PATH
    _check_ncks_path(ncks)

    return startdt, fhr, lsm, routing

#------------------------------------------------------------------------------
def _check_ncks_path(ncks):
    """Check if ncks works"""
    if not os.path.isfile(ncks):
        print(f"[ERR] Binary {ncks} does not exist!")
        print(f"[ERR] Modify {sys.argv[0]} to correct the path to ncks!")
        sys.exit(1)
    if not os.access(ncks, os.X_OK):
        print(f"[ERR] {ncks} cannot be executed by current user!")
        print(f"[ERR] Modify {sys.argv[0]} to correct the path to ncks!")
        sys.exit(1)

#------------------------------------------------------------------------------
def _get_nc_mean_files(startdt, fhr, lsm, routing):
    """Collect netCDF mean files"""
    key = f"{lsm}_3HR"
    invocation_list = _INVOCATIONS[key]

    mean_nc_infiles = {}

    # Collect input files
    for invocation in invocation_list:

        path = f"OUTPUT/STATS.{invocation}.3hr.str"
        path += "/PS.557WW"
        path += "_SC.U"
        path += "_DI.C"
        path += f"_GP.LIS-NRT-{lsm}-{routing}"
        path += "_GR.C0P09DEG"
        path += "_AR.GLOBAL"
        path += "_PA.SURFACEMODEL"
        path += f"_DD.{startdt.year:04}{startdt.month:02}{startdt.day:02}"
        path += f"_CY.{startdt.hour:02}"
        path += f"_FH.{fhr:03}"
        path += f"_DF.NC"
        mean_path = path
        if not os.path.exists(mean_path):
            print(f"[ERR], {mean_path} does not exist!")
            sys.exit(1)

        mean_nc_infiles[invocation] = mean_path

    # Get output file
    path = f"OUTPUT/STATS_merged_3hr.str"
    if not os.path.exists(path):
        os.mkdir(path)
    path += "/PS.557WW"
    path += "_SC.U"
    path += "_DI.C"
    path += f"_GP.LIS-NRT-{lsm}-{routing}"
    path += "_GR.C0P09DEG"
    path += "_AR.GLOBAL"
    path += "_PA.SURFACEMODEL"
    path += f"_DD.{startdt.year:04}{startdt.month:02}{startdt.day:02}"
    path += f"_CY.{startdt.hour:02}"
    path += f"_FH.{fhr:03}"
    path += f"_DF.NC"

    mean_nc_outfile = path

    # All done
    return mean_nc_infiles, mean_nc_outfile

#------------------------------------------------------------------------------
def _merge_nc_files(lsm, routing, nc_infiles,
                    nc_outfile):
    """Use ncks to merge netCDF fields together"""

    ncks = _NCKS_PATH

    key = f"{lsm}_3HR"

    # Start with ensemble mean
    cmd = f"cp {nc_infiles[_INVOCATIONS[key][0]]} {nc_outfile}"
    print(cmd)
    err = subprocess.call(cmd, shell=True)
    if err != 0:
        print("[ERR] Problem with cp!")
        sys.exit(1)

    invocations = _INVOCATIONS[key][1:]
    for invocation in invocations:

        variables = _LIS_VARIABLES[key][invocation]

        for variable in variables:
            cmd = f"{ncks} -A -v {variable} {nc_infiles[invocation]} "
            cmd += f"{nc_outfile}"
            print(cmd)
            err = subprocess.call(cmd, shell=True)
            if err != 0:
                print("[ERR] Problem with ncks!")
                sys.exit(1)

#------------------------------------------------------------------------------
# Main Driver.
def _main():
    """Main driver"""
    # Process command line arguments
    startdt, fhr, lsm, routing = _read_cmd_args()

    # Collect netCDF files
    (mean_nc_infiles, mean_nc_outfile) = \
        _get_nc_mean_files(startdt, fhr, lsm, routing)

    # Merge the input netCDF files together
    _merge_nc_files(lsm, routing, mean_nc_infiles, mean_nc_outfile)

#------------------------------------------------------------------------------
if __name__ == "__main__":
    _main()
