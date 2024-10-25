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
# SCRIPT: run_ncks.py
#
# PURPOSE:  Runs the ncks utility (part of the netCDF Operator utilities)
# to consolidate LVT output netCDF files into two master files (one for
# ensemble mean, one for ensemble spread). The bulk of the consolidation
# (that is, variable, dimension, and attribute copying) is done by ncks,
# while this script keeps track of the input and output files and the
# required variables.
#
# REQUIREMENTS:
# * Python 3
# * NetCDF Operator (NCO) utilities
#
# REVISION HISTORY:
# 07 Feb 2018:  Eric Kemp (SSAI), first version.
# 23 Feb 2018:  Eric Kemp (SSAI), tweaked Noah invocation list.
# 01 Mar 2018:  Eric Kemp (SSAI), reorganized Noah invocation list.
# 15 Mar 2018:  Eric Kemp (SSAI), added 24-hour concatenation
# 28 Mar 2018:  Eric Kemp (SSAI), added JULES support.
# 30 Mar 2018:  Eric Kemp (SSAI), added several JULES variables.
# 03 Apr 2018:  Eric Kemp (SSAI), path to ncks now hardwired to better comply
#               with Air Force security requirements.
# 11 Apr 2018:  Eric Kemp (SSAI), add option to skip ensemble spread
# 04 Dec 2018:  Eric Kemp (SSAI), add mean 24hr Tair.
# 07 Nov 2019:  Eric Kemp (SSAI), removed Soiltype_inst for JULES.  Added
#               support for NoahMP
# 03 Dec 2019:  Eric Kemp (SSAI), added Greenness_inst for Noah and NoahMP.
#               Not included for JULES since that LSM doesn't use it.
# 09 Jan 2020:  Eric Kemp (SSAI), added Tair_f_min for JULES for 3hr.
# 05 Aug 2020:  Eric Kemp (SSAI), added Albedo_tavg and SmLiqFrac_inst.
# 25 Sep 2020:  Eric Kemp (SSAI), tweaked comments for Python version. Also
#               added path for NCKS on Koehr.
# 14 Oct 2020:  Eric Kemp (SSAI), updated NCKS path on Discover.
# 05 Feb 2021:  Eric Kemp (SSAI), added JULES multi-layer snow variables.
# 23 Mar 2021:  Eric Kemp (SSAI), revised JULES multi-layer snow variables
#               for PS41 physics.
# 05 Dec 2022:  Eric Kemp (SSAI), revised to improve pylint score.
# 24 Jan 2023:  Eric Kemp (SSAI), updated filenames.
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
_NCKS_PATH = "/usr/local/other/nco/4.8.1/bin/ncks" # On Discover

# Supported LIS LSMs
_LIS_LSMS = ["NOAH", "NOAHMP", "JULES"]

# The LVT invocations for Noah LSM output.  Each invocation handles a subset
# of the total variable list due to memory limitations.
_LVT_NOAH_INVOCATIONS_3HR = ['Albedo_tavg',
                             'AvgSurfT_inst', 'AvgSurfT_tavg',
                             'CanopInt_inst', 'Elevation_inst', 'Evap_tavg',
                             'Greenness_inst',
                             'LWdown_f_inst', 'LWdown_f_tavg',
                             'Landcover_inst', 'Landmask_inst', 'PotEvap_tavg',
                             'Psurf_f_inst', 'Psurf_f_tavg',
                             'Qair_f_inst', 'Qair_f_tavg',
                             'Qg_tavg', 'Qh_tavg', 'Qle_tavg',
                             'Qs_acc', 'Qsb_acc',
                             'RelSMC_inst', 'RHMin_inst',
                             'SWE_inst',
                             'SWdown_f_inst', 'SWdown_f_tavg',
                             'SmLiqFrac_inst', 'SnowDepth_inst',
                             'Snowcover_inst',
                             'SoilMoist_inst', 'SoilMoist_tavg',
                             'SoilTemp_inst', 'SoilTemp_tavg',
                             'Soiltype_inst',
                             'Tair_f_inst', 'Tair_f_max',
                             'Tair_f_tavg',
                             'TotalPrecip_acc', 'Wind_f_inst', 'Wind_f_tavg']

_LVT_NOAH_INVOCATIONS_24HR = ['Evap_tavg', 'LWdown_f_tavg', 'PotEvap_tavg',
                              'RHMin_inst',
                              'SoilMoist_tavg', 'SoilTemp_tavg',
                              'SWdown_f_tavg', 'Tair_f_max',
                              'Tair_f_tavg',
                              'TotalPrecip_acc', 'Wind_f_tavg']

# The 24-hr postprocessing should include the latest 3-hr snow depth and SWE.
_LVT_NOAH_INVOCATIONS_24HR_LATEST = ['SnowDepth_inst', 'SWE_inst']

# The LVT invocations for NoahMP LSM output.  Each invocation handles a subset
# of the total variable list due to memory limitations.
_LVT_NOAHMP_INVOCATIONS_3HR = ['Albedo_tavg',
                               'AvgSurfT_inst', 'AvgSurfT_tavg',
                               'CanopInt_inst', 'Elevation_inst', 'Evap_tavg',
                               'Greenness_inst',
                               'LWdown_f_inst', 'LWdown_f_tavg',
                               'Landcover_inst', 'Landmask_inst',
                               'Psurf_f_inst', 'Psurf_f_tavg',
                               'Qair_f_inst', 'Qair_f_tavg',
                               'Qg_tavg', 'Qh_tavg', 'Qle_tavg',
                               'Qs_acc', 'Qsb_acc',
                               'RelSMC_inst', 'RHMin_inst',
                               'SWE_inst',
                               'SWdown_f_inst', 'SWdown_f_tavg',
                               'SmLiqFrac_inst', 'SnowDepth_inst',
                               'Snowcover_inst',
                               'SoilMoist_inst', 'SoilMoist_tavg',
                               'SoilTemp_inst', 'SoilTemp_tavg',
                               'Soiltype_inst',
                               'Tair_f_inst', 'Tair_f_max',
                               'Tair_f_tavg',
                               'TotalPrecip_acc', 'Wind_f_inst', 'Wind_f_tavg']

_LVT_NOAHMP_INVOCATIONS_24HR = ['Evap_tavg', 'LWdown_f_tavg',
                                'RHMin_inst',
                                'SoilMoist_tavg', 'SoilTemp_tavg',
                                'SWdown_f_tavg', 'Tair_f_max',
                                'Tair_f_tavg',
                                'TotalPrecip_acc', 'Wind_f_tavg']

# The 24-hr postprocessing should include the latest 3-hr snow depth and SWE.
_LVT_NOAHMP_INVOCATIONS_24HR_LATEST = ['SnowDepth_inst', 'SWE_inst']

# The LVT invocations for JULES LSM output.
_LVT_JULES_INVOCATIONS_3HR = ['Albedo_tavg',
                              'AvgSurfT_inst', 'AvgSurfT_tavg',
                              'CanopInt_inst',
                              'Elevation_inst', 'Evap_tavg',
                              'LWdown_f_inst', 'LWdown_f_tavg',
                              'Landcover_inst', 'Landmask_inst',
                              'Psurf_f_inst', 'Psurf_f_tavg',
                              'Qair_f_inst', 'Qair_f_tavg',
                              'Qh_tavg', 'Qle_tavg',
                              'Qs_acc', 'Qsb_acc',
                              'RelSMC_inst', 'RHMin_inst',
                              'SWE_inst',
                              'SWdown_f_inst', 'SWdown_f_tavg',
                              'SnowDepth_inst',
                              'SmLiqFrac_inst',
                              'SoilMoist_inst', 'SoilMoist_tavg',
                              'SoilTemp_inst', 'SoilTemp_tavg',
                              'Tair_f_inst', 'Tair_f_max',
                              'Tair_f_tavg',
                              'TotalPrecip_acc', 'Wind_f_inst', 'Wind_f_tavg',
                              'ActSnowNL_inst', 'GrndSnow_inst',
                              'LayerSnowDensity_inst', 'LayerSnowDepth_inst',
                              'LayerSnowGrain_inst', 'SnowDensity_inst',
                              'SnowGrain_inst', 'SnowIce_inst',
                              'SnowLiq_inst',
                              'SnowTProf_inst', 'SurftSnow_inst']

# EMK for RECON
_LVT_JULES_INVOCATIONS_3HR = ["AvgSurfT_inst",
                              "SoilMoist_inst","SoilTemp_inst",
                              "PS41Snow_inst"]

# JULES PS41 snow variables are in a unique netCDF file.
_LVT_JULES_PS41_SNOW_3HR = ["SnowDepth_inst", "SWE_inst",
                            'ActSnowNL_inst', 'GrndSnow_inst',
                            'LayerSnowDensity_inst', 'LayerSnowDepth_inst',
                            'LayerSnowGrain_inst', 'SnowDensity_inst',
                            'SnowGrain_inst', 'SnowIce_inst',
                            'SnowLiq_inst',
                            'SnowTProf_inst', 'SurftSnow_inst']

_LVT_JULES_INVOCATIONS_24HR = ['Evap_tavg', 'LWdown_f_tavg',
                               'RHMin_inst',
                               'SoilMoist_tavg', 'SoilTemp_tavg',
                               'SWdown_f_tavg', 'Tair_f_max',
                               'Tair_f_tavg',
                               'TotalPrecip_acc', 'Wind_f_tavg']

# The 24-hr postprocessing should include the latest 3-hr snow depth and SWE.
_LVT_JULES_INVOCATIONS_24HR_LATEST = ['SnowDepth_inst', 'SWE_inst']

# The combined invocation directory for all supported LSMs.
_INVOCATIONS = {
    "NOAH_3HR": _LVT_NOAH_INVOCATIONS_3HR,
    "NOAH_24HR": _LVT_NOAH_INVOCATIONS_24HR,
    "NOAH_24HR_LATEST": _LVT_NOAH_INVOCATIONS_24HR_LATEST,
    "NOAHMP_3HR": _LVT_NOAHMP_INVOCATIONS_3HR,
    "NOAHMP_24HR": _LVT_NOAHMP_INVOCATIONS_24HR,
    "NOAHMP_24HR_LATEST": _LVT_NOAHMP_INVOCATIONS_24HR_LATEST,
    "JULES_3HR": _LVT_JULES_INVOCATIONS_3HR,
    "JULES_24HR": _LVT_JULES_INVOCATIONS_24HR,
    "JULES_24HR_LATEST": _LVT_JULES_INVOCATIONS_24HR_LATEST,
}

# The Noah variables handled by each LVT invocation.
_LIS_NOAH_VARIABLES_3HR = {}
for var in _LVT_NOAH_INVOCATIONS_3HR:
    if var == "RHMin_inst":
        _LIS_NOAH_VARIABLES_3HR[var] = [var, "Tair_f_min"]
    else:
        _LIS_NOAH_VARIABLES_3HR[var] = [var]

# The Noah variables handled by each LVT invocation.
_LIS_NOAH_VARIABLES_24HR = {}
for var in _LVT_NOAH_INVOCATIONS_24HR:
    if var == "RHMin_inst":
        _LIS_NOAH_VARIABLES_24HR[var] = [var, "Tair_f_min"]
    else:
        _LIS_NOAH_VARIABLES_24HR[var] = [var]

_LIS_NOAH_VARIABLES_24HR_LATEST = {}
for var in _LVT_NOAH_INVOCATIONS_24HR_LATEST:
    _LIS_NOAH_VARIABLES_24HR_LATEST[var] = [var]

# The NoahMP variables handled by each LVT invocation.
_LIS_NOAHMP_VARIABLES_3HR = {}
for var in _LVT_NOAHMP_INVOCATIONS_3HR:
    if var == "RHMin_inst":
        _LIS_NOAHMP_VARIABLES_3HR[var] = [var, "Tair_f_min"]
    else:
        _LIS_NOAHMP_VARIABLES_3HR[var] = [var]

# The NoahMP variables handled by each LVT invocation.
_LIS_NOAHMP_VARIABLES_24HR = {}
for var in _LVT_NOAHMP_INVOCATIONS_24HR:
    if var == "RHMin_inst":
        _LIS_NOAHMP_VARIABLES_24HR[var] = [var, "Tair_f_min"]
    else:
        _LIS_NOAHMP_VARIABLES_24HR[var] = [var]

_LIS_NOAHMP_VARIABLES_24HR_LATEST = {}
for var in _LVT_NOAHMP_INVOCATIONS_24HR_LATEST:
    _LIS_NOAHMP_VARIABLES_24HR_LATEST[var] = [var]

# The JULES variables handled by each LVT invocation.
_LIS_JULES_VARIABLES_3HR = {}
for var in _LVT_JULES_INVOCATIONS_3HR:
    if var == "RHMin_inst":
        _LIS_JULES_VARIABLES_3HR[var] = [var, "Tair_f_min"]
    else:
        _LIS_JULES_VARIABLES_3HR[var] = [var]

# The JULES variables handled by each LVT invocation.
_LIS_JULES_VARIABLES_24HR = {}
for var in _LVT_JULES_INVOCATIONS_24HR:
    if var == "RHMin_inst":
        _LIS_JULES_VARIABLES_24HR[var] = [var, "Tair_f_min"]
    else:
        _LIS_JULES_VARIABLES_24HR[var] = [var]

_LIS_JULES_VARIABLES_24HR_LATEST = {}
for var in _LVT_JULES_INVOCATIONS_24HR_LATEST:
    _LIS_JULES_VARIABLES_24HR_LATEST[var] = [var]

# Combined breakdown of all variables handled by LSM and LVT invocation.
_LIS_VARIABLES = {
    "NOAH_3HR": _LIS_NOAH_VARIABLES_3HR,
    "NOAH_24HR": _LIS_NOAH_VARIABLES_24HR,
    "NOAH_24HR_LATEST": _LIS_NOAH_VARIABLES_24HR_LATEST,
    "NOAHMP_3HR": _LIS_NOAHMP_VARIABLES_3HR,
    "NOAHMP_24HR": _LIS_NOAHMP_VARIABLES_24HR,
    "NOAHMP_24HR_LATEST": _LIS_NOAHMP_VARIABLES_24HR_LATEST,
    "JULES_3HR": _LIS_JULES_VARIABLES_3HR,
    "JULES_24HR": _LIS_JULES_VARIABLES_24HR,
    "JULES_24HR_LATEST": _LIS_JULES_VARIABLES_24HR_LATEST,
}

# These variables are processed by all invocations, and are either generated
# automatically or are read in from US Navy GOFS files.
_OTHER_VARIABLES = ["latitude", "longitude",
                    "time", "water_temp", "aice", "hi"]

#------------------------------------------------------------------------------
def _usage():
    """Print command line usage"""
    print(f"Usage: {sys.argv[0]} yyyymmddhh lsm period [--nospread]")
    print("   where:")
    print("           yyyymmddhh is valid year/month/day/hour in UTC")
    print("           lsm is name of land surface model used by LIS")
    print("           period is processing time length in hours (3 or 24)")
    print("           --nospread is optional flag to skip ensemble spread")

#------------------------------------------------------------------------------
def _read_cmd_args():
    """Read command line arguments"""
    # Check if argument count is correct
    if len(sys.argv) not in [4, 5]:
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
        validdt = datetime.datetime(year, month, day, hour)
    except ValueError:
        print("[ERR] Cannot process yyyymmddhh argument!")
        _usage()
        sys.exit(1)

    # Get lsm name
    lsm = None
    if sys.argv[2] in _LIS_LSMS:
        lsm = sys.argv[2]
    if lsm is None:
        print("[ERR] Invalid lsm selection!")
        print(f" lsm value is {sys.argv[2]}")
        text = " Supported lsms:"
        for lsm in _LIS_LSMS:
            text += f" {lsm}"
        print(text)
        sys.exit(1)

    # Get period
    period_options = [3, 24]
    period = None
    tmp_int = int(sys.argv[3])
    if tmp_int in period_options:
        period = tmp_int
    if period is None:
        print("[ERR] Invalid time period selected!")
        print(f"  Read in {sys.argv[3]}")
        print("  Only supports 3 or 24!")
        sys.exit(1)

    # Check if we are skipping ensemble spread
    skip_ens_spread = False
    if len(sys.argv) == 5:
        if sys.argv[4] == "--nospread":
            skip_ens_spread = True
        else:
            print(f"[ERR] Invalid argument {sys.argv[4]}")
            _usage()

    # See if ncks exists and is executable by current user (the script)
    # This used to be specified on the command line, but is now hardwired
    # to better comply with Air Force security requirements.
    ncks = _NCKS_PATH
    _check_ncks_path(ncks)

    return validdt, lsm, period, skip_ens_spread

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
def _get_nc_mean_files(validdt, lsm, period):
    """Collect netCDF mean files"""
    key = f"{lsm}_{period}HR"
    invocation_list = _INVOCATIONS[key]

    mean_nc_infiles = {}

    # Collect input files
    for invocation in invocation_list:

        path = f"OUTPUT/STATS.{invocation}.{period}hr"

        path += f"/PS.557WW_SC.U_DI.C_GP.LIS-{lsm}_GR.C0P09DEG_AR.GLOBAL_PA"
        if period == 24:
            path += ".LIS24_DD."
        else:
            path += ".LIS_DD."
        path += f"{validdt.year:04}{validdt.month:02}{validdt.day:02}_DT"
        path += f".{validdt.hour:02}00_DF"

        mean_path = path + ".nc"
        if not os.path.exists(mean_path):
            print(f"[ERR], {mean_path} does not exist!")
            sys.exit(1)

        mean_nc_infiles[invocation] = mean_path

    # Get output file
    path = f"OUTPUT/STATS_merged_{period}hr"
    if not os.path.exists(path):
        os.mkdir(path)
    path += f"/PS.557WW_SC.U_DI.C_GP.LIS-{lsm}_GR.C0P09DEG_AR.GLOBAL_PA"
    if period == 24:
        path += ".LIS24_DD."
    else:
        path += ".LIS_DD."
    path += f"{validdt.year:04}{validdt.month:02}{validdt.day:02}_DT"
    path += f".{validdt.hour:02}00_DF"

    mean_nc_outfile = path + ".nc"

    # All done
    return mean_nc_infiles, mean_nc_outfile

#------------------------------------------------------------------------------
def _get_nc_ssdev_files(validdt, lsm, period):
    """Collect netCDF ssdev files"""
    key = f"{lsm}_{period}HR"
    invocation_list = _INVOCATIONS[key]

    ssdev_nc_infiles = {}

    # Collect input files
    for invocation in invocation_list:
        path = f"OUTPUT/STATS.{invocation}.{period}hr"
        path += f"/PS.557WW_SC.U_DI.C_GP.LIS-{lsm}_GR.C0P09DEG_AR.GLOBAL_PA"
        if period == 24:
            path += ".LIS24-SSDEV_DD."
        else:
            path += ".SSDEV_DD."
        path += f"{validdt.year:04}{validdt.month:02}{validdt.day:02}_DT"
        path += f".{validdt.hour:02}00_DF"

        ssdev_path = path + ".nc"
        if not os.path.exists(ssdev_path):
            print(f"[ERR], {ssdev_path} does not exist!")
            sys.exit(1)

        ssdev_nc_infiles[invocation] = ssdev_path

    # Get output file
    path = f"OUTPUT/STATS_merged_{period}hr"
    if not os.path.exists(path):
        os.mkdir(path)
    path += f"/PS.557WW_SC.U_DI.C_GP.LIS-{lsm}_GR.C0P09DEG_AR.GLOBAL_PA"
    if period == 24:
        path += ".LIS24-SSDEV_DD."
    else:
        path += ".SSDEV_DD."
    path += f"{validdt.year:04}{validdt.month:02}{validdt.day:02}_DT"
    path += f".{validdt.hour:02}00_DF"

    ssdev_nc_outfile = path + ".nc"

    # All done
    return ssdev_nc_infiles, ssdev_nc_outfile

#------------------------------------------------------------------------------
def _get_nc_latest_files(validdt, lsm):
    """Collect netCDF latest files"""
    key = f"{lsm}_24HR_LATEST"
    invocation_list = _INVOCATIONS[key]

    latest_nc_infiles = {}

    # Collect input files
    for invocation in invocation_list:
        path = f"OUTPUT/STATS.{invocation}.3hr"# Always use 3hr processing
        path += f"/PS.557WW_SC.U_DI.C_GP.LIS-{lsm}_GR.C0P09DEG_AR.GLOBAL_PA"
        path += ".LIS_DD."  # Always use 3-hr output for latest fields
        path += f"{validdt.year:04}{validdt.month:02}{validdt.day:02}_DT"
        path += f".{validdt.hour:02}00_DF"

        latest_path = path + ".nc"
        if not os.path.exists(latest_path):
            print(f"[ERR], {latest_path} does not exist!")
            sys.exit(1)

        latest_nc_infiles[invocation] = latest_path

    # All done
    return latest_nc_infiles

#------------------------------------------------------------------------------
def _merge_nc_files(lsm, period, nc_infiles,
                   nc_outfile, latest_nc_infiles=None):
    """Use ncks to merge netCDF fields together"""

    ncks = _NCKS_PATH

    key = f"{lsm}_{period}HR"

    # Start with ensemble mean
    cmd = f"cp {nc_infiles[_INVOCATIONS[key][0]]} {nc_outfile}"
    print(cmd)
    err = subprocess.call(cmd, shell=True)
    if err != 0:
        print("[ERR] Problem with cp!")
        sys.exit(1)

    invocations = _INVOCATIONS[key][1:]
    for invocation in invocations:

        # Special handing of JULES PS41 snow variables
        if invocation == "PS41Snow_inst":
            variables = _LVT_JULES_PS41_SNOW_3HR[:]
        else:
            variables = _LIS_VARIABLES[key][invocation]

        for variable in variables:
            cmd = f"{ncks} -A -v {variable} {nc_infiles[invocation]} "
            cmd += f"{nc_outfile}"
            print(cmd)
            err = subprocess.call(cmd, shell=True)
            if err != 0:
                print("[ERR] Problem with ncks!")
                sys.exit(1)

    # For 24-hr postprocessing, we also must concatenate several 3-hr fields
    if latest_nc_infiles is not None:
        key = f"{lsm}_24HR_LATEST"
        invocations = _INVOCATIONS[key][:]
        for invocation in invocations:
            variables = _LIS_VARIABLES[key][invocation]
            for variable in variables:
                cmd = f"{ncks} -A -v "
                cmd += f"{variable} {latest_nc_infiles[invocation]} "
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
    validdt, lsm, period, skip_ens_spread = _read_cmd_args()

    # Collect netCDF files
    (mean_nc_infiles, mean_nc_outfile) = \
        _get_nc_mean_files(validdt, lsm, period)

    # 3-hr postprocessing includes ensemble spread files
    if period == 3 and not skip_ens_spread:
        (ssdev_nc_infiles, ssdev_nc_outfile) = \
            _get_nc_ssdev_files(validdt, lsm, period)

    # 24-hr postprocessing includes several latest 3-hr fields
    if period == 24:
        latest_nc_infiles = _get_nc_latest_files(validdt, lsm)

    # Merge the input netCDF files together
    if period == 3:
        _merge_nc_files(lsm, period, mean_nc_infiles, mean_nc_outfile)
        if not skip_ens_spread:
            _merge_nc_files(lsm, period,
                            ssdev_nc_infiles, ssdev_nc_outfile)
    else:
        # 24-hr processing
        _merge_nc_files(lsm, period, mean_nc_infiles,
                        mean_nc_outfile, latest_nc_infiles)

if __name__ == "__main__":
    _main()
