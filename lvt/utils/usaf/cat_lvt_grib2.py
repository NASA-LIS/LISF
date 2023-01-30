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
#------------------------------------------------------------------------------
#
# SCRIPT: cat_lvt_grib2.py
#
# PURPOSE:  Runs Linux 'cat' command to combine multiple GRIB2 output files
# produced by LVT.
#
# REVISION HISTORY:
# 02 Mar 2018:  Eric Kemp (SSAI), first version.
# 13 Mar 2018:  Eric Kemp (SSAI), added 24-hr processing.
# 14 Mar 2018:  Eric Kemp (SSAI), added appending latest 3-hr SnowDepth_inst
#               and SWE_inst.
# 28 Mar 2018:  Eric Kemp (SSAI), added JULES support.  Updated 24-hr file
#               naming convention.
# 30 Mar 2018:  Eric Kemp (SSAI), added more JULES variables.
# 11 Apr 2018:  Eric Kemp (SSAI), added optional flag to skip ensemble spreads.
# 16 Nov 2018:  Eric Kemp (SSAI), added Greenness_inst for 3hr.
# 19 Nov 2018:  Eric Kemp (SSAI), added Tair_tavg for 24hr.
# 07 Nov 2019:  Eric Kemp (SSAI), removed Soiltype_inst and Greenness_inst
#               for JULES.  Added support for NoahMP.
# 05 Aug 2020:  Eric Kemp (SSAI), added Albedo_tavg and SmLiqFrac_inst for
#               JULES.
# 05 Dec 2022:  Eric Kemp (SSAI), updates to improve pylint score.
# 24 Jan 2023:  Eric Kemp (SSAI), updates to GRIB file names.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import subprocess
import sys

#------------------------------------------------------------------------------

# Supported LIS LSMs
_LIS_LSMS = ["NOAH", "NOAHMP", "JULES"]

# The LVT invocations for Noah LSM output.  Each invocation handles a subset
# of the total variable list due to memory limitations.
# EXCEPTION:  RHMin_inst must be processed with Tair_f_min, so Tair_f_min
# is not included in the list below.
_LVT_NOAH_INVOCATIONS_3HR = ['Albedo_tavg', 'AvgSurfT_inst', 'AvgSurfT_tavg',
                             'CanopInt_inst', 'Elevation_inst', 'Evap_tavg',
                             'Greenness_inst',
                             'LWdown_f_inst', 'LWdown_f_tavg',
                             'Landcover_inst', 'Landmask_inst', 'PotEvap_tavg',
                             'Psurf_f_inst', 'Psurf_f_tavg',
                             'Qair_f_inst', 'Qair_f_tavg',
                             'Qg_tavg', 'Qh_tavg', 'Qle_tavg', 'Qs_acc',
                             'Qsb_acc', 'RHMin_inst', 'RelSMC_inst',
                             'SWE_inst', 'SWdown_f_inst', 'SWdown_f_tavg',
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

# The LVT invocation for NOAHMP LSM output.
_LVT_NOAHMP_INVOCATIONS_3HR = ['Albedo_tavg',
                               'AvgSurfT_inst', 'AvgSurfT_tavg',
                               'CanopInt_inst', 'Elevation_inst', 'Evap_tavg',
                               'Greenness_inst',
                               'LWdown_f_inst', 'LWdown_f_tavg',
                               'Landcover_inst', 'Landmask_inst',
                               'Psurf_f_inst', 'Psurf_f_tavg',
                               'Qair_f_inst', 'Qair_f_tavg',
                               'Qg_tavg', 'Qh_tavg', 'Qle_tavg', 'Qs_acc',
                               'Qsb_acc', 'RHMin_inst', 'RelSMC_inst',
                               'SWE_inst', 'SWdown_f_inst', 'SWdown_f_tavg',
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
                              'RHMin_inst', 'RelSMC_inst',
                              'SWE_inst',
                              'SWdown_f_inst', 'SWdown_f_tavg',
                              'SmLiqFrac_inst',
                              'SnowDepth_inst',
                              'SoilMoist_inst', 'SoilMoist_tavg',
                              'SoilTemp_inst', 'SoilTemp_tavg',
                              'Tair_f_inst', 'Tair_f_max',
                              'Tair_f_tavg',
                              'TotalPrecip_acc', 'Wind_f_inst', 'Wind_f_tavg']

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

# -----------------------------------------------------------------------------
def _usage():
    """Print command line usage"""
    print(f"Usage: {sys.argv[0]} yyyymmddhh lsm period [--nospread]")
    print("   where:")
    print("        yyyymmddhh is valid year/month/day/hour in UTC")
    print("        lsm is name of land surface model used by LIS")
    print("        period is time period (hours) for postprocessing (3 or 24)")
    print("        --nospread is optional flag to skip ensemble spread")

# -----------------------------------------------------------------------------
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

    # Get processing hour
    period_options = [3, 24]
    period = None
    tmp_int = int(sys.argv[3])
    if tmp_int in period_options:
        period = tmp_int
    if period is None:
        print("[ERR] Invalid period selection!")
        print(f" period value is {sys.argv[3]}")
        print(" Supported time periods are: 3 and 24")
        sys.exit(1)

    # Check if ensemble spread should be skipped
    skip_ens_spread = False
    if len(sys.argv) == 5:
        if sys.argv[4] == "--nospread":
            skip_ens_spread = True
        else:
            print(f"[ERR] Invalid argument {sys.argv[4]}")
            _usage()

    return validdt, lsm, period, skip_ens_spread

# -----------------------------------------------------------------------------
def _get_gr2_mean_files(validdt, lsm, period):
    """Collect GRIB2 mean files"""
    key = f"{lsm}_{period}HR"
    invocation_list = _INVOCATIONS[key]

    mean_gr2_infiles = {}

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

        mean_path = path + ".GR2"
        if not os.path.exists(mean_path):
            print(f"[ERR], {mean_path} does not exist!")
            sys.exit(1)
        mean_gr2_infiles[invocation] = mean_path

    # Get output files
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

    mean_gr2_outfile = path + ".GR2"

    # All done
    return mean_gr2_infiles, mean_gr2_outfile

# -----------------------------------------------------------------------------
def _get_gr2_ssdev_files(validdt, lsm, period):
    """Collect GRIB2 ssdev files"""
    key = f"{lsm}_{period}HR"
    invocation_list = _INVOCATIONS[key]

    ssdev_gr2_infiles = {}

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

        ssdev_path = path + ".GR2"
        if not os.path.exists(ssdev_path):
            print(f"[ERR], {ssdev_path} does not exist!")
            sys.exit(1)
        ssdev_gr2_infiles[invocation] = ssdev_path

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

    ssdev_gr2_outfile = path + ".GR2"

    # All done
    return ssdev_gr2_infiles, ssdev_gr2_outfile

# -----------------------------------------------------------------------------
def _get_gr2_latest_files(validdt, lsm):
    """Collect GRIB2 latest files"""

    key = f"{lsm}_24HR_LATEST"
    invocation_list = _INVOCATIONS[key]

    latest_gr2_infiles = {}

    # Collect input files
    for invocation in invocation_list:
        path = f"OUTPUT/STATS.{invocation}.3hr" # Always use 3hr processing
        path += f"/PS.557WW_SC.U_DI.C_GP.LIS-{lsm}_GR.C0P09DEG_AR.GLOBAL_PA"
        path += ".LIS_DD."
        path += f"{validdt.year:04}{validdt.month:02}{validdt.day:02}_DT"
        path += f".{validdt.hour:02}00_DF"

        latest_path = path + ".GR2"
        if not os.path.exists(latest_path):
            print(f"[ERR], {latest_path} does not exist!")
            sys.exit(1)
        latest_gr2_infiles[invocation] = latest_path

    # All done
    return latest_gr2_infiles

# -----------------------------------------------------------------------------
def _merge_gr2_files(lsm, period, gr2_infiles, gr2_outfile,
                    latest_gr2_infiles=None):
    """Use cat to merge GRIB2 fields together"""
    key = f"{lsm}_{period}HR"
    invocations = _INVOCATIONS[key][0:]
    cmd = "cat"
    for invocation in invocations:
        cmd += f" {gr2_infiles[invocation]}"
    # For 24-hr postprocessing, we also must concatenate several 3-hr fields
    if latest_gr2_infiles is not None:
        key = f"{lsm}_24HR_LATEST"
        invocations = _INVOCATIONS[key][:]
        for invocation in invocations:
            cmd += f" {latest_gr2_infiles[invocation]}"
    cmd += f" > {gr2_outfile}"

    print(cmd)
    err = subprocess.call(cmd, shell=True)
    if err != 0:
        print("[ERR] Problem with cat!")
        sys.exit(1)

# -----------------------------------------------------------------------------
# Main Driver.

def _main():
    """Main driver"""
    # Process command line arguments
    validdt, lsm, period, skip_ens_spread = _read_cmd_args()

    # Collect GRIB2 files
    (mean_gr2_infiles, mean_gr2_outfile) = \
        _get_gr2_mean_files(validdt, lsm, period)
    # 3-hr postprocessing includes ensemble spread files
    if period == 3 and not skip_ens_spread:
        (ssdev_gr2_infiles, ssdev_gr2_outfile) = \
            _get_gr2_ssdev_files(validdt, lsm, period)
    # 24-hr postprocessing includes several latest 3-hr fields
    if period == 24:
        latest_gr2_infiles = _get_gr2_latest_files(validdt, lsm)

    # Merge the input GRIB2 files together
    if period == 3:
        _merge_gr2_files(lsm, period, mean_gr2_infiles, mean_gr2_outfile)
        if not skip_ens_spread:
            _merge_gr2_files(lsm, period, ssdev_gr2_infiles, ssdev_gr2_outfile)
    else:
        # 24-hr processing
        _merge_gr2_files(lsm, period, mean_gr2_infiles,
                         mean_gr2_outfile, latest_gr2_infiles)

if __name__ == "__main__":
    _main()
