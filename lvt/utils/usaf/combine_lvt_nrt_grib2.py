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
# SCRIPT: combine_lvt_nrt_grib2.py
#
# PURPOSE:  Runs Linux 'cat' command to combine multiple GRIB2 output files
# produced by LVT.
#
# REVISION HISTORY:
# 27 Jan 2026:  Eric Kemp (SSAI), based on cat_lvt_grib2.py.  Updates for
#   merged NRT/NRT-Streamflow.  NoahMP now includes PotEvap.  No JULES
#   support.  No HYMAP support.
# 18 Aug 2026:  Eric Kemp (SSAI), further logic updates, and renamed to
#   better describe function.
# 22 Sep 2026:  Eric Kemp (SSAI), replaced LSM and Routing command line
#   arguments with single generating_process argument.  Revised _INVOCATION
#   dictionary.  Revised order of command line arguments to simplify
#   manual processing of multiple valid times.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import subprocess
import sys

#------------------------------------------------------------------------------

# Supported generating processes
_GENERATING_PROCESSES = ["LIS-NRT-NOAH", "LIS-NRT-NOAHMP",
                         "LIS-NRT-NOAH-RAPID", "LIS-NRT-NOAHMP-RAPID"]

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
# NOTE:  It could be argued that this can be consolidated with the NOAH lists.
# However, NOAHMP may have additional variables in the future, so we will keep
# them separate.
_LVT_NOAHMP_INVOCATIONS_3HR = ['Albedo_tavg',
                               'AvgSurfT_inst', 'AvgSurfT_tavg',
                               'CanopInt_inst', 'Elevation_inst', 'Evap_tavg',
                               'Greenness_inst',
                               'LWdown_f_inst', 'LWdown_f_tavg',
                               'Landcover_inst', 'Landmask_inst',
                               'PotEvap_tavg',
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

_LVT_NOAHMP_INVOCATIONS_24HR = ['Evap_tavg', 'LWdown_f_tavg', 'PotEvap_tavg',
                                'RHMin_inst',
                                'SoilMoist_tavg', 'SoilTemp_tavg',
                                'SWdown_f_tavg', 'Tair_f_max',
                                'Tair_f_tavg',
                                'TotalPrecip_acc', 'Wind_f_tavg']

# The 24-hr postprocessing should include the latest 3-hr snow depth and SWE.
_LVT_NOAHMP_INVOCATIONS_24HR_LATEST = ['SnowDepth_inst', 'SWE_inst']

# The combined invocation dictionary for all supported generating processes.
_INVOCATIONS = {
    "LIS-NRT-NOAH_3HR": _LVT_NOAH_INVOCATIONS_3HR,
    "LIS-NRT-NOAH-RAPID_3HR": _LVT_NOAH_INVOCATIONS_3HR,
    "LIS-NRT-NOAH_24HR": _LVT_NOAH_INVOCATIONS_24HR,
    "LIS-NRT-NOAH-RAPID_24HR": _LVT_NOAH_INVOCATIONS_24HR,
    "LIS-NRT-NOAH_24HR_LATEST": _LVT_NOAH_INVOCATIONS_24HR_LATEST,
    "LIS-NRT-NOAH-RAPID_24HR_LATEST": _LVT_NOAH_INVOCATIONS_24HR_LATEST,
    "LIS-NRT-NOAHMP_3HR": _LVT_NOAHMP_INVOCATIONS_3HR,
    "LIS-NRT-NOAHMP-RAPID_3HR": _LVT_NOAHMP_INVOCATIONS_3HR,
    "LIS-NRT-NOAHMP_24HR": _LVT_NOAHMP_INVOCATIONS_24HR,
    "LIS-NRT-NOAHMP-RAPID_24HR": _LVT_NOAHMP_INVOCATIONS_24HR,
    "LIS-NRT-NOAHMP_24HR_LATEST": _LVT_NOAHMP_INVOCATIONS_24HR_LATEST,
    "LIS-NRT-NOAHMP-RAPID_24HR_LATEST": _LVT_NOAHMP_INVOCATIONS_24HR_LATEST,
}

# -----------------------------------------------------------------------------
def _usage():
    """Print command line usage"""
    print(f"Usage: {sys.argv[0]} generating_process period yyyymmddhh fhh" + \
          " [--nospread]")
    print("   where:")
    print("        generating_process is GP section of LIS output filename")
    print("        period is time period (hours) for postprocessing (3 or 24)")
    print("        yyyymmddhh is LIS start year/month/day/hour in UTC")
    print("        fhh is LIS forecast hour (hh) in UTC")
    print("        --nospread is optional flag to skip ensemble spread")

# -----------------------------------------------------------------------------
def _read_cmd_args():
    """Read command line arguments"""
    # Check if argument count is correct
    if len(sys.argv) not in [5, 6]:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Get generating_process
    generating_process = None
    if sys.argv[1] in _GENERATING_PROCESSES:
        generating_process = sys.argv[1]
    if generating_process is None:
        print("[ERR] Invalid generating_process selection!")
        print(f" generated_processes value is {sys.argv[1]}")
        text = " Supported generated_processes:"
        for generating_process in _GENERATING_PROCESSES:
            text += f" {generating_process}"
        print(text)
        sys.exit(1)

    # Get processing period
    period_options = [3, 24]
    period = None
    tmp_int = int(sys.argv[2])
    if tmp_int in period_options:
        period = tmp_int
    if period is None:
        print("[ERR] Invalid period selection!")
        print(f" period value is {sys.argv[2]}")
        print(" Supported time periods are: 3 and 24")
        sys.exit(1)

    # Get the start date, the start hour (cycle hour), and the forecast hour
    yyyymmddhh = sys.argv[3]
    fhh = sys.argv[4]
    try:
        year = int(yyyymmddhh[0:4])
        month = int(yyyymmddhh[4:6])
        day = int(yyyymmddhh[6:8])
        hour = int(yyyymmddhh[8:10])
        startdt = datetime.datetime(year, month, day, hour)
        forecast_hour = int(fhh)
        validdt = startdt + datetime.timedelta(hours=forecast_hour)
    except ValueError:
        print("[ERR] Cannot process valid time arguments!")
        _usage()
        sys.exit(1)

    # Check if ensemble spread should be skipped
    skip_ens_spread = False
    if len(sys.argv) == 6:
        if sys.argv[5] == "--nospread":
            skip_ens_spread = True
        else:
            print(f"[ERR] Invalid argument {sys.argv[5]}")
            _usage()

    return validdt, forecast_hour, generating_process, period, skip_ens_spread

# -----------------------------------------------------------------------------
def _get_gr2_mean_files(validdt, forecast_hour, generating_process, period):
    """Collect GRIB2 mean files"""
    key = f"{generating_process}_{period}HR"
    invocation_list = _INVOCATIONS[key]

    mean_gr2_infiles = {}

    startdt = validdt - datetime.timedelta(hours=forecast_hour)

    basename = "PS.557WW"
    basename += "_SC.U"
    basename += "_DI.C"
    basename += f"_GP.{generating_process}"
    basename += "_GR.C0P09DEG"
    basename += "_AR.GLOBAL"
    if period == 24:
        basename += "_PA.LIS24"
    else:
        basename += "_PA.LIS"
    basename += f"_DD.{startdt.year:04}{startdt.month:02}{startdt.day:02}"
    basename += f"_CY.{startdt.hour:02}"
    if period == 24:
        basename += "_FH.012"
    else:
        basename += f"_FH.{forecast_hour:03}"
    basename += "_DF.GR2"

    # Collect input files
    for invocation in invocation_list:
        topdir = f"OUTPUT/STATS.{invocation}.{period}hr"
        mean_path = f"{topdir}/{basename}"
        if not os.path.exists(mean_path):
            print(f"[ERR], {mean_path} does not exist!")
            sys.exit(1)
        mean_gr2_infiles[invocation] = mean_path

    # Get output files
    topdir = f"OUTPUT/STATS_merged_{period}hr"
    if not os.path.exists(topdir):
        os.mkdir(topdir)
    mean_gr2_outfile = f"{topdir}/{basename}"

    # All done
    return mean_gr2_infiles, mean_gr2_outfile

# -----------------------------------------------------------------------------
def _get_gr2_ssdev_files(validdt, forecast_hour, generating_process, period):
    """Collect GRIB2 ssdev files"""
    key = f"{generating_process}_{period}HR"
    invocation_list = _INVOCATIONS[key]

    ssdev_gr2_infiles = {}

    startdt = validdt - datetime.timedelta(hours=forecast_hour)

    basename = "PS.557WW"
    basename += "_SC.U"
    basename += "_DI.C"
    basename += f"_GP.{generating_process}"
    basename += "_GR.C0P09DEG"
    basename += "_AR.GLOBAL"
    if period == 24:
        basename += "_PA.LIS24-SSDEV"
    else:
        basename += "_PA.SSDEV"
    basename += f"_DD.{startdt.year:04}{startdt.month:02}{startdt.day:02}"
    basename += f"_CY.{startdt.hour:02}"
    if period == 24:
        basename += "_FH.012"
    else:
        basename += f"_FH.{forecast_hour:03}"
    basename += "_DF.GR2"

    # Collect input files
    for invocation in invocation_list:
        topdir = f"OUTPUT/STATS.{invocation}.{period}hr"
        ssdev_path = f"{topdir}/{basename}"
        if not os.path.exists(ssdev_path):
            print(f"[ERR], {ssdev_path} does not exist!")
            sys.exit(1)
        ssdev_gr2_infiles[invocation] = ssdev_path

    # Get output file
    topdir = f"OUTPUT/STATS_merged_{period}hr"
    if not os.path.exists(topdir):
        os.mkdir(topdir)

    ssdev_gr2_outfile = f"{topdir}/{basename}"

    # All done
    return ssdev_gr2_infiles, ssdev_gr2_outfile

# -----------------------------------------------------------------------------
def _get_gr2_latest_files(validdt, forecast_hour, generating_process):
    """Collect GRIB2 latest files"""

    key = f"{generating_process}_24HR_LATEST"
    invocation_list = _INVOCATIONS[key]

    latest_gr2_infiles = {}

    startdt = validdt - datetime.timedelta(hours=forecast_hour)

    basename = "PS.557WW"
    basename += "_SC.U"
    basename += "_DI.C"
    basename += f"_GP.{generating_process}"
    basename += "_GR.C0P09DEG"
    basename += "_AR.GLOBAL"
    basename += "_PA.LIS"
    basename += f"_DD.{startdt.year:04}{startdt.month:02}{startdt.day:02}"
    basename += f"_CY.{startdt.hour:02}"
    basename += f"_FH.{forecast_hour:03}"
    basename += "_DF.GR2"

    # Collect input files
    for invocation in invocation_list:
        topdir = f"OUTPUT/STATS.{invocation}.3hr" # Always use 3hr processing
        latest_path = f"{topdir}/{basename}"
        if not os.path.exists(latest_path):
            print(f"[ERR], {latest_path} does not exist!")
            sys.exit(1)
        latest_gr2_infiles[invocation] = latest_path

    # All done
    return latest_gr2_infiles

# -----------------------------------------------------------------------------
def _merge_gr2_files(generating_process, period, gr2_infiles, gr2_outfile,
                     latest_gr2_infiles=None):
    """Use cat to merge GRIB2 fields together"""
    key = f"{generating_process}_{period}HR"
    invocations = _INVOCATIONS[key][0:]
    cmd = "cat"
    for invocation in invocations:
        cmd += f" {gr2_infiles[invocation]}"
    # For 24-hr postprocessing, we also must concatenate several 3-hr fields
    if latest_gr2_infiles is not None:
        key = f"{generating_process}_24HR_LATEST"
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
    validdt, forecast_hour, generating_process, period, skip_ens_spread = \
        _read_cmd_args()

    # Collect GRIB2 files
    (mean_gr2_infiles, mean_gr2_outfile) = \
        _get_gr2_mean_files(validdt, forecast_hour, generating_process, period)
    # 3-hr postprocessing includes ensemble spread files
    if period == 3 and not skip_ens_spread:
        (ssdev_gr2_infiles, ssdev_gr2_outfile) = \
            _get_gr2_ssdev_files(validdt, forecast_hour, generating_process,
                                 period)
    # 24-hr postprocessing includes several latest 3-hr fields
    if period == 24:
        latest_gr2_infiles = _get_gr2_latest_files(validdt, forecast_hour,
                                                   generating_process)

    # Merge the input GRIB2 files together
    if period == 3:
        _merge_gr2_files(generating_process, period, mean_gr2_infiles, \
                         mean_gr2_outfile)
        if not skip_ens_spread:
            _merge_gr2_files(generating_process, period, ssdev_gr2_infiles, \
                             ssdev_gr2_outfile)
    else:
        # 24-hr processing
        _merge_gr2_files(generating_process, period, mean_gr2_infiles,
                         mean_gr2_outfile, latest_gr2_infiles)

if __name__ == "__main__":
    _main()
