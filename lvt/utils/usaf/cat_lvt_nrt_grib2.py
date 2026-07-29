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
# SCRIPT: cat_lvt_nrt_grib2.py
#
# PURPOSE:  Runs Linux 'cat' command to combine multiple GRIB2 output files
# produced by LVT.
#
# REVISION HISTORY:
# 27 Jan 2026:  Eric Kemp (SSAI), based on cat_lvt_grib2.py.  Updates for
#   merged NRT/NRT-Streamflow.  NoahMP now includes PotEvap.  No JULES
#   support.  No HYMAP support.
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
_LIS_LSMS = ["NOAH", "NOAHMP"]

# Supported LIS Routing models
_LIS_ROUTINGS = ["RAPID"]

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


# The combined invocation directory for all supported LSMs and Routing models.
_INVOCATIONS = {
    "NOAH_RAPID_3HR": _LVT_NOAH_INVOCATIONS_3HR,
    "NOAH_RAPID_24HR": _LVT_NOAH_INVOCATIONS_24HR,
    "NOAH_RAPID_24HR_LATEST": _LVT_NOAH_INVOCATIONS_24HR_LATEST,
    "NOAHMP_RAPID_3HR": _LVT_NOAHMP_INVOCATIONS_3HR,
    "NOAHMP_RAPID_24HR": _LVT_NOAHMP_INVOCATIONS_24HR,
    "NOAHMP_RAPID_24HR_LATEST": _LVT_NOAHMP_INVOCATIONS_24HR_LATEST,
}

# -----------------------------------------------------------------------------
def _usage():
    """Print command line usage"""
    print(f"Usage: {sys.argv[0]} yyyymmdd chh fhh lsm routing period " + \
          "[--nospread]")
    print("   where:")
    print("        yyyymmdd is LIS start year/month/day in UTC")
    print("        chh is LIS start hour (hh) in UTC")
    print("        fhh is LIS forecast hour (hh) in UTC")
    print("        lsm is name of land surface model used by LIS")
    print("        routing is name of routing model used by LIS")
    print("        period is time period (hours) for postprocessing (3 or 24)")
    print("        --nospread is optional flag to skip ensemble spread")

# -----------------------------------------------------------------------------
def _read_cmd_args():
    """Read command line arguments"""
    # Check if argument count is correct
    if len(sys.argv) not in [7, 8]:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Get the start date, the start hour (cycle hour), and the forecast hour
    yyyymmdd = sys.argv[1]
    chh = sys.argv[2]
    fhh = sys.argv[3]
    try:
        year = int(yyyymmdd[0:4])
        month = int(yyyymmdd[4:6])
        day = int(yyyymmdd[6:8])
        hour = int(chh)
        startdt = datetime.datetime(year, month, day, hour)
        forecast_hour = int(fhh)
        validdt = startdt + datetime.timedelta(hours=forecast_hour)
    except ValueError:
        print("[ERR] Cannot process valid time arguments!")
        _usage()
        sys.exit(1)

    # Get lsm name
    lsm = None
    if sys.argv[4] in _LIS_LSMS:
        lsm = sys.argv[4]
    if lsm is None:
        print("[ERR] Invalid lsm selection!")
        print(f" lsm value is {sys.argv[4]}")
        text = " Supported lsms:"
        for lsm in _LIS_LSMS:
            text += f" {lsm}"
        print(text)
        sys.exit(1)

    # Get routing name
    routing = None
    if sys.argv[5] in _LIS_ROUTINGS:
        routing = sys.argv[5]
    if routing is None:
        print("[ERR] Invalid routing selection!")
        print(f" routing value is {sys.argv[5]}")
        text = " Supported routing models:"
        for routing in _LIS_ROUTINGS:
            text += f" {routing}"
        print(text)
        sys.exit(1)

    # Get processing hour
    period_options = [3, 24]
    period = None
    tmp_int = int(sys.argv[6])
    if tmp_int in period_options:
        period = tmp_int
    if period is None:
        print("[ERR] Invalid period selection!")
        print(f" period value is {sys.argv[6]}")
        print(" Supported time periods are: 3 and 24")
        sys.exit(1)

    # Check if ensemble spread should be skipped
    skip_ens_spread = False
    if len(sys.argv) == 8:
        if sys.argv[7] == "--nospread":
            skip_ens_spread = True
        else:
            print(f"[ERR] Invalid argument {sys.argv[7]}")
            _usage()

    return validdt, forecast_hour, lsm, routing, period, skip_ens_spread

# -----------------------------------------------------------------------------
def _get_gr2_mean_files(validdt, forecast_hour, lsm, routing, period):
    """Collect GRIB2 mean files"""
    key = f"{lsm}_{routing}_{period}HR"
    invocation_list = _INVOCATIONS[key]

    mean_gr2_infiles = {}

    startdt = validdt - datetime.timedelta(hours=forecast_hour)

    basename = "PS.557WW"
    basename += "_SC.U"
    basename += "_DI.C"
    basename += f"_GP.LIS-NRT-{lsm}-{routing}"
    basename += "_GR.C0P09DEG"
    basename += "_AR.GLOBAL"
    if period == 24:
        basename += "_PA.LIS24"
    else:
        basename += "_PA.LIS"
    basename += f"_DD.{startdt.year:04}{startdt.month:02}{startdt.day:02}"
    basename += f"_CY.{startdt.hour:02}"
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
def _get_gr2_ssdev_files(validdt, forecast_hour, lsm, routing, period):
    """Collect GRIB2 ssdev files"""
    key = f"{lsm}_{routing}_{period}HR"
    invocation_list = _INVOCATIONS[key]

    ssdev_gr2_infiles = {}

    startdt = validdt - datetime.timedelta(hours=forecast_hour)

    basename = "PS.557WW"
    basename += "_SC.U"
    basename += "_DI.C"
    basename += f"_GP.LIS-NRT-{lsm}-{routing}"
    basename += "_GR.C0P09DEG"
    basename += "_AR.GLOBAL"
    if period == 24:
        basename += "_PA.LIS24-SSDEV"
    else:
        basename += "_PA.SSDEV"
    basename += f"_DD.{startdt.year:04}{startdt.month:02}{startdt.day:02}"
    basename += f"_CY.{startdt.hour:02}"
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
def _get_gr2_latest_files(validdt, forecast_hour, lsm, routing):
    """Collect GRIB2 latest files"""

    key = f"{lsm}_{routing}_24HR_LATEST"
    invocation_list = _INVOCATIONS[key]

    latest_gr2_infiles = {}

    startdt = validdt - datetime.timedelta(hours=forecast_hour)

    basename = "PS.557WW"
    basename += "_SC.U"
    basename += "_DI.C"
    basename += f"_GP.LIS-NRT-{lsm}-{routing}"
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
def _merge_gr2_files(lsm, routing, period, gr2_infiles, gr2_outfile,
                    latest_gr2_infiles=None):
    """Use cat to merge GRIB2 fields together"""
    key = f"{lsm}_{routing}_{period}HR"
    invocations = _INVOCATIONS[key][0:]
    cmd = "cat"
    for invocation in invocations:
        cmd += f" {gr2_infiles[invocation]}"
    # For 24-hr postprocessing, we also must concatenate several 3-hr fields
    if latest_gr2_infiles is not None:
        key = f"{lsm}_{routing}_24HR_LATEST"
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
    validdt, forecast_hour, lsm, routing, period, skip_ens_spread = \
        _read_cmd_args()

    # Collect GRIB2 files
    (mean_gr2_infiles, mean_gr2_outfile) = \
        _get_gr2_mean_files(validdt, forecast_hour, lsm, routing, period)
    # 3-hr postprocessing includes ensemble spread files
    if period == 3 and not skip_ens_spread:
        (ssdev_gr2_infiles, ssdev_gr2_outfile) = \
            _get_gr2_ssdev_files(validdt, forecast_hour, lsm, routing, period)
    # 24-hr postprocessing includes several latest 3-hr fields
    if period == 24:
        latest_gr2_infiles = _get_gr2_latest_files(validdt, forecast_hour,
                                                   lsm, routing)

    # Merge the input GRIB2 files together
    if period == 3:
        _merge_gr2_files(lsm, routing, period, mean_gr2_infiles, \
                         mean_gr2_outfile)
        if not skip_ens_spread:
            _merge_gr2_files(lsm, routing, period, ssdev_gr2_infiles, \
                             ssdev_gr2_outfile)
    else:
        # 24-hr processing
        _merge_gr2_files(lsm, routing, period, mean_gr2_infiles,
                         mean_gr2_outfile, latest_gr2_infiles)

if __name__ == "__main__":
    _main()
