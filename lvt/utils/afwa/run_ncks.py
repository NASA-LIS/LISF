#!/usr/bin/env python
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
# * Python 2.6 or 2.7
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
#
#------------------------------------------------------------------------------

# Standard modules
import datetime
import os
import subprocess
import sys

#------------------------------------------------------------------------------

# Path to NCO ncks program
_NCKS_PATH = "/app/nco/4.5.2-gnu/bin/ncks" # On Conrad
#_NCKS_PATH = "/usr/local/other/SLES11.1/nco/4.4.4/intel-12.1.0.233/bin/ncks"

# Supported LIS LSMs
_LIS_LSMS = ["NOAH","JULES"] 

# The LVT invocations for Noah LSM output.  Each invocation handles a subset
# of the total variable list due to memory limitations.  
_LVT_NOAH_INVOCATIONS_3HR = ['Albedo_tavg', 
                             'AvgSurfT_inst', 'AvgSurfT_tavg', 
                             'CanopInt_inst', 'Elevation_inst', 'Evap_tavg', 
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
                              'SWdown_f_tavg','Tair_f_max',
                              'Tair_f_tavg',
                              'TotalPrecip_acc','Wind_f_tavg']

# The 24-hr postprocessing should include the latest 3-hr snow depth and SWE.
_LVT_NOAH_INVOCATIONS_24HR_LATEST = ['SnowDepth_inst','SWE_inst']

# The LVT invocations for JULES LSM output.
_LVT_JULES_INVOCATIONS_3HR = ['AvgSurfT_inst', 'AvgSurfT_tavg', 
                              'CanopInt_inst',
                              'Elevation_inst', 'Evap_tavg', 
                              'LWdown_f_inst', 'LWdown_f_tavg',
                              'Landcover_inst', 'Landmask_inst',
                              'Psurf_f_inst', 'Psurf_f_tavg', 
                              'Qair_f_inst', 'Qair_f_tavg',
                              'Qh_tavg', 'Qle_tavg', 
                              'Qs_acc','Qsb_acc',
                              'RelSMC_inst','RHMin_inst',
                              'SWE_inst',
                              'SWdown_f_inst', 'SWdown_f_tavg', 
                              'SnowDepth_inst', 
                              'SoilMoist_inst', 'SoilMoist_tavg',
                              'SoilTemp_inst', 'SoilTemp_tavg',
                              'Soiltype_inst', 
                              'Tair_f_inst', 'Tair_f_max', 
                              'Tair_f_tavg',
                              'TotalPrecip_acc', 'Wind_f_inst', 'Wind_f_tavg']

_LVT_JULES_INVOCATIONS_24HR = ['Evap_tavg', 'LWdown_f_tavg',
                               'RHMin_inst',
                               'SoilMoist_tavg', 'SoilTemp_tavg',
                               'SWdown_f_tavg','Tair_f_max',
                               'Tair_f_tavg',
                               'TotalPrecip_acc','Wind_f_tavg']

# The 24-hr postprocessing should include the latest 3-hr snow depth and SWE.
_LVT_JULES_INVOCATIONS_24HR_LATEST = ['SnowDepth_inst','SWE_inst']

# The combined invocation directory for all supported LSMs.
_INVOCATIONS = {
    "NOAH_3HR" : _LVT_NOAH_INVOCATIONS_3HR,
    "NOAH_24HR" : _LVT_NOAH_INVOCATIONS_24HR,
    "NOAH_24HR_LATEST" : _LVT_NOAH_INVOCATIONS_24HR_LATEST,
    "JULES_3HR" : _LVT_JULES_INVOCATIONS_3HR,
    "JULES_24HR" : _LVT_JULES_INVOCATIONS_24HR,
    "JULES_24HR_LATEST" : _LVT_JULES_INVOCATIONS_24HR_LATEST,
}

# The Noah variables handled by each LVT invocation.
_LIS_NOAH_VARIABLES_3HR = {}
for var in _LVT_NOAH_INVOCATIONS_3HR:
    if var == "RHMin_inst":
        _LIS_NOAH_VARIABLES_3HR[var] = [var,"Tair_f_min"]
    else:
        _LIS_NOAH_VARIABLES_3HR[var] = [var]

# The Noah variables handled by each LVT invocation.
_LIS_NOAH_VARIABLES_24HR = {}
for var in _LVT_NOAH_INVOCATIONS_24HR:
    if var == "RHMin_inst":
        _LIS_NOAH_VARIABLES_24HR[var] = [var,"Tair_f_min"]
    else:
        _LIS_NOAH_VARIABLES_24HR[var] = [var]

_LIS_NOAH_VARIABLES_24HR_LATEST = {}
for var in _LVT_NOAH_INVOCATIONS_24HR_LATEST:
    _LIS_NOAH_VARIABLES_24HR_LATEST[var] = [var]

# The JULES variables handled by each LVT invocation.
_LIS_JULES_VARIABLES_3HR = {}
for var in _LVT_JULES_INVOCATIONS_3HR:
    _LIS_JULES_VARIABLES_3HR[var] = [var]

# The JULES variables handled by each LVT invocation.
_LIS_JULES_VARIABLES_24HR = {}
for var in _LVT_JULES_INVOCATIONS_24HR:
    if var == "RHMin_inst":
        _LIS_JULES_VARIABLES_24HR[var] = [var,"Tair_f_min"]
    else:
        _LIS_JULES_VARIABLES_24HR[var] = [var]

_LIS_JULES_VARIABLES_24HR_LATEST = {}
for var in _LVT_JULES_INVOCATIONS_24HR_LATEST:
    _LIS_JULES_VARIABLES_24HR_LATEST[var] = [var]

# Combined breakdown of all variables handled by LSM and LVT invocation.
_LIS_VARIABLES = {
    "NOAH_3HR" : _LIS_NOAH_VARIABLES_3HR,
    "NOAH_24HR" : _LIS_NOAH_VARIABLES_24HR,
    "NOAH_24HR_LATEST" : _LIS_NOAH_VARIABLES_24HR_LATEST,
    "JULES_3HR" : _LIS_JULES_VARIABLES_3HR,
    "JULES_24HR" : _LIS_JULES_VARIABLES_24HR,
    "JULES_24HR_LATEST" : _LIS_JULES_VARIABLES_24HR_LATEST,
}

# These variables are processed by all invocations, and are either generated
# automatically or are read in from US Navy GOFS files.
_OTHER_VARIABLES = ["latitude","longitude","time","water_temp","aice","hi"]

#------------------------------------------------------------------------------
# Print command line usage
def usage():
    print "Usage: %s yyyymmddhh lsm period [--nospread]" %(sys.argv[0])
    print "   where:"
    print "           yyyymmddhh is valid year/month/day/hour in UTC"
    print "           lsm is name of land surface model used by LIS"
    print "           period is processing time length in hours (3 or 24)"
    print "           --nospread is optional flag to skip ensemble spread"

#------------------------------------------------------------------------------
# Read command line arguments
def read_cmd_args():

    # Check if argument count is correct
    if len(sys.argv) not in [4,5]:
        print "[ERR] Invalid number of command line arguments!"
        usage()
        sys.exit(1)

    # Convert yyyymmddhh argument to a datetime object
    yyyymmddhh = sys.argv[1]
    try:
        year   = int(yyyymmddhh[0:4])
        month  = int(yyyymmddhh[4:6])
        day    = int(yyyymmddhh[6:8])
        hour   = int(yyyymmddhh[8:10])
        validdt = datetime.datetime(year,month,day,hour)
    except:
        print "[ERR] Cannot process yyyymmddhh argument!"
        usage()
        sys.exit(1)

    # Get lsm name
    lsm = None
    for i in range(0,len(_LIS_LSMS)):
        if sys.argv[2] == _LIS_LSMS[i]:
            lsm = _LIS_LSMS[i]
            break
    if lsm == None:
        print "[ERR] Invalid lsm selection!"
        print " lsm value is %s" %(sys.argv[2])
        text = " Supported lsms:" 
        for lsm in _LIS_LSMS:
            text += " %s" %(lsm)
        print text
        sys.exit(1)

    # Get period
    period_options = [3,24]
    period = None
    for i in range(0,len(period_options)):
        if int(sys.argv[3]) == period_options[i]:
            period = period_options[i]
            break
    if period == None:
        print "[ERR] Invalid time period selected!"
        print "  Read in %s" %(sys.argv[3])
        print "  Only supports 3 or 24!"
        sys.exit(1)

    # Check if we are skipping ensemble spread
    skip_ens_spread = False
    if len(sys.argv) == 5:
        if sys.argv[4] == "--nospread":
            skip_ens_spread = True
        else:
            print "[ERR] Invalid argument %s" %(sys.argv[4])
            usage()
    
    # See if ncks exists and is executable by current user (the script)
    # This used to be specified on the command line, but is now hardwired
    # to better comply with Air Force security requirements.
    ncks = _NCKS_PATH
    if not os.path.isfile(ncks):
        print "[ERR] Binary %s does not exist!" %(ncks)
        print "[ERR] Modify %s to correct the path to ncks!" %(sys.argv[0])
        sys.exit(1)
    if not os.access(ncks,os.X_OK):
        print "[ERR] %s cannot be executed by current user!" %(ncks)
        print "[ERR] Modify %s to correct the path to ncks!" %(sys.argv[0])
        sys.exit(1)

    return validdt,lsm,ncks,period,skip_ens_spread

#------------------------------------------------------------------------------
# Collect netCDF mean files
def get_nc_mean_files(validdt,lsm,period):

    key = "%s_%dHR" %(lsm,period)
    invocation_list = _INVOCATIONS[key]

    mean_nc_infiles = {}

    # Collect input files
    for invocation in invocation_list:
        # FIXME -- Let user configure output directory prefix
        path = "STATS.%s.%shr" %(invocation,period)

        # FIXME -- Let user configure file name
        path += "/PS.557WW_SC.U_DI.C_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA"
        if period == 24:
            path += ".LIS24_DD."
        else:
            path += ".LIS_DD."
        path += "%4.4d%2.2d%2.2d_DT" %(validdt.year,
                                       validdt.month,
                                       validdt.day)
        path += ".%2.2d00_DF" %(validdt.hour)

        mean_path = path + ".nc"
        if not os.path.exists(mean_path):
            print "[ERR], %s does not exist!" %(mean_path)
            sys.exit(1)

        mean_nc_infiles[invocation] = mean_path

    # Get output file
    # FIXME -- Let user configure output directory prefix
    path = "STATS_merged_%shr" %(period)
    if not os.path.exists(path):
        os.mkdir(path)

    # FIXME -- Let user configure file name
    path += "/PS.557WW_SC.U_DI.C_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA"
    if period == 24:
        path += ".LIS24_DD."
    else:
        path += ".LIS_DD."
    path += "%4.4d%2.2d%2.2d_DT" %(validdt.year,
                                   validdt.month,
                                   validdt.day)
    path += ".%2.2d00_DF" %(validdt.hour)

    mean_nc_outfile = path + ".nc"

    # All done
    return mean_nc_infiles, mean_nc_outfile

#------------------------------------------------------------------------------
# Collect netCDF ssdev files
def get_nc_ssdev_files(validdt,lsm,period):

    key = "%s_%dHR" %(lsm,period)
    invocation_list = _INVOCATIONS[key]

    ssdev_nc_infiles = {}

    # Collect input files
    for invocation in invocation_list:
        # FIXME -- Let user configure output directory prefix
        path = "STATS.%s.%shr" %(invocation,period)

        # FIXME -- Let user configure file name
        path += "/PS.557WW_SC.U_DI.C_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA"
        if period == 24:
            path += ".LIS24_DD."
        else:
            path += ".LIS_DD." 
        path += "%4.4d%2.2d%2.2d_DT" %(validdt.year,
                                       validdt.month,
                                       validdt.day)
        path += ".%2.2d00_DF" %(validdt.hour)

        ssdev_path = path + "_SSDEV.nc"
        if not os.path.exists(ssdev_path):
            print "[ERR], %s does not exist!" %(ssdev_path)
            sys.exit(1)

        ssdev_nc_infiles[invocation] = ssdev_path

    # Get output file
    # FIXME -- Let user configure output directory prefix
    path = "STATS_merged_%shr" %(period)
    if not os.path.exists(path):
        os.mkdir(path)

    # FIXME -- Let user configure file name
    path += "/PS.557WW_SC.U_DI.C_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA"
    if period == 24:
        path += ".LIS24_DD."
    else:
        path += ".LIS_DD."
    path += "%4.4d%2.2d%2.2d_DT" %(validdt.year,
                                   validdt.month,
                                   validdt.day)
    path += ".%2.2d00_DF" %(validdt.hour)

    ssdev_nc_outfile = path + "_SSDEV.nc"

    # All done
    return ssdev_nc_infiles, ssdev_nc_outfile

#------------------------------------------------------------------------------
# Collect netCDF latest files
def get_nc_latest_files(validdt,lsm):

    key = "%s_24HR_LATEST" %(lsm)
    invocation_list = _INVOCATIONS[key]

    latest_nc_infiles = {}

    # Collect input files
    for invocation in invocation_list:
        # FIXME -- Let user configure output directory prefix
        path = "STATS.%s.3hr" %(invocation) # Always use 3hr processing

        # FIXME -- Let user configure file name
        path += "/PS.557WW_SC.U_DI.C_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA"
        path += ".LIS_DD." # Always use 3-hr output for latest fields
        path += "%4.4d%2.2d%2.2d_DT" %(validdt.year,
                                       validdt.month,
                                       validdt.day)
        path += ".%2.2d00_DF" %(validdt.hour)

        latest_path = path + ".nc"
        if not os.path.exists(latest_path):
            print "[ERR], %s does not exist!" %(latest_path)
            sys.exit(1)

        latest_nc_infiles[invocation] = latest_path

    # All done
    return latest_nc_infiles

#------------------------------------------------------------------------------
# Use ncks to merge netCDF fields together
def merge_nc_files(lsm,ncks,period,nc_infiles, \
                       nc_outfile, latest_nc_infiles=None):

    key = "%s_%dHR" %(lsm,period)

    # Start with ensemble mean
    cmd = "cp %s %s" %(nc_infiles[_INVOCATIONS[key][0]], \
                           nc_outfile)
    print cmd
    rc = subprocess.call(cmd,shell=True)
    if rc != 0:
        print "[ERR] Problem with cp!"
        sys.exit(1)

    invocations = _INVOCATIONS[key][1:]
    for invocation in invocations:
        variables = _LIS_VARIABLES[key][invocation]
        for variable in variables:
            cmd = "%s -A -v %s %s %s" \
                %(ncks,variable,nc_infiles[invocation],nc_outfile)
            print cmd
            rc = subprocess.call(cmd,shell=True)
            if rc != 0:
                print "[ERR] Problem with ncks!"
                sys.exit(1)

    # For 24-hr postprocessing, we also must concatenate several 3-hr fields
    if latest_nc_infiles != None:
        key = "%s_24HR_LATEST" %(lsm)
        invocations = _INVOCATIONS[key][:]
        for invocation in invocations:
            variables = _LIS_VARIABLES[key][invocation]
            for variable in variables:
                cmd = "%s -A -v %s %s %s" \
                    %(ncks,variable,
                      latest_nc_infiles[invocation],nc_outfile)
                print cmd
                rc = subprocess.call(cmd,shell=True)
                if rc != 0:
                    print "[ERR] Problem with ncks!"
                    sys.exit(1)

#------------------------------------------------------------------------------
# Main Driver.

if __name__ == "__main__":

    # Process command line arguments
    validdt,lsm,ncks,period,skip_ens_spread = read_cmd_args()
    
    # Collect netCDF files
    (mean_nc_infiles, mean_nc_outfile) = \
        get_nc_mean_files(validdt,lsm,period)
    # 3-hr postprocessing includes ensemble spread files
    if period == 3 and not skip_ens_spread:
        (ssdev_nc_infiles, ssdev_nc_outfile) = \
            get_nc_ssdev_files(validdt,lsm,period)
    # 24-hr postprocessing includes several latest 3-hr fields
    if period == 24:
        latest_nc_infiles = get_nc_latest_files(validdt,lsm)

    # Merge the input netCDF files together
    if period == 3:
        merge_nc_files(lsm,ncks,period,mean_nc_infiles,mean_nc_outfile)
        if not skip_ens_spread:
            merge_nc_files(lsm,ncks,period,ssdev_nc_infiles,ssdev_nc_outfile)
    else:
        # 24-hr processing
        merge_nc_files(lsm,ncks,period,mean_nc_infiles,
                       mean_nc_outfile,latest_nc_infiles)


    
    
