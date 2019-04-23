#!/usr/bin/env python
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
#
#------------------------------------------------------------------------------

# Standard modules
import datetime
import os
import subprocess
import sys

#------------------------------------------------------------------------------

# Supported LIS LSMs
_LIS_LSMS = ["NOAH","JULES"]

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
                              'SWdown_f_tavg','Tair_f_max',
                              'Tair_f_tavg',
                              'TotalPrecip_acc','Wind_f_tavg']

# The 24-hr postprocessing should include the latest 3-hr snow depth and SWE.
_LVT_NOAH_INVOCATIONS_24HR_LATEST = ['SnowDepth_inst','SWE_inst']

# The LVT invocations for JULES LSM output.
_LVT_JULES_INVOCATIONS_3HR = ['AvgSurfT_inst', 'AvgSurfT_tavg', 
                              'CanopInt_inst',
                              'Elevation_inst', 'Evap_tavg', 
                              'Greenness_inst',
                              'LWdown_f_inst', 'LWdown_f_tavg',
                              'Landcover_inst', 'Landmask_inst',
                              'Psurf_f_inst', 'Psurf_f_tavg', 
                              'Qair_f_inst', 'Qair_f_tavg',
                              'Qh_tavg', 'Qle_tavg', 
                              'Qs_acc','Qsb_acc',
                              'RHMin_inst','RelSMC_inst',
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

#------------------------------------------------------------------------------
# Print command line usage
def usage():
    print "Usage: %s yyyymmddhh lsm period [--nospread]" %(sys.argv[0])
    print "   where:"
    print "         yyyymmddhh is valid year/month/day/hour in UTC"
    print "         lsm is name of land surface model used by LIS"
    print "         period is time period (hours) for postprocessing (3 or 24)"
    print "         --nospread is optional flag to skip ensemble spread"

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
      
    # Get processing hour
    period_options = [3,24]
    period = None
    tmp_int = sys.argv[3]
    for i in range(0,len(period_options)):
        if int(sys.argv[3]) == period_options[i]:
            period = period_options[i]
            break
    if period == None:    
        print "[ERR] Invalid period selection!"
        print " period value is %s" %(sys.argv[3])
        print " Supported time periods are: 3 and 24"
        sys.exit(1)

    # Check if ensemble spread should be skipped
    skip_ens_spread = False
    if len(sys.argv) == 5:
        if sys.argv[4] == "--nospread":
            skip_ens_spread = True
        else:
            print "[ERR] Invalid argument %s" %(sys.argv[4])
            usage()

    return validdt,lsm,period,skip_ens_spread

#------------------------------------------------------------------------------
# Collect GRIB2 mean files
def get_gr2_mean_files(validdt,lsm,period):

    key = "%s_%dHR" %(lsm,period)
    invocation_list = _INVOCATIONS[key]

    mean_gr2_infiles = {}
    
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

        mean_path = path + ".GR2"
        if not os.path.exists(mean_path):
            print "[ERR], %s does not exist!" %(mean_path)
            sys.exit(1)
        mean_gr2_infiles[invocation] = mean_path
        
    # Get output files
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

    mean_gr2_outfile = path + ".GR2"

    # All done
    return mean_gr2_infiles, mean_gr2_outfile

#------------------------------------------------------------------------------
# Collect GRIB2 ssdev files
def get_gr2_ssdev_files(validdt,lsm,period):

    key = "%s_%dHR" %(lsm,period)
    invocation_list = _INVOCATIONS[key]

    ssdev_gr2_infiles = {}
    
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

        ssdev_path = path + "_SSDEV.GR2"
        if not os.path.exists(ssdev_path):
            print "[ERR], %s does not exist!" %(ssdev_path)
            sys.exit(1)
        ssdev_gr2_infiles[invocation] = ssdev_path

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

    ssdev_gr2_outfile = path + "_SSDEV.GR2"

    # All done
    return ssdev_gr2_infiles, ssdev_gr2_outfile

#------------------------------------------------------------------------------
# Collect GRIB2 latest files
def get_gr2_latest_files(validdt,lsm):

    key = "%s_24HR_LATEST" %(lsm) 
    invocation_list = _INVOCATIONS[key]

    latest_gr2_infiles = {}
    
    # Collect input files
    for invocation in invocation_list:
        # FIXME -- Let user configure output directory prefix
        path = "STATS.%s.3hr" %(invocation) # Always use 3hr processing

        # FIXME -- Let user configure file name
        path += "/PS.557WW_SC.U_DI.C_GP.LIS_GR.C0P09DEG_AR.GLOBAL_PA"
        path += ".LIS_DD."
        path += "%4.4d%2.2d%2.2d_DT" %(validdt.year,
                                       validdt.month,
                                       validdt.day)
        path += ".%2.2d00_DF" %(validdt.hour)

        latest_path = path + ".GR2"
        if not os.path.exists(latest_path):
            print "[ERR], %s does not exist!" %(latest_path)
            sys.exit(1)
        latest_gr2_infiles[invocation] = latest_path
        
    # All done
    return latest_gr2_infiles

#------------------------------------------------------------------------------
# Use cat to merge GRIB2 fields together
def merge_gr2_files(lsm, period,gr2_infiles, gr2_outfile, 
                    latest_gr2_infiles=None):

    key = "%s_%dHR" %(lsm,period)
    invocations = _INVOCATIONS[key][0:]    
    cmd = "cat"
    for invocation in invocations:
        cmd += " %s" %(gr2_infiles[invocation])
    # For 24-hr postprocessing, we also must concatenate several 3-hr fields
    if latest_gr2_infiles != None:
        key = "%s_24HR_LATEST" %(lsm)
        invocations = _INVOCATIONS[key][:]
        for invocation in invocations:
            cmd += " %s" %(latest_gr2_infiles[invocation])
    cmd += " > %s" %(gr2_outfile)

    print cmd
    rc = subprocess.call(cmd,shell=True)
    if rc != 0:
        print "[ERR] Problem with cat!"
        sys.exit(1)

#------------------------------------------------------------------------------
# Main Driver.

if __name__ == "__main__":

    # Process command line arguments
    validdt,lsm,period,skip_ens_spread = read_cmd_args()

    # Collect GRIB2 files
    (mean_gr2_infiles, mean_gr2_outfile) = \
        get_gr2_mean_files(validdt,lsm,period)
    # 3-hr postprocessing includes ensemble spread files
    if period == 3 and not skip_ens_spread:
        (ssdev_gr2_infiles, ssdev_gr2_outfile) = \
            get_gr2_ssdev_files(validdt,lsm,period)
    # 24-hr postprocessing includes several latest 3-hr fields
    if period == 24:
        latest_gr2_infiles = get_gr2_latest_files(validdt,lsm)
        
    # Merge the input GRIB2 files together
    if period == 3:
        merge_gr2_files(lsm,period,mean_gr2_infiles,mean_gr2_outfile) 
        if not skip_ens_spread:
            merge_gr2_files(lsm,period,ssdev_gr2_infiles,ssdev_gr2_outfile)
    else:
        # 24-hr processing
        merge_gr2_files(lsm,period,mean_gr2_infiles,
                        mean_gr2_outfile,latest_gr2_infiles)        

