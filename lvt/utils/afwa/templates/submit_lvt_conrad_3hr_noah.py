#!/usr/bin/env python

import os
import subprocess
import sys
import time

vars = ['RelSMC_inst','SmLiqFrac_inst',
        'SoilMoist_inst', 'SoilMoist_tavg',
        'SoilTemp_inst', 'SoilTemp_tavg',
        'RHMin_inst',
        'Albedo_tavg', 'AvgSurfT_inst', 'AvgSurfT_tavg', 
        'CanopInt_inst', 'Elevation_inst', 'Evap_tavg', 
        'Greenness_inst', 
        'LWdown_f_inst', 'LWdown_f_tavg',
        'Landcover_inst', 'Landmask_inst', 'PotEvap_tavg',
        'Psurf_f_inst', 'Psurf_f_tavg', 
        'Qair_f_inst', 'Qair_f_tavg',
        'Qg_tavg', 'Qh_tavg', 'Qle_tavg', 'Qs_acc', 
        'Qsb_acc', 'SWE_inst',
        'SWdown_f_inst', 'SWdown_f_tavg', 
        'SnowDepth_inst', 'Snowcover_inst',
        'Soiltype_inst', 
        'Tair_f_inst', 'Tair_f_max', 
        'Tair_f_tavg',
        'TotalPrecip_acc', 'Wind_f_inst', 'Wind_f_tavg']

# Handle command line
def usage():
    print "Usage:  %s chargecode queue" %(sys.argv[0])
    print "  where chargecode is PBS project_code"
    print "  and   queue is PBS queue OR reservation number"
if len(sys.argv) != 3:
    print "ERROR, problem with command line arguments!"
    usage()
    sys.exit(1)
project_code = sys.argv[1]
reservation = sys.argv[2]

# Make sure LVT executable is in place before launching jobs
if not os.path.exists("LVT"):
    print "ERROR, LVT executable does not exist!"
    sys.exit(1)

# Loop through each invocation, create a batch script, and launch the
# batch script.
for var in vars:
    scriptname = "run_lvt.%s_3hr.sh" %(var)
    f = open(scriptname,"w")
    line = """#!/bin/sh
#PBS -A %s\n""" %(project_code)
    line += """#PBS -j oe
#PBS -l walltime=0:15:00
#PBS -l select=1:ncpus=32
#PBS -N %s.3hr\n""" %(var)
    line += """#PBS -q %s\n""" %(reservation)
    line +="""#PBS -W sandbox=PRIVATE
#PBS -V

module use --append ~jim/README
module load lis_7_intel_17_0_2_174
ulimit -c unlimited
ulimit -m unlimited
ulimit -s unlimited

cd "$PBS_O_WORKDIR" || exit 1
echo `pwd`

if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

if [ ! -e lvt.config.%s.3hr ] ; then
   echo "ERROR, lvt.config.%s.3hr does not exist!" && exit 1
fi

aprun -n 1 -j 1 ./LVT lvt.config.%s.3hr || exit 1

exit 0
""" %(var,var,var)
    f.write(line)
    f.close()

    cmd = "qsub %s" %(scriptname)
    print cmd
    rc = subprocess.call(cmd,shell=True)
    if rc != 0:
        print "[ERR] Problem with qsub!"
        sys.exit(1)
    time.sleep(1) # Don't overwhelm PBS!


