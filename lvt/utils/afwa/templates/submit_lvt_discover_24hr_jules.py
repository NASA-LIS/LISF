#!/usr/bin/env python

import os
import subprocess
import sys
import time

vars = ["SoilMoist_tavg", "SoilTemp_tavg",
        "RHMin_inst",
        "Evap_tavg", "LWdown_f_tavg", 
        "SWdown_f_tavg",
        "Tair_f_max",
        "Tair_f_tavg",
        "TotalPrecip_acc", "Wind_f_tavg"]

if not os.path.exists("LVT"):
    print "ERROR, LVT executable does not exist!"
    sys.exit(1)

for var in vars:
    scriptname = "run_lvt.%s_24hr.sh" %(var)
    f = open(scriptname,"w")
    line = """#!/bin/sh
#SBATCH --job-name=%s.24hr
#SBATCH --time=0:15:00
#SBATCH --account s0942
#SBATCH --output %s.24hr.slurm.out
#Adjust node, core, and hardware constraints here
#SBATCH --ntasks=1 --constraint=hasw
#Set quality of service, if needed.
#SBATCH --qos=high

if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

#source /discover/nobackup/emkemp/AFWA/modules.sh
module purge
module load ~jvgeiger/privatemodules/lis_7_intel_18_0_3_222


if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

if [ ! -e lvt.config.%s.24hr ] ; then
   echo "ERROR, lvt.config.%s.24hr does not exist!" && exit 1
fi

mkdir -p STATS.%s.24hr || exit 1

time ./LVT lvt.config.%s.24hr || exit 1

exit 0
""" %(var,var,var,var,var,var)
    f.write(line)
    f.close()

    cmd = "sbatch %s" %(scriptname)
    print cmd
    rc = subprocess.call(cmd,shell=True)
    if rc != 0:
        print "[ERR] Problem with sbatch!"
        sys.exit(1)
    time.sleep(1) # Don't overwhelm SLURM


