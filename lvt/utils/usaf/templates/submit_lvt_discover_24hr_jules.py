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
Sample script for submitting LVT postprocessing batch jobs on Discover for
jules50 for 557WW.
"""

import os
import subprocess
import sys
import time

_VARS = ["SoilMoist_tavg", "SoilTemp_tavg",
        "RHMin_inst",
        "Evap_tavg", "LWdown_f_tavg",
        "SWdown_f_tavg",
        "Tair_f_max",
        "Tair_f_tavg",
        "TotalPrecip_acc", "Wind_f_tavg"]

def _main():
    """Main driver"""

    if not os.path.exists("LVT"):
        print("ERROR, LVT executable does not exist!")
        sys.exit(1)

    for var in _VARS:
        scriptname = f"run_lvt.{var}_24hr.sh"
        with open(scriptname, "w", encoding="ascii") as file:
            line = f"""#!/bin/sh
#SBATCH --job-name={var}.24hr
#SBATCH --time=1:00:00
#SBATCH --account s1189
#SBATCH --output {var}.24hr.slurm.out
#Adjust node, core, and hardware constraints here
#SBATCH --ntasks=1 --constraint="[cas|sky]"

if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

module purge
module use --append ~/privatemodules
#module load lisf_7.5_intel_2021.4.0
module load lisf_7.5_intel_2021.4.0_s2s

if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

if [ ! -e configs/lvt.config.{var}.24hr ] ; then
   echo "ERROR, configs/lvt.config.{var}.24hr does not exist!" && exit 1
fi

mkdir -p OUTPUT/STATS.{var}.24hr || exit 1

time mpirun -np 1 ./LVT configs/lvt.config.{var}.24hr || exit 1

exit 0
"""
            file.write(line)

        cmd = f"sbatch {scriptname}"
        print(cmd)
        err = subprocess.call(cmd, shell=True)
        if err != 0:
            print("[ERR] Problem with sbatch!")
            sys.exit(1)
        time.sleep(1)  # Don't overwhelm SLURM

if __name__ == "__main__":
    _main()
