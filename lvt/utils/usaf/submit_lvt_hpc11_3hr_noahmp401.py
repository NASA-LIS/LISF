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
Sample script for submitting LVT postprocessing batch jobs on Discover for
noahmp401 for 557WW.
"""

import os
import subprocess
import sys
import time

_VARS = ['RelSMC_inst', 'SmLiqFrac_inst',
        'SoilMoist_inst', 'SoilMoist_tavg',
        'SoilTemp_inst', 'SoilTemp_tavg',
        'RHMin_inst',
        'Albedo_tavg', 'AvgSurfT_inst', 'AvgSurfT_tavg',
        'CanopInt_inst', 'Elevation_inst', 'Evap_tavg',
        'Greenness_inst',
        'LWdown_f_inst', 'LWdown_f_tavg',
        'Landcover_inst', 'Landmask_inst',
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

def _main():
    """Main driver"""

    if not os.path.exists("LVT"):
        print("ERROR, LVT executable does not exist!")
        sys.exit(1)

    for var in _VARS:
        scriptname = f"run_lvt.{var}_3hr.sh"
        with open(scriptname, "w", encoding="ascii") as file:
            line = f"""#!/bin/sh
#SBATCH --job-name={var}.3hr
#SBATCH --time=1:00:00
#SBATCH --account nwp601
#SBATCH --output {var}.3hr.slurm.out
#Set quality of service, if needed.
#SBATCH --ntasks=1

if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

#module purge
module use --append ~/hpc11/privatemodules
module load lisf_74_prgenv_cray_8.2.0


if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

if [ ! -e configs/lvt.config.{var}.3hr ] ; then
   echo "ERROR, configs/lvt.config.{var}.3hr does not exist!" && exit 1
fi
srun -n 1 ./LVT configs/lvt.config.{var}.3hr || exit 1

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
