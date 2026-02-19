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
#--------------------------------------------------------------------------
#
# SCRIPT: submit_lvt_snippost_hpc11.py
#
# PURPOSE: Constructs a batch script for running LVT in SNIPpost runmode,
# and submits batch job to queue on HPC11.
#
# REQUIREMENTS as of 04 Aug 2025:
# * Python 3.11
#
# REVISION HISTORY:
# 28 Oct 2024: Eric Kemp (SSAI), first version.
# 04 Aug 2025: Eric Kemp (SSAI), SNIP version.
#
#--------------------------------------------------------------------------
"""

# Standard modules
import os
import subprocess
import sys

def _usage():
    """Prints usage statement for script."""
    print(f"Usage: {sys.argv[0]} chargecode")
    print("  where:")
    print("    chargecode is SLURM account")

def _main():
    """Main driver"""

    # Check command-line arguments
    if len(sys.argv) != 2:
        print("[ERR], problem with command line arguments!")
        _usage()
        sys.exit(1)

    account = sys.argv[1]

    # Make sure LVT executable is in place before launching job
    if not os.path.exists("LVT"):
        print("[ERR], LVT executable does not exist!")
        sys.exit(1)

    # Create a batch script.
    scriptname = "run_lvt.snippost.sh"
    with open(scriptname, "w", encoding="ascii") as file:
        line = f"""#!/bin/sh
#SBATCH --job-name=snippost
#SBATCH --time=0:05:00
#SBATCH --account {account}
#SBATCH --output snippost.slurm.out
#SBATCH --ntasks=1
#SBATCH --cluster-constraint=blue
#SBATCH --exclusive
#SBATCH --mem=0

if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# Environment
module use --append /ccs/home/emkemp/hpc11/privatemodules
module load lisf_7.6_prgenv_cray_8.5.0_cpe_23.12
module load afw-python/3.11-202406

if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

lvtconfig=lvt.config.foc.snippost.76

if [ ! -e $lvtconfig ] ; then
   echo "ERROR, $lvtconfig does not exist!" && exit 1
fi

mpirun -np 1 ./LVT $lvtconfig || exit 1

exit 0

"""
        file.write(line)

    # Submit the batch job to SLURM
    cmd = f"sbatch {scriptname}"
    print(cmd)
    err = subprocess.call(cmd, shell=True)
    if err != 0:
        print("[ERR] Problem with sbatch!")
        sys.exit(1)

# Main driver
if __name__ == "__main__":
    _main()
