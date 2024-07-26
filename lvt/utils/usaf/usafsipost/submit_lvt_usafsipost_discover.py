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
# SCRIPT: submit_lvt_usafsipost_discover.py
#
# PURPOSE: Constructs a batch script for running LVT in USAFSIpost runmode,
# and submits batch job to queue on NASA Discover supercomputer.
#
# REQUIREMENTS as of 15 Jul 2021:
# * Python 3.8
#
# REVISION HISTORY:
# 15 Jul 2021: Eric Kemp (SSAI), first version.
# 19 Jan 2022: Eric Kemp (SSAI), Discover updates.
# 08 Dec 2022: Eric Kemp (SSAI), refactored to increase pylint score.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import os
import subprocess
import sys

def _usage():
    """Prints usage statement for script."""
    print(f"Usage: {sys.argv[0]} chargecode qos")
    print("  where:")
    print("    chargecode is SLURM account")
    print("    qos is the SLURM quality-of-service")

def _main():
    """Main driver"""

    # Check command-line arguments
    if len(sys.argv) != 3:
        print("[ERR], problem with command line arguments!")
        _usage()
        sys.exit(1)

    account = sys.argv[1]
    qos = sys.argv[2]

    # Make sure LVT executable is in place before launching job
    if not os.path.exists("LVT"):
        print("[ERR], LVT executable does not exist!")
        sys.exit(1)

    # Create a batch script.
    scriptname = "run_lvt.usafsipost.sh"
    with open(scriptname, "w", encoding="ascii") as file:
        line = f"""#!/bin/sh
#SBATCH --account {account}
#SBATCH --constraint="hasw|sky|cas"
#SBATCH --job-name=usafsipost
#SBATCH --ntasks=1
#SBATCH --output usafsipost.slurm.out
#SBATCH --qos={qos}
#SBATCH --time=0:05:00

if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# NOTE: This privatemodule can be found in LISF/env/discover
module purge
module use --append ~/privatemodules
module load lisf_7.5_intel_2021.4.0

if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

if [ ! -e lvt.config.usafsipost ] ; then
   echo "ERROR, lvt.config.usafsipost does not exist!" && exit 1
fi

mpirun -np 1 ./LVT lvt.config.usafsipost || exit 1

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
