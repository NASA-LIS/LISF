#!/usr/bin/env python3
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
#
#------------------------------------------------------------------------------
"""

# Standard modules
import os
import subprocess
import sys

def _usage():
    """Prints usage statement for script."""
    print("Usage: %s chargecode qos" %(sys.argv[0]))
    print("  where:")
    print("    chargecode is SLURM account")
    print("    qos is the SLURM quality-of-service")

# Main driver
if __name__ == "__main__":

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
    SCRIPTNAME = "run_lvt.usafsipost.sh"
    f = open(SCRIPTNAME, "w")
    line = """#!/bin/sh
#SBATCH --account %s
#SBATCH --constraint="hasw|sky|cas"
#SBATCH --job-name=usafsipost
#SBATCH --ntasks=1
#SBATCH --output usafsipost.slurm.out
#SBATCH --qos=%s
#SBATCH --time=0:05:00

if [ ! -z $SLURM_SUBMIT_DIR ] ; then
    cd $SLURM_SUBMIT_DIR || exit 1
fi

# NOTE: This privatemodule can be found in LISF/env/discover
module purge
module use --append ~/privatemodules
module load lisf_7_intel_2021.4.0_s2s

if [ ! -e ./LVT ] ; then
   echo "ERROR, LVT does not exist!" && exit 1
fi

if [ ! -e lvt.config.usafsipost ] ; then
   echo "ERROR, lvt.config.usafsipost does not exist!" && exit 1
fi

mpirun -np 1 ./LVT lvt.config.usafsipost || exit 1

exit 0

""" %(account, qos)
    f.write(line)
    f.close()

    # Submit the batch job to SLURM
    # NOTE: pylint is insistent on treating CMD and RC as constants, and
    # thus requiring UPPER_CASE naming style.
    CMD = "sbatch %s" %(SCRIPTNAME)
    print(CMD)
    RC = subprocess.call(CMD, shell=True)
    if RC != 0:
        print("[ERR] Problem with sbatch!")
        sys.exit(1)
