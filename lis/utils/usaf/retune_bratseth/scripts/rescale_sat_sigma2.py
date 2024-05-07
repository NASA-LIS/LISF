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
SCRIPT: rescale_sat_sigma2.py

Rescales estimate of satellite precip error variance based on gage comparison
for use with NWP.

REVISION HISTORY:
03 Nov 2020:  Eric Kemp.  Initial specification.
"""

import os
import sys

def usage():
    """Print usage message for this script."""
    print(f"Usage: {sys.argv[0]} SATDATATYPE")
    print("   SATDATATYPE is source of satellite precip estimates")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        usage()
        sys.exit(1)

    # What satellite data source are we using?
    SATDATA = sys.argv[1]
    if SATDATA not in ["cmorph", "imerg", "geoprecip", "ssmi"]:
        print(f"[ERR] Invalid satellite data type {SATDATA}")
        sys.exit(1)

    # Pull the gage error variance estimated by comparison with NWP.
    FILE1 = "gage_nwp.param"
    if not os.path.exists(FILE1):
        print(f"[ERR] {FILE1} does not exist!")
        sys.exit(1)
    with open(FILE1, "r", encoding="ascii") as file:
        lines = file.readlines()
    SIGMA2_GAGE_NWP = -9999
    for line in lines:
        if "SIGMA2_obs:" in line:
            SIGMA2_GAGE_NWP = float(line.split()[-1])
            break
    if SIGMA2_GAGE_NWP < 0:
        print("[ERR] No valid value found for SIGMA2_GAGE_NWP!")
        sys.exit(1)

    # Pull the error values estimated by gage vs satellite comparison
    FILE2 = f"gage_{SATDATA}.param"
    if not os.path.exists(FILE2):
        print(f"[ERR] {FILE2} does not exist!")
        sys.exit(1)
    with open(FILE2, "r", encoding="ascii") as file:
        lines = file.readlines()
    SIGMA2_GAGE_SAT = -9999
    SIGMA2_SAT = -9999
    L_SAT = -9999
    for line in lines:
        if "SIGMA2_obs:" in line:
            SIGMA2_GAGE_SAT = float(line.split()[-1])
            continue
        if "SIGMA2_back:" in line:
            SIGMA2_SAT = float(line.split()[-1])
            continue
        if "L_back:" in line:
            L_SAT = float(line.split()[-1])
            continue
    if SIGMA2_GAGE_SAT < 0:
        print("[ERR] No valid value found for SIGMA2_GAGE_SAT!")
        sys.exit(1)
    if SIGMA2_SAT < 0:
        print("[ERR] No valid value found for SIGMA2_SAT!")
        sys.exit(1)
    if L_SAT < 0:
        print("[ERR] No valid value found for L_SAT!")
        sys.exit(1)

    # Rescale the satellite error variance
    SIGMA2_SAT = (SIGMA2_GAGE_NWP / SIGMA2_GAGE_SAT) * SIGMA2_SAT

    # Now write the new param file
    with open(f"gage_{SATDATA}_rescaled.param", "w", encoding="ascii") \
         as outfile:
        outfile.write(f"SIGMA2_obs: {SIGMA2_GAGE_NWP}\n")
        outfile.write(f"SIGMA2_back: {SIGMA2_SAT}\n")
        outfile.write(f"L_back: {L_SAT}\n")
