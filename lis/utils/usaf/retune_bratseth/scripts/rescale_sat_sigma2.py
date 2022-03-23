#!/usr/bin/env python3
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
    print("Usage: %s SATDATATYPE" %(sys.argv[0]))
    print("   SATDATATYPE is source of satellite precip estimates")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        usage()
        sys.exit(1)

    # What satellite data source are we using?
    satdata = sys.argv[1]
    if satdata not in ["cmorph", "imerg", "geoprecip", "ssmi"]:
        print("[ERR] Invalid satellite data type %s" %(satdata))
        sys.exit(1)

    # Pull the gage error variance estimated by comparison with NWP.
    file1 = "gage_nwp.param"
    if not os.path.exists(file1):
        print("[ERR] %s does not exist!" %(file1))
        sys.exit(1)
    lines = open(file1, "r").readlines()
    sigma2_gage_nwp = -9999
    for line in lines:
        if "SIGMA2_obs:" in line:
            sigma2_gage_nwp = float(line.split()[-1])
            break
    if sigma2_gage_nwp < 0:
        print("[ERR] No valid value found for sigma2_gage_nwp!")
        sys.exit(1)

    # Pull the error values estimated by gage vs satellite comparison
    file2 = "gage_%s.param" %(satdata)
    if not os.path.exists(file2):
        print("[ERR] %s does not exist!" %(file2))
        sys.exit(1)
    lines = open(file2, "r").readlines()
    sigma2_gage_sat = -9999
    sigma2_sat = -9999
    L_sat = -9999
    for line in lines:
        if "SIGMA2_obs:" in line:
            sigma2_gage_sat = float(line.split()[-1])
            continue
        if "SIGMA2_back:" in line:
            sigma2_sat = float(line.split()[-1])
            continue
        if "L_back:" in line:
            L_sat = float(line.split()[-1])
            continue
    if sigma2_gage_sat < 0:
        print("[ERR] No valid value found for sigma2_gage_sat!")
        sys.exit(1)
    if sigma2_sat < 0:
        print("[ERR] No valid value found for sigma2_sat!")
        sys.exit(1)
    if L_sat < 0:
        print("[ERR] No valid value found for L_sat!")
        sys.exit(1)

    # Rescale the satellite error variance
    sigma2_sat = (sigma2_gage_nwp / sigma2_gage_sat) * sigma2_sat

    # Now write the new param file
    fd = open("gage_%s_rescaled.param" %(satdata), "w")
    fd.write("SIGMA2_obs: %s\n" %(sigma2_gage_nwp))
    fd.write("SIGMA2_back: %s\n" %(sigma2_sat))
    fd.write("L_back: %s\n" %(L_sat))
    fd.close()
