#!/usr/bin/env python3
"""
SCRIPT: create_blacklist.py

Script for creating a blacklist.  Will flag stations with unacceptably high or
low mean innovations (observation-minus-background values).  Will also flag
stations with too few reports, and stations that have more than one reported
network.  Satellite observations are ignored.

REVISION HISTORY:
26 Oct 2020: Eric Kemp. Initial Specification.
02 Nov 2020: Eric Kemp. Added datetime constraints.
16 Dec 2020: Eric Kemp. Reorganized most driver code into new function that
             can be called directly (via Python module import) instead of
             via shell call.
"""

# Standard modules
import configparser
import datetime
import glob
import os
import sys

#------------------------------------------------------------------------------
def usage():
    """Print usage statement to standard out"""
    print("Usage: %s CFGFILENAME BLACKLISTFILENAME ENDDATETIME DAYRANGE "
          %(sys.argv[0]))

#------------------------------------------------------------------------------
def is_sat(my_platform):
    """Check if platform is satellite-based"""
    test = False
    if my_platform in ["SSMI", "GEOPRECIP", "CMORPH", "IMERG"]:
        test = True
    return test

#------------------------------------------------------------------------------
def create_blacklist(cfgfile, blacklistfilename, yyyymmddhh, dayrange):
    """Main driver for creating a blacklist file."""

    # Read config file
    if not os.path.exists(cfgfile):
        print("[ERR] Config file %s does not exist!" %(cfgfile))
        sys.exit(1)
    config = configparser.ConfigParser()
    config.read(cfgfile)

    # Process the threshold settings
    network_thresh = config.get('Input', 'network_thresh')
    count_thresh = config.get('Input', 'count_thresh')
    omb_thresh = config.get('Input', 'omb_thresh')

    # Process the filename information.
    data_directory = config.get('Input', 'data_directory')
    if not os.path.exists(data_directory):
        print("[ERR] Directory %s does not exist!" %(data_directory))
    data_frequency = config.get('Input', 'data_frequency')
    data_hours = config.get('Input', 'data_hours')
    data_hour_list = data_hours.split(",")

    data = {}

    # Set start and end datetimes
    year = int(yyyymmddhh[0:4])
    month = int(yyyymmddhh[4:6])
    day = int(yyyymmddhh[6:8])
    hour = int(yyyymmddhh[8:10])
    enddt = datetime.datetime(year=year,
                              month=month,
                              day=day,
                              hour=hour)
    delta = datetime.timedelta(days=int(dayrange))
    startdt = enddt - delta

    # Open the new blacklist file
    fd = open(blacklistfilename, "w")

    # Build list of OBA files
    files = glob.glob("%s/oba_*_%s.txt" %(data_directory,
                                          data_frequency))
    files.sort()

    # Loop through each file, and store the O and B for each platform
    for file in files:

        # Check the valid hour to see if we want to use it.
        chunk = file.split("/")[1][12:14]
        if chunk not in data_hour_list:
            continue

        # Make sure the valid datetime is in the desired range
        yyyymmddhh = file.split("/")[1][4:14]
        year = int(yyyymmddhh[0:4])
        month = int(yyyymmddhh[4:6])
        day = int(yyyymmddhh[6:8])
        hour = int(yyyymmddhh[8:10])
        curdt = datetime.datetime(year=year,
                                  month=month,
                                  day=day,
                                  hour=hour)
        if curdt <= startdt:
            continue
        if curdt >= enddt:
            continue

        # At this point, we trust the file.  Read it.
        lines = open(file, "r").readlines()
        for line in lines[1:]:
            network = line.split()[0]
            platform = line.split()[1]
            if is_sat(platform):
                continue
            if platform not in data:
                data[platform] = []
            ob = line.split()[4]
            back = line.split()[5]
            data[platform].append("%s:%s:%s" %(network, ob, back))

    # Now, loop through each station and calculate the mean OMB.  Create
    # a blacklist
    fd.write("# Rejecting stations with more than %s network\n"
             %(network_thresh))
    fd.write("# Rejecting stations with less than %s observations\n"
             %(count_thresh))
    fd.write("# Rejecting stations with absolute mean OMB beyond %s\n"
             %(omb_thresh))

    platforms = list(data.keys())
    platforms.sort()
    for platform in platforms:
        n = {}
        OMB = {}
        entries = data[platform]
        for entry in entries:
            network = entry.split(":")[0]
            ob = float(entry.split(":")[1])
            back = float(entry.split(":")[2])
            if network not in n:
                n[network] = 1
            else:
                n[network] += 1
            if network not in OMB:
                OMB[network] = 0.
            OMB[network] += (ob - back)
        length = len(OMB.keys())
        if length > int(network_thresh):
            fd.write("%s # Found in %s networks\n" %(platform, length))
            continue
        for network in OMB:
            if n[network] < int(count_thresh):
                fd.write("%s # Only have %s observation(s)\n"
                         %(platform, n[network]))
                continue
            OMB[network] = OMB[network] / float(n[network])
            if abs(OMB[network]) > float(omb_thresh):
                fd.write("%s # Mean OMB is %s\n" %(platform, OMB[network]))
                continue
    fd.close()

#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Check command line
    if len(sys.argv) != 5:
        print("[ERR] Bad command line arguments!")
        usage()
        sys.exit(1)

    # Get command line args
    cfgfile = sys.argv[1]
    blacklistfilename = sys.argv[2]
    yyyymmddhh = sys.argv[3]
    dayrange = sys.argv[4]
