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
    txt = \
        f"Usage: {sys.argv[0]} CFGFILENAME BLACKLISTFILENAME ENDDATETIME " + \
        "DAYRANGE"
    print(txt)

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
        print(f"[ERR] Config file {cfgfile} does not exist!")
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
        print(f"[ERR] Directory {data_directory} does not exist!")
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
    with open(blacklistfilename, "w", encoding="ascii") as outfile:

        # Build list of OBA files
        txt = f"{data_directory}/oba_*_{data_frequency}.txt"
        files = glob.glob(txt)
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
            with open(file, "r", encoding="ascii") as infile:
                lines = infile.readlines()
            for line in lines[1:]:
                network = line.split()[0]
                platform = line.split()[1]
                if is_sat(platform):
                    continue
                if platform not in data:
                    data[platform] = []
                obs = line.split()[4]
                back = line.split()[5]
                data[platform].append(f"{network}:{obs}:{back}")

        # Now, loop through each station and calculate the mean OMB.
        # Create a blacklist
        txt = "# Rejecting stations with more than " + \
            f"{network_thresh} network\n"
        outfile.write(txt)
        txt = "# Rejecting stations with less than " + \
            f"{count_thresh} observations\n"
        outfile.write(txt)
        txt = "# Rejecting stations with absolute mean OMB beyond " + \
            f"{omb_thresh}\n"
        outfile.write(txt)

        platforms = list(data.keys())
        platforms.sort()
        for platform in platforms:
            num = {}
            omb = {}
            entries = data[platform]
            for entry in entries:
                network = entry.split(":")[0]
                obs = float(entry.split(":")[1])
                back = float(entry.split(":")[2])
                if network not in num:
                    num[network] = 1
                else:
                    num[network] += 1
                if network not in omb:
                    omb[network] = 0.
                omb[network] += (obs - back)
            length = len(omb.keys())
            if length > int(network_thresh):
                outfile.write(f"{platform} # Found in {length} networks\n")
                continue
            for network in omb:
                if num[network] < int(count_thresh):
                    txt = \
                      f"{platform} # Only have {num[network]} observation(s)\n"
                    outfile.write(txt)
                    continue
                omb[network] = omb[network] / float(num[network])
                if abs(omb[network]) > float(omb_thresh):
                    txt = f"{platform} # Mean OMB is {omb[network]}\n"
                    outfile.write(txt)
                    continue

#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Check command line
    if len(sys.argv) != 5:
        print("[ERR] Bad command line arguments!")
        usage()
        sys.exit(1)
