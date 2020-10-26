#!/usr/bin/env python

import glob
import sys

network_thresh = 1
count_thresh = 30
omb_thresh = 3

data = {}

files = glob.glob("spd10m_OBA/oba_*_01.txt")
files.sort()

# Loop through each file, and store the O and B for each platform
for file in files:
    chunk = file.split("/")[1][12:14]
    if chunk not in ["00", "06", "12", "18"]:
        continue
    #print("Processing %s" %(file))
    lines = open(file, "r").readlines()
    for line in lines[1:]:
        network = line.split()[0]
        platform = line.split()[1]
        if platform not in data:
            data[platform] = []
        ob = line.split()[4]
        back = line.split()[5]
        data[platform].append("%s:%s:%s" %(network, ob, back))

# Now, loop through each station and calculate the mean OMB.  Create
# a blacklist
print("# Creating blacklist")
print("# Rejecting stations with more than %s network" %(network_thresh))
print("# Rejecting stations with less than %s observations" %(count_thresh))
print("# Rejection stations with absolute mean OMB beyond %s" %(omb_thresh))

platforms = data.keys()
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
    if length > network_thresh:
        #print("NOTE:  %s networks found for %s!" %(length, platform))
        print("%s # Too many networks" %(platform))
        continue
    for network in OMB:
        if n[network] < count_thresh:
            print("%s # Only have %s observation(s)" %(platform, n[network]))
            continue
        OMB[network] = OMB[network] / float(n[network])
        if abs(OMB[network]) > omb_thresh:
            #print("mean OMB for %s %s: %s, sample size: %s" 
            #      %(platform, network, OMB[network], n[network]))
            print("%s # Mean OMB is %s" %(platform, OMB[network]))
            continue
