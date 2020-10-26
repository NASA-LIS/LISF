#!/usr/bin/env python

filenames = {
    "galwem10" : "/discover/nobackup/projects/lis_aist17/emkemp/AFWA/retune_galwem10/galwem10_imerger/blacklist.txt",
    "galwem17" : "/discover/nobackup/projects/lis_aist17/emkemp/AFWA/retune_galwem10/galwem17_imerger/blacklist.txt",
    "gfsfv3" : "/discover/nobackup/projects/lis_aist17/emkemp/AFWA/retune_galwem10/gfsfv3_imerger/blacklist.txt",
}


# Collect the flagged biased stations
bad_omb_stations = {}
keys = list(filenames)
keys.sort()
for key in keys:
    lines = open(filenames[key], "r").readlines()
    for line in lines:
        if "Mean OMB" in line:
            station = line.split()[0]
            if station not in bad_omb_stations:
                bad_omb_stations[station] = [key]
            else:
                bad_omb_stations[station].append(key)

# Construct consolidated OMB-based blacklist, including background ids
print("# Consolidated blacklist based on large mean innovations")
stations = list(bad_omb_stations)
stations.sort()
for station in stations:
    msg = "%s #" %(station)
    for bkg in bad_omb_stations[station]:
        msg += " %s" %(bkg)
    print(msg)



