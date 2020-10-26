#!/usr/bin/env python

# Standard libraries
import configparser
import datetime
import glob
import math
import os
import sys

# Other libraries
import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------------
def usage():
    print("Usage: %s config.cfg" %(sys.argv[0]))

#------------------------------------------------------------------------------
# Welford algorithm for calculating mean and variance
def update(existingAggregate, newValue):
    (count, mean, M2) = existingAggregate
    count = count + 1
    delta = newValue - mean
    mean = mean + delta / count
    delta2 = newValue - mean
    M2 = M2 + delta * delta2
    return (count, mean, M2)
def finalize(existingAggregate):
    (count, mean, M2) = existingAggregate
    if count < 2:
        sampleStdDev = float('nan')
    else:
        (mean, sampleVariance) = (mean, M2/(count-1))
        sampleStdDev = np.sqrt(sampleVariance)
    return (mean, sampleStdDev)

#------------------------------------------------------------------------------
def is_gage(network):
    answer = False
    if network in ["AMIL", "CANA", "FAA", "ICAO", "WMO", "MOBL", "SUPERGAGE"]:
        answer = True
    return answer
def is_ssmi(network):
    answer = False
    if network in ["SSMI"]:
        answer = True
    return answer
def is_geoprecip(network):
    answer = False
    if network in ["GEOPRECIP"]:
        answer = True
    return answer
def is_cmorph(network):
    answer = False
    if network in ["CMORPH"]:
        answer = True
    return answer
def is_imerg(network):
    answer = False
    if network in ["IMERG"]:
        answer = True
    return answer

obtype_func_dict = {
    "Gage" : is_gage,
    "SSMI"  : is_ssmi,
    "GEOPRECIP" : is_geoprecip,
    "CMORPH" : is_cmorph,
    "IMERG" : is_imerg,
}

#------------------------------------------------------------------------------
# Check command line
if len(sys.argv) != 2:
    print("[ERR] Bad command line arguments!")
    usage()
    sys.exit(1)

# Read config file
cfgfile = sys.argv[1]
if not os.path.exists(cfgfile):
    print("[ERR] Config file %s does not exist!" %(cfgfile))
    sys.exit(1)
config = configparser.ConfigParser()
config.read(cfgfile)

# Handle the Input section
try:
    datadir = config.get('Input', 'datadir')
    obtype = config.get('Input', 'obtype')
except:
    print("[ERR] Problem reading from config file!")
    raise
if not os.path.exists(datadir):
    print("[ERR] %s does not exist!" %(datadir))
    sys.exit(1)

obtypes = list(obtype_func_dict.keys())
obtypes.sort()
if not obtype in obtypes:
    print("[ERR] Observation type %s is not supported!" %(obtype))
    print("Currently only the following observations are supported:")
    for type in obtypes:
        print("   %s" %(type))
    sys.exit(1)
obtype_func = obtype_func_dict[obtype]

try:
    time_gap_hrs = config.getint('Input', 'time_gap_hrs')
except:
    print("[ERR] Problem reading from config file!")
    raise
if time_gap_hrs < 3:
    print("[ERR] time_gap_hrs must be at least 3!")
    sys.exit(1)

# Handle the Fit section
try:
    function_type = config.get('Fit', 'function_type')
    degree = config.getint('Fit', 'degree')
    deltax = config.getfloat('Fit', 'deltax')
    mincount = config.getint('Fit', 'mincount')
except:
    print("[ERR] Problem reading from config file!")
    raise
function_types = ["Polynomial"]
if function_type not in function_types:
    print("[ERR] Function type %s is not supported!" %(function_type))
    print("Currently only the following functions are supported:")
    for type in function_types:
        print("   %" %(type))
    sys.exit(1)
if degree < 1:
    print("[ERR] Polynomial degree must be positive!")
    sys.exit(1)
if deltax <= 0:
    print("[ERR] deltax must be positive!")
    sys.exit(1)
if mincount < 1:
    print("[ERR] mincount must be positive!")
    sys.exit(1)

# Handle the Plot section
try:
    title = config.get('Plot', 'title')
except:
    print("[ERR] Problem reading from config file!")
    raise

# Handle the Blacklist section
try:
    apply_blacklist = config.get('Blacklist', 'apply_blacklist')
    if apply_blacklist == 1:
        blacklist = config.get('Blacklist', 'filename')
except:
    print("[ERR] Problem reading from config file!")
    raise
if apply_blacklist == 1:
    if not os.path.exists(blacklist):
        print("[ERR] %s does not exist!" %(blacklist))
        sys.exit(1)
bad_stations = []
if apply_blacklist == 1:
    lines = open(blacklist, "r").readlines()
    for line in lines:
        if line[0] == "#":
            continue
        bad_stations.append(line.split()[0])

# Initializations
OVector = []
BVector = []
meanOBVector = []
meanOB_max = 0 # meanOB_min must be -3
B_max = 0
O_max = 0

# Loop through precip OBA files.
files = glob.glob("%s/*_12.txt" %(datadir))
files.sort()
first = True
file_count = 0
for file in files:
    yyyymmddhh = file.split("/")[-1].split("_")[1]
    yyyy = int(yyyymmddhh[0:4])
    mm   = int(yyyymmddhh[4:6])
    dd   = int(yyyymmddhh[6:8])
    hh   = int(yyyymmddhh[8:10])
    curdt = datetime.datetime(yyyy, mm, dd, hh)
    process = False
    if first:
        first = False
        process = True
    else:
        if curdt == enddt:
            process = True
    if process:
        print("Processing %s" %(file))
        startdt = curdt
        enddt = startdt + datetime.timedelta(hours=time_gap_hrs)
        file_count += 1
    else:
        continue
    lines = open(file).readlines()
    for line in lines:
        if "#" in line:
            continue
        network = line.split()[0]
        #if not is_gage(network):
        if not obtype_func(network):
            continue
        if apply_blacklist == 1:
            if line.split()[1] in bad_stations:
                print("Skipping blacklisted station ", line.split()[1])
                continue
        # Save the O and B values, and the average of the two.  Also update
        # the maximums
        O = float(line.split()[4])
        B = float(line.split()[5])
        OVector.append(O)
        BVector.append(B)
        #meanOBVector.append(0.5*(O + 3 + B + 3)) # Shift so no precip is zero
        meanOBVector.append(0.5*(O + B)) # Use raw precip
        if meanOBVector[-1] > meanOB_max:
            meanOB_max = meanOBVector[-1]
        if BVector[-1] > B_max:
            B_max = BVector[-1]
        if OVector[-1] > O_max:
            O_max = OVector[-1]

print("len(OVector) = ", len(OVector))
print("Maximum Background value: ", B_max)
print("Maximum Observed value: ", O_max)
print("Maximum average of shifted O and B: ", meanOB_max)


#------------------------------------------------------------------------------
# From Lopez
# Collect the innovations (O-B) as a function of mean(O+B)
# Some binning is used here.
#xrange = np.arange(-3, 20 , 0.25)
xrange = np.arange(0.5*deltax, math.ceil(meanOB_max), deltax)
length = len(xrange)
#print("xrange: ", xrange)

meanOMBs = np.zeros(shape=(length)) # Mean innovations
m2OMBs   = np.zeros(shape=(length))
stdDevOMBs = np.zeros(shape=(length)) # Standard deviations of innovations
icounts = np.zeros(shape=(length))
print("len(meanOBVector) = ", len(meanOBVector))
for i in range(0, len(meanOBVector)):
    OMB = OVector[i] - BVector[i]  # Innovation
    meanOB = meanOBVector[i]           # Average of shifted O + B
    # Find the appropriate bin for mean(O+B)
    if meanOB == 0:
        #index = 0
        continue
    else:
        index = int(np.ceil(meanOB/deltax)) - 1
    try:
#        meanOMBs[index] += OMB
#        STDevOMBs[index] += OMB
#        icounts[index] += 1
        # Use Welford algorithm
        aggregate = (icounts[index], meanOMBs[index], m2OMBs[index])
        (icounts[index], meanOMBs[index], m2OMBs[index]) = \
            update(aggregate, OMB)
    except:
        continue

# Now find the mean innovation for each mean(O+B)
for i in range(0, length):
    #if icounts[i] > 0:
    #    meanOMBs[i] = meanOMBs[i] / float(icounts[i])
    #else:
    #    meanOMBs[i] = np.nan
    # Use Welford algorithm
    aggregate = (icounts[i], meanOMBs[i], m2OMBs[i])
    meanOMBs[i], stdDevOMBs[i] = finalize(aggregate)
    if icounts[i] < mincount:
        meanOMBs[i] = np.nan
        stdDevOMBs[i] = np.nan
        continue
        #if not np.isnan(stdDevOMBs[i]):
        #    stderr = stdDevOMBs[i]/math.sqrt(float(icounts[i]))
        #    if (meanOMBs[i] + 1.96*stderr) > 0 and \
        #       (meanOMBs[i] - 1.96*stderr) < 0:
        #        meanOMBs[i] = np.nan
        #        stdDevOMBs[i] = np.nan

# Count number of bins with actual values
count_real = 0
for i in range(0, length):
    if icounts[i] > 0:
        print(i, icounts[i], xrange[i], meanOMBs[i], \
            stdDevOMBs[i]/math.sqrt(icounts[i]))
    if not np.isnan(meanOMBs[i]):
        count_real += 1
n = len(meanOBVector)
print("Sum of icounts: ", sum(icounts))
print("Full sample size: ", n)
print("max meanOB: ", meanOB_max)
print("Number of files used: ", file_count)

# Collect actual values in new arrays for fitting and plotting
meanOMBs_nonan = np.zeros(shape=(count_real))
xrange_nonan   = np.zeros(shape=(count_real))
j = 0
for i in range(0, length):
    if not np.isnan(meanOMBs[i]):
        meanOMBs_nonan[j] = meanOMBs[i]
        xrange_nonan[j] = xrange[i]
        j += 1

print("xrange_nonan = ", xrange_nonan)

#for i in range(length-1,-1,-1):
#    if not np.isnan(meanOMBs[i]):
#        break
#print i

# Plot the mean innovations as a function of mean(O+B)
if function_type == "Polynomial":
    if degree == 3:
        # Force polynomial to have zero intercept
        x = xrange_nonan
        y = meanOMBs_nonan
        #XX = np.vstack((x**3,x**2,x,np.ones_like(x))).T
        #p = np.linalg.lstsq(XX[:, :-1], y)[0]
        #y_fit = np.dot(p,XX[:,:-1].T)
        XX = np.vstack((x**3, x**2, x)).T
        p = np.linalg.lstsq(XX, y)[0]
        y_fit = np.dot(p, XX.T)
        coeffs = np.zeros(shape=degree+1)
        print("p = ", p)
        for i in range(0, degree):
            coeffs[i] = p[i]
        coeffs[-1] = 0
        p = np.poly1d(coeffs)
        print(p)
        xdummy = np.arange(0, x[-1], 0.01)
        #plt.plot(x,y    ,'bD',
        #         xdummy,p(xdummy),'r-')
        plt.plot(x,      y        , 'bD-',
                 xdummy, p(xdummy), 'r-')
        params = r"$y = %fx^3 + (%f)x^2 + (%f)x$" \
            %(coeffs[0], coeffs[1], coeffs[2])
        plt.figtext(0.35, 0.15, params)

    elif degree == 4:
        # Force polynomial to have zero intercept
        x = xrange_nonan
        y = meanOMBs_nonan
        #XX = np.vstack((x**4,x**3,x**2,x,np.ones_like(x))).T
        #p = np.linalg.lstsq(XX[:, :-1], y)[0]
        #y_fit = np.dot(p,XX[:,:-1].T)
        XX = np.vstack((x**4, x**3, x**2, x)).T
        p = np.linalg.lstsq(XX, y)[0]
        y_fit = np.dot(p, XX.T)
        coeffs = np.zeros(shape=degree+1)
        for i in range(0,degree):
            coeffs[i] = p[i]
        coeffs[-1] = 0
        p = np.poly1d(coeffs)
        print(p)
        xdummy = np.arange(0, x[-1], 0.01)
        plt.plot(x,      y        , 'bD-',
                 xdummy, p(xdummy), 'r-')
        params = r"$y = %fx^4 + (%f)x^3 + (%f)x^2 + (%f)x$" \
            %(coeffs[0], coeffs[1], coeffs[2], coeffs[3])
        #plt.figtext(0.35, 0.15, params)
        plt.figtext(0.05, 0.95, params)

    elif degree == 5:
        # Force polynomial to have zero intercept
        x = xrange_nonan
        y = meanOMBs_nonan
        #XX = np.vstack((x**5,x**4,x**3,x**2,x,np.ones_like(x))).T
        #p = np.linalg.lstsq(XX[:, :-1], y)[0]
        #y_fit = np.dot(p,XX[:,:-1].T)
        XX = np.vstack((x**5, x**4, x**3, x**2, x)).T
        p = np.linalg.lstsq(XX, y)[0]
        y_fit = np.dot(p, XX.T)
        coeffs = np.zeros(shape=degree+1)
        for i in range(0, degree):
            coeffs[i] = p[i]
        coeffs[-1] = 0
        p = np.poly1d(coeffs)
        print(p)
        xdummy = np.arange(0, x[-1], 0.01)
        plt.plot(x,      y        ,'bD-',
                 xdummy, p(xdummy),'r-')
        params = r"$y = %fx^5 + (%f)x^4 + (%f)x^3 + (%f)x^2 + (%f)x$" \
            %(coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4])
        plt.figtext(0.35, 0.15, params)

    else:
        print("WARNING, using generic polynomial fitting for degree ", degree)
        p = np.poly1d(np.polyfit(xrange_nonan, meanOMBs_nonan, degree))
        print(p)
        plt.plot(xrange_nonan, meanOMBs_nonan, 'bD-', \
                     xrange_nonan, p(xrange_nonan),'r-')

# Finish the plot
#params = \
#    r"$x = T + 3, T = t^{-1} (P^{t} - 1), t = 1/3$"
#plt.figtext(0.35, 0.10, params)
plt.title(title)
#plt.xlabel('Mean of O and B Box-Cox Shifted Transformed Precip (T + 3)')
plt.xlabel('Mean of O and B Precipitation (mm)')

#plt.xticks([x for x in range(0, int(math.ceil(meanOB_max))+1)])
#plt.xticks([x for x in range(0, int(math.ceil(xrange_nonan[-1]))+1)])

#plt.ylabel('Mean Innovation (Box-Cox Transformed Precip)')
plt.ylabel('Mean Innovation (mm)')
plt.grid(True)
plt.show()

