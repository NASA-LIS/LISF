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
SCRIPT: calc_gfsgalwem_biasratios_multiyear.py

Calculates bias ratio fields using monthly IMERG-FR V07A, GFS, and GALWEM
data, all already interpolated to the NAFPA grid and summed to monthly
totals via LVT.  Outputs 12 fields total, each covering a month based
on multiple years of data.

REVISION HISTORY:
04 Oct 2023:  Eric Kemp:  Initial specification.
06 Oct 2023:  Eric Kemp:  Added log10 transform to output for plotting.
                          Refactored to please pylint.
"""

import configparser
import datetime
import os
import sys

# Disable false alarm in pylint (Dataset is not visible in netCDF4 module
# because it is not implemented in python)
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
import numpy as np

def _usage():
    """Print usage message for this script."""
    print(f"Usage: {sys.argv[0]} CFGFILE BACKSOURCE STARTDATE ENDDATE")
    print("   CFGFILE is path to config file")
    print("   BACKSOURCE is GFS or GALWEM")
    print("   STARTDATE is start date of calculation (YYYYMM)")
    print("   ENDDATE is end date of calculateion (YYYYMM)")

def _process_cmd_line():
    """Process command line arguments."""
    if len(sys.argv) != 5:
        _usage()
        sys.exit(1)

    cfgfile = sys.argv[1]
    if not os.path.exists(cfgfile):
        print(f"[ERR] {cfgfile} does not exist!")
        sys.exit(1)

    backsource = sys.argv[2]
    if backsource not in ["GFS", "GALWEM"]:
        print(f"[ERR] {backsource} is not a valid background source")
        print("[ERR] Only GFS and GALWEM recognized")
        sys.exit(1)

    yyyymm = sys.argv[3]
    year = int(yyyymm[0:4])
    month = int(yyyymm[4:6])
    startdate = datetime.datetime(year, month, 1)

    yyyymm = sys.argv[4]
    year = int(yyyymm[0:4])
    month = int(yyyymm[4:6])
    enddate = datetime.datetime(year, month, 1)

    if startdate > enddate:
        print("[ERR] STARTDATE is beyond ENDDATE!")
        sys.exit(1)

    return cfgfile, backsource, startdate, enddate

def _process_cfg_file(cfgfile, backsource):
    """Processes config file for this script."""
    config = configparser.ConfigParser()
    config.read(cfgfile)

    if backsource == "GFS":
        backindir = config.get('Input', 'gfsdir')
    elif backsource == "GALWEM":
        backindir = config.get('Input', 'galwemdir')
    imergdir = config.get('Input', 'imergdir')
    outdir = config.get('Output', 'outdir')

    return backindir, imergdir, outdir

def _create_arrays():
    """Allocate and zero out numpy arrays"""
    nlat = 1920
    nlon = 2560
    nmon = 12
    sum_blended = np.zeros([nmon,nlat,nlon])
    sum_back = np.zeros([nmon,nlat,nlon])
    precip_ratio = np.zeros([nmon,nlat,nlon])
    return sum_blended, sum_back, precip_ratio

def _read_imerg_file(imergdir, nextdate):
    """Read the interpolated monthly IMERG file."""
    filename = f"{imergdir}/SUM_TS."
    filename += f"{nextdate.year:04d}{nextdate.month:02d}"
    filename += f"{nextdate.day:02d}0000.d01.nc"
    ncid_imerg = nc4_dataset(filename, mode='r', \
                             format="NETCDF4_CLASSIC")
    precip_imerg = ncid_imerg.variables['TotalPrecip'][:,:]
    del ncid_imerg
    return precip_imerg

def _read_imerg_file_latlondims(imergdir, startdate):
    """Read the interpolated monthly IMERG file."""
    nextdate = _get_nextdate(startdate)
    filename = f"{imergdir}/SUM_TS."
    filename += f"{nextdate.year:04d}{nextdate.month:02d}"
    filename += f"{nextdate.day:02d}0000.d01.nc"
    ncid_imerg = nc4_dataset(filename, mode='r', \
                             format="NETCDF4_CLASSIC")
    lats_imerg = ncid_imerg.variables["latitude"][:,:]
    lons_imerg = ncid_imerg.variables["longitude"][:,:]
    north_south = ncid_imerg.dimensions['north_south'].size
    east_west = ncid_imerg.dimensions['east_west'].size
    del ncid_imerg
    return lats_imerg, lons_imerg, north_south, east_west

def _read_background_file(backindir, nextdate):
    """Read interpolated monthly background file."""
    filename = f"{backindir}/SUM_TS."
    filename += f"{nextdate.year:04d}{nextdate.month:02d}"
    filename += f"{nextdate.day:02d}0000.d01.nc"
    ncid_back = nc4_dataset(filename, mode='r', \
                            format="NETCDF4_CLASSIC")
    precip_back = ncid_back.variables['TotalPrecip'][:,:]
    del ncid_back
    return precip_back

def _calc_precip_blended(lats_imerg, precip_imerg, precip_back):
    """Calculate weighted blend of IMERG and background precip"""

    # Use IMERG from 40S to 71N; use linear tapers from 40S to 60S,
    # and 51N to 71N; and don't use IMERG south of 60S and north of
    # 71N.  Rationale:  No IR and little gauge data is available
    # south of 60S, so we don't expect IMERG to be useful here for
    # bias correction.  In the northern hemisphere, there is no IR data
    # north of 60N, but there is good GPCC gauge density in Scandinavia
    # up to 71N.  Finally, the 20 degree latitude linear tapers mimicks
    # MERRA-2 usage of CPCU gauge analyses.
    precip_imerg_weights = np.ones(np.shape(lats_imerg))
    precip_imerg_weights = np.where(lats_imerg < -60,
                                    0,
                                    precip_imerg_weights)
    precip_imerg_weights = np.where( (lats_imerg >= -60) &
                                     (lats_imerg <  -40),
                                     ( (lats_imerg + 60.) / 20.),
                                     precip_imerg_weights)
    precip_imerg_weights = np.where( (lats_imerg >  51) &
                                     (lats_imerg <= 71),
                                     ( (71. - lats_imerg) / 20.),
                                     precip_imerg_weights)
    precip_imerg_weights = np.where(lats_imerg > 71,
                                    0,
                                    precip_imerg_weights)

    # Conservative alternative:  Linear taper from 40N to 60N, and
    # screen out everything north of 60N.  We'll do this as a backup
    # if we find unphysical IMERG patterns north of 60N
    #precip_imerg_weights = np.where( (lats_imerg >  40) &
    #                                (lats_imerg <= 60),
    #                                ( (60. - lats_imerg) / 20.),
    #                                precip_imerg_weights)
    #precip_imerg_weights = np.where(lats_imerg > 60,
    #                              0,
    #                              precip_imerg_weights)

    precip_blended = precip_imerg_weights[:,:]*precip_imerg[:,:] + \
        (1. - precip_imerg_weights[:,:])*precip_back[:,:]

    return precip_blended

def _get_nextdate(curdate):
    """Determine next date (first of next month)"""
    # LVT output files have data from the prior month, so we must
    # advance one month to find the appropriate file, e.g., data
    # for May 2020 will be in SUM_TS.202006010000.d01.nc.  Python's
    # timedelta object isn't smart enough to jump a whole month, so
    # we loop through each day instead.
    timedelta = datetime.timedelta(days=1)
    nextdate = curdate + timedelta
    while nextdate.day != 1:
        nextdate += timedelta
    return nextdate

def _calc_biasratios(imergdir, backindir, startdate, enddate, lats_imerg):
    """Calculate the bias ratios from the input data files."""

    sum_blended, sum_back, precip_ratio = _create_arrays()

    # Loop through each month
    curdate = startdate
    while curdate <= enddate:

        # LVT output files have data from the prior month, so we must
        # advance one month to find the appropriate file, e.g., data
        # for May 2020 will be in SUM_TS.202006010000.d01.nc.
        nextdate = _get_nextdate(curdate)

        # First, the IMERG file
        precip_imerg = \
            _read_imerg_file(imergdir, nextdate)

        # Next, the background data (GFS or GALWEM)
        precip_back = _read_background_file(backindir, nextdate)

        # Calculate blended precipitation (weighted average of IMERG
        # and background, varying by latitude).
        precip_blended = \
            _calc_precip_blended(lats_imerg, precip_imerg, precip_back)

        # Updated precip sums.  Add trace precipitation every month to
        # prevent undefined ratios in deserts.
        sum_blended[(curdate.month-1),:,:] += precip_blended[:,:] + 0.05
        sum_back[(curdate.month-1),:,:] += precip_back[:,:] + 0.05

        # Move on to next month
        curdate = nextdate

    # Finish calculation
    precip_ratio[:,:,:] = sum_blended[:,:,:] / sum_back[:,:,:]
    precip_ratio[:,:,:] = np.where(precip_ratio == 0, 1, precip_ratio)

    return precip_ratio

def _create_output_filename(outdir, backsource, startdate, enddate):
    """Create output netCDF file name"""
    filename = f"{outdir}"
    filename += f"/{backsource}_pcp_biasratios_"
    filename += f"{startdate.year:04d}{startdate.month:02d}_"
    filename += f"{enddate.year:04d}{enddate.month:02d}.nc"
    return filename

def _create_latitude(rootgrp):
    """Create latitude dataset in output file"""
    latitude = rootgrp.createVariable("latitude", "f4", \
                                      ("north_south", "east_west",))
    latitude.units = "degree_north"
    latitude.standard_name = "latitude"
    latitude.long_name = "latitude"
    latitude.scale_factor = np.float32("1.")
    latitude.add_offset = np.float32("0.")
    latitude.missing_value = np.float32("-9999.")
    return latitude

def _create_longitude(rootgrp):
    """Create longitude dataset in output file."""
    longitude = rootgrp.createVariable("longitude", "f4", \
                                       ("north_south", "east_west",))
    longitude.units = "degree_east"
    longitude.standard_name = "longitude"
    longitude.long_name = "longitude"
    longitude.scale_factor = np.float32("1.")
    longitude.add_offset = np.float32("0.")
    longitude.missing_value = np.float32("-9999.")
    return longitude

def _create_bias_ratio(rootgrp, backsource):
    """Create the biasRatio dataset in the output file."""
    bias_ratio = rootgrp.createVariable("biasRatio", "f4", \
                              ("months", "north_south", "east_west",))
    bias_ratio.units = "-"
    bias_ratio.long_name = \
        f"bias_ratio_for_{backsource}_precipitation"
    bias_ratio.scale_factor = np.float32("1.")
    bias_ratio.add_offset = np.float32("0.")
    bias_ratio.missing_value = np.float32("-9999.")
    return bias_ratio

def _create_log10_bias_ratio(rootgrp, backsource):
    """Create the log10BiasRatio dataset in the output file."""
    log10_bias_ratio = rootgrp.createVariable("log10BiasRatio", "f4", \
                                ("months", "north_south", "east_west",))
    log10_bias_ratio.units = "-"
    log10_bias_ratio.long_name = \
        f"log10_bias_ratio_for_{backsource}_precipitation"
    log10_bias_ratio.scale_factor = np.float32("1.")
    log10_bias_ratio.add_offset = np.float32("0.")
    log10_bias_ratio.missing_value = np.float32("-9999.")
    return log10_bias_ratio

def _write_biasratios(args):
    """Write out bias ratios to netCDF file"""

    os.makedirs(args['outdir'], exist_ok=True)

    now = datetime.datetime.utcnow()

    history = "created on date: "
    history += f"{now.year:04d}-{now.month:02d}-{now.day:02d}"
    history += f"T{now.hour:02d}:{now.minute:02d}:{now.second:02d}"

    rootgrp = nc4_dataset(args['outfile'], "w", format="NETCDF4")
    rootgrp.missing_value = np.float32("-9999.")
    rootgrp.title = f"Monthly bias ratio for IMERG / {args['backsource']}"
    rootgrp.institution = "NASA GSFC"

    rootgrp.history = f"created on date: {history}"
    rootgrp.comment = "website: http://lis.gsfc.nasa.gov/"
    rootgrp.MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL"
    rootgrp.SOUTH_WEST_CORNER_LAT = np.float32("-89.95312")
    rootgrp.SOUTH_WEST_CORNER_LON = np.float32("-179.9297")
    rootgrp.DX = np.float32("0.140625")
    rootgrp.DY = np.float32("0.09375")

    # Define dimensions
    rootgrp.createDimension("months", 12)
    rootgrp.createDimension("north_south", \
                            args['north_south'])
    rootgrp.createDimension("east_west", \
                            args['east_west'])

    # Define output variables
    latitude = _create_latitude(rootgrp)
    longitude = _create_longitude(rootgrp)
    bias_ratio = _create_bias_ratio(rootgrp, args['backsource'])
    log10_bias_ratio = _create_log10_bias_ratio(rootgrp,
                                                args['backsource'])

    latitude[:,:] = args['lats_imerg'][:,:]
    longitude[:,:] = args['lons_imerg'][:,:]
    bias_ratio[:,:,:] = args['precip_ratio'][:,:,:]
    log10_bias_ratio[:,:,:] = np.log10(args['precip_ratio'][:,:,:])

    rootgrp.close()

def _main():
    """Main driver"""

    cfgfile, backsource, startdate, enddate = _process_cmd_line()
    backindir, imergdir, outdir = \
        _process_cfg_file(cfgfile, backsource)
    lats_imerg, lons_imerg, north_south, east_west = \
        _read_imerg_file_latlondims(imergdir, startdate)
    precip_ratio = \
        _calc_biasratios(imergdir, backindir, startdate, enddate, \
                         lats_imerg)
    outfile = _create_output_filename(outdir, backsource,
                                      startdate, enddate)
    # To satisfy pylint....
    args = {
        "outfile" : outfile,
        "outdir" : outdir,
        "backsource" : backsource,
        "north_south" : north_south,
        "east_west" : east_west,
        "lats_imerg" : lats_imerg,
        "lons_imerg" : lons_imerg,
        "precip_ratio" : precip_ratio,
    }
    _write_biasratios(args)

if __name__ == "__main__":
    _main()
