#!/usr/bin/env python3

"""
SCRIPT: calc_gfsgalwem_biasratios_multiyear.py

Calculates bias ratio fields using monthly IMERG-FR V07A, GFS, and GALWEM
data, all already interpolated to the NAFPA grid and summed to monthly
totals via LVT.  Outputs 12 fields total, each covering a month based
on multiple years of data.

REVISION HISTORY:
04 Oct 2023:  Eric Kemp:  Initial specification.
"""

import configparser
import datetime
import os
import sys

from netCDF4 import Dataset as nc4_dataset
import numpy as np

# Start and end date of GFS period for NAFPA
back_source = "GFS"
#startdt = datetime.datetime(2008, 1, 1, 0, 0, 0) # GFS-Spectral
#enddt   = datetime.datetime(2019, 6, 1, 0, 0, 0) # GFS-Spectral
#startdt = datetime.datetime(2019, 7, 1, 0, 0, 0) # GFS-FV3
#enddt   = datetime.datetime(2023, 4, 1, 0, 0, 0) # GFS-FV3
#startdt = datetime.datetime(2008, 1, 1, 0, 0, 0) # All GFS
#enddt   = datetime.datetime(2023, 4, 1, 0, 0, 0) # All GFS
#startdt = datetime.datetime(2008, 1, 1, 0, 0, 0) # Pre-GALWEM
#enddt   = datetime.datetime(2017, 9, 1, 0, 0, 0) # Pre-GALWEM
startdt = datetime.datetime(2017, 10, 1, 0, 0, 0) # GALWEM Era
enddt   = datetime.datetime(2023,  4, 1, 0, 0, 0) # GALWEM Era

topdir_back = "/discover/nobackup/projects/usaf_lis/emkemp/AFWA/lis76_imergf_biascorr/data/GFS_NAFPA_Monthly_all"

# Start and end date of GALWEM period for NAFPA
#back_source = "GALWEM"
#startdt = datetime.datetime(2017, 10, 1, 0, 0, 0)
#enddt = datetime.datetime(2023, 4, 1, 0, 0, 0)
#topdir_back = "/discover/nobackup/projects/usaf_lis/emkemp/AFWA/lis76_imergf_biascorr/data/GALWEM_NAFPA_Monthly_all"

timedelta = datetime.timedelta(days=1)

# Other data
topdir_imergf = "/discover/nobackup/projects/usaf_lis/emkemp/AFWA/lis76_imergf_biascorr/data/IMERGF_V07A_NAFPA_Monthly"
topdir_biasratio = f"testdir_{back_source}"

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
        backoutdir = config.get('Output', 'gfsdir')
    elif backsource == "GALWEM":
        backindir = config.get('Input', 'galwemdir')
        backoutdir = config.get('Output', 'galwemdir')
    imergdir = config.get('Input', 'imergdir')

    return backindir, backoutdir, imergdir

def _calc_biasratios(imergdir, backindir, startdate, enddate):
    """Calculate the bias ratios from the input data files."""

    ny = 1920
    nx = 2560
    nmon = 12

    sum_blended = np.zeros([nmon,ny,nx])
    sum_back = np.zeros([nmon,ny,nx])
    precip_ratio = np.zeros([nmon,ny,nx])

    timedelta = datetime.timedelta(days=1)

    # Loop through each month
    curdate = startdate
    while curdate <= enddate:

        # LVT output files have data from the prior month, so we must
        # advance one month to find the appropriate file, e.g., data
        # for May 2020 will be in SUM_TS.202006010000.d01.nc.  Python's
        # timedelta object isn't smart enough to jump a whole month, so
        # we loop through each day instead.
        nextdate = curdate + timedelta
        while nextdate.day != 1:
            nextdate += timedelta

        # First, the IMERG file
        filename = f"{imergdir}/SUM_TS."
        filename += f"{nextdate.year:04d}{nextdate.month:02d}"
        filename += f"{nextdate.day:02d}0000.d01.nc"
        ncid_imerg = nc4_dataset(filename, mode='r', \
                                 format="NETCDF4_CLASSIC")
        precip_imerg = ncid_imerg.variables['TotalPrecip'][:,:]
        lats_imerg = ncid_imerg.variables["latitude"][:,:]
        lons_imerg = ncid_imerg.variables["longitude"][:,:]
        north_south = ncid_imerg.dimensions['north_south'].size
        east_west = ncid_imerg.dimensions['east_west'].size
        del ncid_imerg

        # Next, the background data (GFS or GALWEM)
        filename = f"{backindir}/SUM_TS."
        filename += f"{nextdate.year:04d}{nextdate.month:02d}"
        filename += f"{nextdate.day:02d}0000.d01.nc"
        ncid_back = nc4_dataset(filename, mode='r', \
                               format="NETCDF4_CLASSIC")
        precip_back = ncid_back.variables['TotalPrecip'][:,:]
        del ncid_back

        # Use IMERG from 40S to 71N; use linear tapers from 40S to 60S,
        # and 51N to 71N; and don't use IMERG south of 60S and north of
        # 71N.  Rationale:  No IR and little gauge data is available
        # south of 60S, so we don't expect IMERG to be useful here.
        # In the northern hemisphere, there is no IR data north of 60N,
        # but there is good GPCC gauge density in Scandinavia up to 71N.
        # Finally, the 20 degree latitude linear taper mimicks MERRA-2
        # usage of CPCU gauge analyses.
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

        # Add trace precipitation every month to prevent undefined
        # ratios in deserts.
        sum_blended[(curdate.month-1),:,:] += precip_blended[:,:] + 0.05
        sum_back[(curdate.month-1),:,:] += precip_back[:,:] + 0.05

        # Move on to next month
        curdate = nextdate

    # Finish calculation
    precip_ratio[:,:,:] = sum_blended[:,:,:] / sum_back[:,:,:]
    precip_ratio[:,:,:] = np.where(precip_ratio == 0, 1, precip_ratio)

    return lats_imerg, lons_imerg, precip_ratio, east_west, north_south

def _create_output_filename(backoutdir, backsource, startdate, enddate):
    """Create output netCDF file name"""
    filename = f"{backoutdir}"
    filename += f"/{backsource}_pcp_biasratios_"
    filename += f"{startdate.year:04d}{startdate.month:02d}_"
    filename += f"{enddate.year:04d}{enddate.month:02d}.nc"
    return filename

def _write_biasratios(args):
    """Write out bias ratios to netCDF file"""

    os.makedirs(args['backoutdir'], exist_ok=True)

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
    months = rootgrp.createDimension("months", 12)
    north_south = rootgrp.createDimension("north_south", \
                                          args['north_south'])
    east_west = rootgrp.createDimension("east_west", \
                                        args['east_west'])

    # Define latitude
    latitude = rootgrp.createVariable("latitude", "f4", \
                                      ("north_south", "east_west",))
    latitude.units = "degree_north"
    latitude.standard_name = "latitude"
    latitude.long_name = "latitude"
    latitude.scale_factor = np.float32("1.")
    latitude.add_offset = np.float32("0.")
    latitude.missing_value = np.float32("-9999.")

    # Define longitude
    longitude = rootgrp.createVariable("longitude", "f4", \
                                       ("north_south", "east_west",))
    longitude.units = "degree_east"
    longitude.standard_name = "longitude"
    longitude.long_name = "longitude"
    longitude.scale_factor = np.float32("1.")
    longitude.add_offset = np.float32("0.")
    longitude.missing_value = np.float32("-9999.")

    # Define biasRatio
    biasRatio = rootgrp.createVariable("biasRatio", "f4", \
                              ("months", "north_south", "east_west",))
    biasRatio.units = "-"
    biasRatio.long_name = \
        f"bias_ratio_for_{args['backsource']}_precipitation"
    biasRatio.scale_factor = np.float32("1.")
    biasRatio.add_offset = np.float32("0.")
    biasRatio.missing_value = np.float32("-9999.")

    latitude[:,:] = args['lats_imerg'][:,:]
    longitude[:,:] = args['lons_imerg'][:,:]
    biasRatio[:,:,:] = args['precip_ratio'][:,:,:]

    rootgrp.close()

def _main():
    """Main driver"""

    cfgfile, backsource, startdate, enddate = _process_cmd_line()
    backindir, backoutdir, imergdir = \
        _process_cfg_file(cfgfile, backsource)
    lats_imerg, lons_imerg, precip_ratio, east_west, north_south = \
        _calc_biasratios(imergdir, backindir, startdate, enddate)
    outfile = _create_output_filename(backoutdir, backsource,
                                      startdate, enddate)
    # To satisfy pylint....
    args = {
        "outfile" : outfile,
        "backoutdir" : backoutdir,
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
