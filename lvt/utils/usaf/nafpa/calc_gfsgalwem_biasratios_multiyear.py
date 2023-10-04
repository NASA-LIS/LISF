#!/usr/bin/env python3

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

def _main():
    """Main driver"""

    ny = 1920
    nx = 2560
    nmon = 12

    sum_blended = np.zeros([nmon,ny,nx])
    sum_back = np.zeros([nmon,ny,nx])
    precip_ratio = np.zeros([nmon,ny,nx])

    # Loop through each month
    curdt = startdt
    while curdt <= enddt:

        # LVT output files have data from the prior month, so we must
        # advance one month to find the appropriate file.  Python's
        # timedelta object isn't smart enough to jump a whole month, so
        # we loop through each day instead.
        nextdt = curdt + timedelta
        while nextdt.day != 1:
            nextdt += timedelta

        # First, the IMERG file
        filename = f"{topdir_imergf}/SUM_TS."
        filename += f"{nextdt.year:04d}{nextdt.month:02d}{nextdt.day:02d}"
        filename += "0000.d01.nc"
        ncid_imerg = nc4_dataset(filename, mode='r', \
                                 format="NETCDF4_CLASSIC")
        precip_imerg = ncid_imerg.variables['TotalPrecip'][:,:]
        lats_imerg = ncid_imerg.variables["latitude"][:,:]
        lons_imerg = ncid_imerg.variables["longitude"][:,:]
        north_south = ncid_imerg.dimensions['north_south'].size
        east_west = ncid_imerg.dimensions['east_west'].size
        del ncid_imerg

        # Next, the GFS or GALWEM file
        filename = f"{topdir_back}/SUM_TS."
        filename += f"{nextdt.year:04d}{nextdt.month:02d}{nextdt.day:02d}"
        filename += "0000.d01.nc"
        ncid_back = nc4_dataset(filename, mode='r', \
                               format="NETCDF4_CLASSIC")
        precip_back = ncid_back.variables['TotalPrecip'][:,:]
        del ncid_back

        # Use IMERG from 40S to 40N, use linear tapers from 40S to 60S,
        # and 40N to 60N, and don't use IMERG poleward of 60N and 60S.
        # RATIONALE:  No IR data is available poleward of 60N and 60S,
        # and PMW data has diminished performance over frozen surfaces.
        # Also, GPCC gage coverage tends to thin out in poleward regions,
        # (no GPCC gages at all in Antarctica), so we screen out these
        # areas for simplicity.  Finally, the 20 degree linear taper
        # mimicks MERRA-2 usage of CPCU gauge analyses.
        precip_imerg_weights = np.ones(np.shape(lats_imerg))
        precip_imerg_weights = np.where(lats_imerg < -60,
                                        0,
                                        precip_imerg_weights)
        precip_imerg_weights = np.where( (lats_imerg >= -60) &
                                         (lats_imerg <  -40),
                                         ( (lats_imerg + 60.) / 20.),
                                         precip_imerg_weights)
        #precip_imerg_weights = np.where( (lats_imerg >  40) &
        #                                (lats_imerg <= 60),
        #                                ( (60. - lats_imerg) / 20.),
        #                                precip_imerg_weights)
        #precip_imerg_weights = np.where(lats_imerg > 60,
        #                              0,
        #                              precip_imerg_weights)
        # EMK Alternative...Move the northern linear taper region to
        # 51N to 71N.  Rationale is to leverage the relatively high
        # GPCC gauge density in the Scandinavian Peninsula.
        precip_imerg_weights = np.where( (lats_imerg >  51) &
                                          (lats_imerg <= 71),
                                         ( (71. - lats_imerg) / 20.),
                                          precip_imerg_weights)
        precip_imerg_weights = np.where(lats_imerg > 71,
                                        0,
                                        precip_imerg_weights)

        precip_blended = precip_imerg_weights[:,:]*precip_imerg[:,:] + \
            (1. - precip_imerg_weights[:,:])*precip_back[:,:]

        sum_blended[(curdt.month-1),:,:] += precip_blended[:,:] + 0.05
        sum_back[(curdt.month-1),:,:] += precip_back[:,:] + 0.05

        # Move on to next month
        curdt = nextdt

    # Output to file
    filename = f"{topdir_biasratio}/{back_source}_pcp_biasratio.nc"

    rootgrp = nc4_dataset(filename, "w", format="NETCDF4")
    rootgrp.missing_value = np.float32("-9999.")
    rootgrp.title = f"Monthly bias ratio for IMERG / {back_source}"
    rootgrp.institution = "NASA GSFC"
    now = datetime.datetime.utcnow()
    history = "created on date: "
    history += f"{now.year:04d}-{now.month:02d}-{now.day:02d}"
    history += f"T{now.hour:02d}:{now.minute:02d}:{now.second:02d}"
    rootgrp.history = f"created on date: {history}"
    rootgrp.comment = "website: http://lis.gsfc.nasa.gov/"
    rootgrp.MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL"
    rootgrp.SOUTH_WEST_CORNER_LAT = np.float32("-89.95312")
    rootgrp.SOUTH_WEST_CORNER_LON = np.float32("-179.9297")
    rootgrp.DX = np.float32("0.140625")
    rootgrp.DY = np.float32("0.09375")

    # Define dimensions
    months = rootgrp.createDimension("months", 12)
    north_south = rootgrp.createDimension("north_south", north_south)
    east_west = rootgrp.createDimension("east_west", east_west)

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
    biasRatio.long_name = f"bias_ratio_for_{back_source}_precipitation"
    biasRatio.scale_factor = np.float32("1.")
    biasRatio.add_offset = np.float32("0.")
    biasRatio.missing_value = np.float32("-9999.")

    latitude[:,:] = lats_imerg[:,:]
    longitude[:,:] = lons_imerg[:,:]

    precip_ratio[:,:,:] = sum_blended[:,:,:] / sum_back[:,:,:]
    precip_ratio[:,:,:] = np.where(precip_ratio == 0, 1, precip_ratio)
    biasRatio[:,:,:] = precip_ratio[:,:,:]

    rootgrp.close()


if __name__ == "__main__":
    _main()
