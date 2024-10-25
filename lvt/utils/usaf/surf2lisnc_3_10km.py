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

''' 
Created on June 10, 2019
@author: Shugong Wang
@email: shugong.wang@nasa.gov
'''

import os
import mule
import numpy.ma as ma
import netCDF4 as nc4
import datetime as dt


class surf2lisnc(object):
    '''
    surf2lisnc: taking JULES state variables from UM SURF file
    and write into a NetCDF file of LIS output format 
    '''

    def __init__(self, nc_file, stash_file):
        self.nc_file = nc_file
        self.stashmaster = stash_file
        self.nrow = 1920
        self.ncol = 2560
        self.surf_var = ma.zeros((self.nrow, self.ncol))
        self.lat = ma.zeros((self.nrow, self.ncol))
        self.lon = ma.zeros((self.nrow, self.ncol))

        if os.path.isfile(self.nc_file) == True:
            os.remove(self.nc_file)
            print("File %s has been deleted. A new file will be created." %
                  self.nc_file)
        print("creating LIS style NetCDF file {0}".format(self.nc_file))
        self.rootgrp = nc4.Dataset(self.nc_file, mode="w")
        self.lat = self.rootgrp.createDimension("north_south", 1920)
        self.lon = self.rootgrp.createDimension("east_west", 2560)
        self.smcprof = self.rootgrp.createDimension("SoilMoist_profiles", 4)
        self.time = self.rootgrp.createDimension("time", 1)

    def read_surf_variable(self, section, item, lblev=-1, lbuser5=-1):
        if lblev >= 0 and lbuser5 == -1:
            for field in self.surf.fields:
                if field.stash.section == section and field.stash.item == item and field.lblev == lblev:
                    info = "section=%03d item=%03d lblev=%02d %s" % (
                        field.stash.section, field.stash.item, field.lblev, field.stash.name)
                    print(info)
                    self.surf_var[:, :] = field.get_data()[:, :]
                    break
        elif lblev == -1 and lbuser5 >= 0:
            for field in self.surf.fields:
                if field.stash.section == section and field.stash.item == item and field.lbuser5 == lbuser5:
                    info = "section=%03d item=%03d lbuser5=%02d %s" % (
                        field.stash.section, field.stash.item, field.lbuser5, field.stash.name)
                    print(info)
                    self.surf_var[:, :] = field.get_data()[:, :]
                    break
        elif lblev == -1 and lbuser5 == -1:
            for field in self.surf.fields:
                if field.stash.section == section and field.stash.item == item:
                    info = "section=%03d item=%03d %s" % (
                        field.stash.section, field.stash.item,  field.stash.name)
                    print(info)
                    self.surf_var[:, :] = field.get_data()[:, :]
                    break

        # shift coordinate
        data = self.surf_var
        col1 = self.ncol//2
        col2 = self.ncol
        tmp = ma.zeros((self.nrow, self.ncol))
        tmp[:, :] = data[:, :]

        data[:, 0:col1] = tmp[:, col1:col2]
        data[:, col1:col2] = tmp[:, 0:col1]

    def read_surf_smc(self, surf_file):
        layer_thickness = [100.0, 250.0, 650.0, 2000.0]
        print("reading STASH master file {0}".format(self.stashmaster))
        sm = mule.stashmaster.STASHmaster.from_file(self.stashmaster)

        print("reading UM surf file {0}".format(surf_file))
        self.surf = mule.AncilFile.from_file(surf_file)
        print("attaching the STASH master file {0} to the mule reader".format(
            self.stashmaster))
        self.surf.attach_stashmaster_info(sm)
        # section=000 item=009 lblev=01 lbuser4=0009 SOIL MOISTURE CONTENT IN A LAYER
        smc = ma.zeros((4, self.nrow, self.ncol))
        for k in range(0, 4):
            self.read_surf_variable(0, 9, lblev=k+1)
            smc[k, :, :] = self.surf_var[:, :]/layer_thickness[k]

        print(smc[0, :, :].max())
        # reset no_data value
        smc[smc < -9999] = -9999
        nc_data = self.rootgrp.createVariable(
            "SoilMoist_inst", "f4", ("SoilMoist_profiles", "north_south", "east_west"))
        nc_data[:, :, :] = smc[:, :, :]

    def read_surf_swe(self, surf_file):
        print("reading STASH master file {0}".format(self.stashmaster))
        sm = mule.stashmaster.STASHmaster.from_file(self.stashmaster)

        print("reading UM surf file {0}".format(surf_file))
        self.surf = mule.AncilFile.from_file(surf_file)
        print("attaching the STASH master file {0} to the mule reader".format(
            self.stashmaster))
        self.surf.attach_stashmaster_info(sm)

        # 9 SNOW_MASS_IJ, section=000 item=023 lblev=9999 lbuser4=0023 SNOW AMOUNT OVER LAND AFT TSTP KG/M2
        swe = ma.zeros((self.nrow, self.ncol))
        self.read_surf_variable(0, 23)
        swe[:, :] = self.surf_var[:, :]
        swe[swe < -9999] = -9999
        nc_data = self.rootgrp.createVariable(
            "SWE_inst", "f4", ("north_south", "east_west"))
        nc_data[:, :] = swe[:, :]

    def write_lat_lon(self, lat, lon):
        nc_lat = self.rootgrp.createVariable(
            "lat", "f4", ("north_south", "east_west"))
        nc_lon = self.rootgrp.createVariable(
            "lon", "f4", ("north_south", "east_west"))
        nc_lat[:, :] = lat[:, :]
        nc_lon[:, :] = lon[:, :]

    def write_attributes(self):
        # missing_value = -9999.f ;
        self.rootgrp.missing_value = -9999.0
        #NUM_SOIL_LAYERS = 4 ;
        self.rootgrp.NUM_SOIL_LAYERS = 4
        # SOIL_LAYER_THICKNESSES = 10.f, 25.f, 65.f, 200.f ;
        self.rootgrp.SOIL_LAYER_THICKNESSES = [10.0, 25.0, 65.0, 200.0]
        #title = "LIS land surface model output" ;
        self.rootgrp.title = "LIS land surface model output"
        #institution = "NASA GSFC" ;
        self.rootgrp.institution = "NASA GSFC"
        #source = "+template open water" ;
        self.rootgrp.source = "+template open water"
        #history = "created on date: 2019-06-11T16:44:58.636" ;
        self.rootgrp.history = "created on date: 2019-06-11T16:44:58.636"
        #references = "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007" ;
        self.rootgrp.references = "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"
        #conventions = "CF-1.6" ;
        self.rootgrp.conventions = "CF-1.6"
        #comment = "website: http://lis.gsfc.nasa.gov/" ;
        self.rootgrp.comment = "website: http://lis.gsfc.nasa.gov/"
        #MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL" ;
        self.rootgrp.MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL"
        # SOUTH_WEST_CORNER_LAT = -89.95312f ;
        self.rootgrp.SOUTH_WEST_CORNER_LAT = -89.95312
        # SOUTH_WEST_CORNER_LON = -179.9297f ;
        self.rootgrp.SOUTH_WEST_CORNER_LON = -179.9297
        # DX = 0.140625f ;
        self.rootgrp.DX = 0.140625
        # DY = 0.09375f ;
        self.rootgrp.DY = 0.09375

    def finalize(self):
        self.rootgrp.close()


# main script
if __name__ == "__main__":
    stash_file = "./STASHmaster_A"
    # 20170710T0600Z
    start_date = dt.datetime(2017, 7, 11)
    end_date = dt.datetime(2019, 6, 25)
    td = dt.timedelta(days=1)
    output_dir = "OUTPUT_10KM"
    cur_date = start_date
    lis_example = nc4.Dataset("LIS_HIST_200711080000.d01.nc", "r")
    lat = lis_example["lat"][:, :]
    lon = lis_example["lon"][:, :]
    prev_smc_file = ""
    prev_snow_file = ""
    while cur_date <= end_date:
        out_dir = "%s/SURFACEMODEL/%04d%02d" % (
            output_dir, cur_date.year, cur_date.month)
        if os.path.isdir(out_dir) == False:
            os.mkdir(out_dir)
        # LIS_HIST_200711300000.d01.nc
        nc_file = "./%s/LIS_HIST_%04d%02d%02d0000.d01.nc" % (
            out_dir, cur_date.year, cur_date.month, cur_date.day)
        s2n = surf2lisnc(nc_file, stash_file)

        smc_file = "./WARM_START/%04d%02d%02dT0600Z_glu_smc" % (
            cur_date.year, cur_date.month, cur_date.day)
        # deal with missing data
        if os.path.isfile(smc_file) == False:
            smc_file = prev_smc_file
        s2n.read_surf_smc(smc_file)

        snow_file = "./WARM_START/%04d%02d%02dT0600Z_glu_snow" % (
            cur_date.year, cur_date.month, cur_date.day)
        if os.path.isfile(snow_file) == False:
            snow_file = prev_snow_file
        s2n.read_surf_swe(snow_file)
        # lat lon
        s2n.write_lat_lon(lat, lon)
        s2n.write_attributes()
        s2n.finalize()
        cur_date = cur_date + td
        prev_smc_file = smc_file
        prev_snow_file = snow_file
