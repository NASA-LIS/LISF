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
#------------------------------------------------------------------------------
#
# SCRIPT: make_s2s_metric_geotiff.py
#
# PURPOSE: Generates GeoTIFF files of medians of S2S metrics.  Only a single
# metric is processed per invocation of this script.
#
# REQUIREMENTS as of 28 May 2023:
# * Python 3.9 or higher
# * NumPy Python library 1.23.4 or higher
# * UNIDATA NetCDF4 Python library
# * GDAL Python library (bundled in osgeo package)
#
# REVISION HISTORY:
# * 18 Oct 2021: Eric Kemp/SSAI, first version.
# * 28 Oct 2021: Eric Kemp/SSAI, fixed pylint string complaints.  Changed logic
#   to calculate single median across entire ensemble collection, per metric
#   per month.
# * 30 Oct 2021: Eric Kemp/SSAI, added support for s2smetric config.
# * 02 Jun 2023: K. Arsenault + S. Mahanama, Updated the 557 WW file names.
#
#------------------------------------------------------------------------------
"""


# Standard modules
import os
import sys
import datetime
import glob
import yaml

# Third party modules
# NOTE: pylint cannot see the Dataset class in netCDF4 since that latter is not
# written in Python. We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
import numpy as np
from osgeo import gdal, osr

# Private class
class _MetricGeoTiff:
    """Class for building GeoTIFF files from medians of S2S metrics."""

    def __init__(self, topdir, metric, config):
        """Constructor"""
        self.topdir = topdir
        self.metric = metric
        self.config = config
        self.nmme_metric_files = {}
        self.filename_elements = {}
        self.median_data = {}

    def save_nmme_metric_filenames(self):
        """Searches for S2S metric files for each NMME model."""
        nmme_models = self.config["EXP"]["NMME_models"]
        for nmme in nmme_models:
            #subdir = "%s/%s" %(self.topdir, nmme)
            subdir = f"{self.topdir}"
            if not os.path.exists(subdir):
                txt = \
                    f"[WARN] Cannot find directory {subdir} for S2S metrics!"
                print(txt)
                continue
            regex = f"{subdir}/PS.*GP.LIS-S2S-{nmme.upper()}*DF.NC"
            files = glob.glob(regex)
            if len(files) == 0:
                print(f"[WARN] Cannot find metric file in {subdir}")
                continue
            if len(files) > 1:
                print(f"[WARN] Too many metric files in {subdir}, skipping...")
                continue
            self.nmme_metric_files[nmme] = files[0]

    def save_filename_elements(self):
        """Save elements from first NMME filename."""
        nmme = next(iter(self.nmme_metric_files))
        metric_file = self.nmme_metric_files[nmme]
        metric_file = os.path.basename(metric_file)
        filename_elements = {}
        element_list = metric_file.split("_")
        for element in element_list:
            key = element.split(".")[0]
            value = element.split(".")[1]
            filename_elements[key] = value
        self.filename_elements = filename_elements

    def calc_medians(self):
        """Calculate medians of the metrics."""

        self.median_data = {}
        total_ens_size = 0
        lead = 0
        latitude = 0
        longitude = 0

        # We need to first count the total number of ensembles.
        nmme_models = self.config["EXP"]["NMME_models"]
        for nmme in nmme_models:
            if nmme not in self.nmme_metric_files:
                continue
            metric_file = self.nmme_metric_files[nmme]
            rootgrp = nc4_dataset(metric_file, 'r',
                                  format="NETCDF4_CLASSIC")
            total_ens_size += rootgrp.dimensions["ens"].size
            lead = rootgrp.dimensions["time"].size
            latitude = rootgrp.dimensions["latitude"].size
            longitude = rootgrp.dimensions["longitude"].size

        # Find dimensions for a single month
        dims = (total_ens_size, latitude, longitude)

        # Now loop through each month, load a NMME contribution into the
        # total var array, and calculate the median.
        self.median_data["median"] = []
        for itime in range(0, lead):
            var = np.zeros(dims)
            iens = 0
            start_end_set = False
            for nmme in nmme_models:
                if nmme not in self.nmme_metric_files:
                    continue
                metric_file = self.nmme_metric_files[nmme]
                rootgrp = nc4_dataset(metric_file, 'r',
                                      format="NETCDF4_CLASSIC")
                ens = rootgrp.dimensions["ens"].size
                var[iens:(iens+ens), :, :] = \
                    rootgrp.variables[self.metric][:, itime, :, :].data
                iens += ens
                if not start_end_set:
                    self.median_data["num_months"] = \
                        rootgrp.dimensions["time"].size
                    self.median_data["latitudes"] = \
                        rootgrp.variables["latitude"][:]
                    self.median_data["longitudes"] = \
                        rootgrp.variables["longitude"][:]
                    self.median_data["units"] = \
                        rootgrp.variables[self.metric].units
                    self.median_data["long_name"] = \
                        rootgrp.variables[self.metric].long_name
                    self.set_startdates_enddates_by_month()
                    start_end_set = True
                rootgrp.close()

            # Calculate the median for the current month
            self.median_data["median"].append(np.median(var, axis=0))
            del var

    def set_startdates_enddates_by_month(self):
        """Set the startdates and enddates by month."""
        startdates_by_month = []
        enddates_by_month = []
        num_months = self.get_num_of_months()
        filename_elements = self.get_filename_elements()
        for imonth in range(0, num_months):
            if imonth == 0:
                first_startdate = \
                    _get_first_startdate(filename_elements)
                startdates_by_month.append(first_startdate)
                newdate2 = _set_newdate(first_startdate)
                enddates_by_month.append(newdate2)
            else:
                newdate1 = _set_newdate(startdates_by_month[imonth-1])
                newdate2 = _set_newdate(enddates_by_month[imonth-1])
                startdates_by_month.append(newdate1)
                enddates_by_month.append(newdate2)
        self.median_data["startdates_by_month"] = \
            startdates_by_month
        self.median_data["enddates_by_month"] = \
            enddates_by_month

    def get_geotransform(self):
        """Set affine transformation from image coordinate space to
        georeferenced space. See
        https://gdal.org/tutorials/geotransforms_tut.html"""
        xmin = self.median_data["longitudes"].min()
        xmax = self.median_data["longitudes"].max()
        ymin = self.median_data["latitudes"].min()
        ymax = self.median_data["latitudes"].max()
        nxx = self.median_data["longitudes"].size
        nyy = self.median_data["latitudes"].size
        xres = (xmax - xmin) / float(nxx-1)
        yres = (ymax - ymin) / float(nyy-1)
        xmin_corner = xmin - 0.5*xres
        ymax_corner = ymax + 0.5*yres
        # Below is based on gdal.org/tutorials/geotransforms_tut.html
        # xmin is x-coordinate of upper-left corner of upper-left pixel
        # ymax in y-coordinate of upper-left corner of upper-left pixel
        # Third variable is row rotation, set to zero for north-up image
        # Fifth variable is column rotation, set to zero
        # Sixth variable is n-s pixel resolution (negative for north-up image)
        #return (xmin, xres, 0, ymax, 0, -1*yres)
        return (xmin_corner, xres, 0, ymax_corner, 0, -1*yres)

    def make_geotiff_filename(self, metric, imonth):
        """Make name of new geotiff file."""
        startdate = \
            self.median_data["startdates_by_month"][imonth]
        enddate = self.median_data["enddates_by_month"][imonth]
        filename_elements = self.filename_elements
        filename = f"{self.topdir}"
        filename += f"/PS.{filename_elements['PS']}"
        filename += f"_SC.{filename_elements['SC']}"
        filename += f"_DI.{filename_elements['DI']}"
        metricname = metric.split("_")[1]
        filename += "_GP.LIS-S2S-"+metricname
        filename += f"_GR.{filename_elements['GR']}"
        filename += f"_AR.{filename_elements['AR']}"

        variable = metric.split("_")[0]
        filename += f"_PA.{variable.upper()}"
        filename += f"_DD.{filename_elements['DD']}"
        filename += \
            f"_FP.{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
        filename += \
            f"-{enddate.year:4d}{enddate.month:02d}{enddate.day:02d}"

        filename += "_DF.TIF"
        return filename

    def create_output_raster(self, outfile):
        """Create the output raster file (the GeoTIFF), including map
        projection"""
        nxx = self.median_data["longitudes"].size
        nyy = self.median_data["latitudes"].size
        geotransform = self.get_geotransform()
        options = ["COMPRESS=NONE"]
        output_raster = gdal.GetDriverByName('GTiff').Create(outfile,
                                                             nxx, nyy,
                                                             1,
                                                             gdal.GDT_Float64,
                                                             options=options)
        output_raster.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326) # Corresponds to WGS 84
        output_raster.SetProjection(srs.ExportToWkt())
        return output_raster

    def get_num_of_months(self):
        """Get the number of months for a given metric."""
        return self.median_data["num_months"]

    def get_units(self):
        """Get the units for a given metric."""
        return self.median_data["units"]

    def get_long_name(self):
        """Get the long_name for a given metric."""
        return self.median_data["long_name"]

    def get_startdates_by_month(self):
        """Get the startdates for each month for a given metric."""
        return self.median_data["startdates_by_month"]

    def get_median_values(self):
        """Get the metric values for a given metric."""
        return self.median_data["median"]

    def get_filename_elements(self):
        """Get the filename elements for a given metric."""
        return self.filename_elements

# Private module methods
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} topdir metric configfile"
    print(txt)
    print("[INFO] where:")
    print("[INFO] topdir: top directory with all S2S metrics netCDF files.")
    print("[INFO] metric: name of metric to process.")
    print("[INFO] configfile: path to s2smetric config file.")

def _read_cmd_args():
    """Read command line arguments."""
    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)
    topdir = sys.argv[1]
    if not os.path.exists(topdir):
        print(f"[ERR] {topdir} does not exist!")
        sys.exit(1)
    metric = sys.argv[2]
    configfile = sys.argv[3]
    if not os.path.exists(configfile):
        print(f"[ERR] Config file {configfile} not found!")
        sys.exit(1)
    return topdir, metric, configfile

def _find_nmme_in_filename(config, filename):
    """Find and return NMME model in filename."""
    basename = os.path.basename(filename)
    nmme_models = config["EXP"]["NMME_models"]
    for nmme_model in nmme_models:
        upper = nmme_model.upper()
        if upper in basename:
            return upper
    return None

def _get_first_startdate(filename_elements):
    """Get start date from name of metric file"""
    try:
        yyyymmdd = filename_elements["FP"].split("_")[0]
    except KeyError:
        print("[ERR] Cannot resolve data period from NMME files!")
        sys.exit(1)
    return datetime.datetime(year=int(yyyymmdd[0:4]),
                             month=int(yyyymmdd[4:6]),
                             day=int(yyyymmdd[6:8]))

def _set_newdate(date):
    """Set a new date one month in the future."""
    if date.month == 12:
        newdate = datetime.datetime(year=(date.year + 1),
                                    month=1,
                                    day=1)
    else:
        newdate = datetime.datetime(year=date.year,
                                    month=(date.month + 1),
                                    day=1)
    return newdate

def _driver():
    """Main driver."""
    topdir, metric, configfile = _read_cmd_args()
    with open(configfile, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    mgt = _MetricGeoTiff(topdir, metric, config)
    mgt.save_nmme_metric_filenames()
    mgt.save_filename_elements()
    mgt.calc_medians()
    metadata = {}

    # Write single GeoTIFF file for each month.
    # Only one metric is processed.
    num_months = mgt.get_num_of_months()
    metadata["varname"] = metric
    metadata["units"] = mgt.get_units()
    metadata["long_name"] = mgt.get_long_name()
    metadata["generating_process"] = "LIS-S2S-ANOM"
    startdates_by_month = mgt.get_startdates_by_month()
    median_values = mgt.get_median_values()
    for imonth in range(0, num_months):
        metadata["forecast_month"] = f"{imonth + 1}"
        metadata["valid_year_and_month"] = \
            f"{startdates_by_month[imonth].year:04d}" + \
            f"-{startdates_by_month[imonth].month:02d}"
        var2d = median_values[imonth][::-1, :] # Flip y-axis
        geotiff_filename = \
            mgt.make_geotiff_filename(metric, imonth)
        output_raster = mgt.create_output_raster(geotiff_filename)
        output_raster.SetMetadata(metadata)
        output_raster.GetRasterBand(1).SetNoDataValue(-9999)
        output_raster.GetRasterBand(1).WriteArray(var2d)
        output_raster.GetRasterBand(1).SetMetadata(metadata)
        output_raster.FlushCache() # Write to disk
        del output_raster

# Invoke driver
if __name__ == "__main__":
    _driver()
