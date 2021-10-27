#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: make_s2s_metric_geotiff.py
#
# PURPOSE: Generates GeoTIFF files of medians of S2S metrics.  Only a single
# metric is processed per invocation of this script.
#
# REQUIREMENTS as of 18 Oct 2021:
# * Python 3.8 or higher
# * NumPy Python library 1.19 or higher
# * UNIDATA NetCDF4 Python library
# * GDAL Python library (bundled in osgeo package)
#
# REVISION HISTORY:
# * 18 Oct 2021: Eric Kemp/SSAI, first version.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import glob
import os
import sys

# Third party modules
# NOTE: pylint cannot see the Dataset class in netCDF4 since that latter is not
# written in Python. We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
import numpy as np
from osgeo import gdal, osr

# Private constants
#_NMME_MODELS = ["CCM4", "CCSM4", "CFSv2", "GEOSv2", "GFDL", "GNEMO"]
_NMME_MODELS = ["GEOSv2"]

_METRICS = ["RootZone_SM_ANOM", "RootZone_SM_SANOM",
            "Streamflow_ANOM",  "Streamflow_SANOM",
            "Surface_SM_ANOM",  "Surface_SM_SANOM"]

# Private class
class _MetricGeoTiff:
    """Class for building GeoTIFF files from medians of S2S metrics."""

    def __init__(self, topdir, metric):
        """Constructor"""
        self.topdir = topdir
        self.metric = metric
        self.nmme_metric_files = {}
        self.filename_elements = {}
        self.median_data = {}

    def save_nmme_metric_filenames(self):
        """Searches for S2S metric files for each NMME model."""
        for nmme in _NMME_MODELS:
            #subdir = "%s/%s" %(self.topdir, nmme)
            subdir = f"{self.topdir}"
            if not os.path.exists(subdir):
                txt = \
                    f"[WARN] Cannot find directory {subdir} for S2S metrics!"
                print(txt)
                continue
            regex = f"{subdir}/PS.*GP.LIS-S2S-{nmme.upper()}-ANOM*DF.NC"
            files = glob.glob(regex)
            if len(files) == 0:
                print(f"[WARN] Cannot find metric file in {subdir}")
                continue
            if len(files) > 1:
                print(f"[WARN] Too many metric files in {subdir}, skipping...")
                continue
            self.nmme_metric_files[nmme] = files[0]

    def save_filename_elements(self):
        """Save elements from filenames."""
        keys = self.nmme_metric_files.keys()
        for nmme in keys:
            metric_file = self.nmme_metric_files[nmme]
            metric_file = os.path.basename(metric_file)
            filename_elements = {}
            element_list = metric_file.split("_")
            for element in element_list:
                key = element.split(".")[0]
                value = element.split(".")[1]
                filename_elements[key] = value
            self.filename_elements[nmme] = filename_elements

    def calc_medians(self):
        """Calculate medians of the metrics."""
        self.median_data = {}
        metric = self.metric
        for nmme in _NMME_MODELS:
            if not nmme in self.nmme_metric_files:
                continue
            metric_file = self.nmme_metric_files[nmme]
            self.median_data[nmme] = {}
            rootgrp = nc4_dataset(metric_file, 'r',
                                  format="NETCDF4_CLASSIC")
            self.median_data[nmme]["num_months"] = \
                rootgrp.dimensions["lead"].size
            self.median_data[nmme]["latitudes"] = \
                rootgrp.variables["latitude"][:]
            self.median_data[nmme]["longitudes"] = \
                rootgrp.variables["longitude"][:]
            var = rootgrp.variables[metric][:,:,:,:].data
            self.median_data[nmme]["units"] = \
                rootgrp.variables[metric].units
            self.median_data[nmme]["long_name"] = \
                rootgrp.variables[metric].long_name
            rootgrp.close()
            self.median_data[nmme]["median"] = \
                np.median(var, axis=0)
            self.set_startdates_enddates_by_month(nmme)

    def set_startdates_enddates_by_month(self, nmme):
        """Set the startdates and enddates by month."""
        startdates_by_month = []
        enddates_by_month = []
        num_months = self.get_num_of_months(nmme)
        filename_elements = self.get_filename_elements(nmme)
        metricfile = self.nmme_metric_files[nmme]
        for imonth in range(0, num_months):
            if imonth == 0:
                first_startdate = \
                    _get_first_startdate(filename_elements, metricfile)
                startdates_by_month.append(first_startdate)
                newdate2 = _set_newdate(first_startdate)
                enddates_by_month.append(newdate2)
            else:
                newdate1 = _set_newdate(startdates_by_month[imonth-1])
                newdate2 = _set_newdate(enddates_by_month[imonth-1])
                startdates_by_month.append(newdate1)
                enddates_by_month.append(newdate2)
        self.median_data[nmme]["startdates_by_month"] = \
            startdates_by_month
        self.median_data[nmme]["enddates_by_month"] = \
            enddates_by_month

    def get_geotransform(self, nmme):
        """Set affine transformation from image coordinate space to
        georeferenced space. See
        https://gdal.org/tutorials/geotransforms_tut.html"""
        xmin = self.median_data[nmme]["longitudes"].min()
        xmax = self.median_data[nmme]["longitudes"].max()
        ymin = self.median_data[nmme]["latitudes"].min()
        ymax = self.median_data[nmme]["latitudes"].max()
        nxx = self.median_data[nmme]["longitudes"].size
        nyy = self.median_data[nmme]["latitudes"].size
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

    def make_geotiff_filename(self, metric, nmme, imonth):
        """Make name of new geotiff file."""
        startdate = \
            self.median_data[nmme]["startdates_by_month"][imonth]
        enddate = self.median_data[nmme]["enddates_by_month"][imonth]
        filename_elements = self.filename_elements[nmme]
        filename = f"{self.topdir}"
        filename += f"/PS.{filename_elements['PS']}"
        filename += f"_SC.{filename_elements['SC']}"
        filename += f"_DI.{filename_elements['DI']}"
        filename += f"_GP.{filename_elements['GP']}"
        filename += f"_GR.{filename_elements['GR']}"
        filename += f"_AR.{filename_elements['AR']}"
        filename += \
            f"_PA.LIS-S2S-{metric.replace('_','-').upper()}-ENS-MEDIAN"
        filename += \
            f"_DP.{startdate.year:4d}{startdate.month:2d}{startdate.day:2d}"
        filename += \
            f"-{enddate.year:4d}{enddate.month:2d}{enddate.day:2d}"
        filename += f"_TP.{filename_elements['TP']}"
        filename += "_DF.TIF"
        return filename

    def create_output_raster(self, nmme, outfile):
        """Create the output raster file (the GeoTIFF), including map
        projection"""
        nxx = self.median_data[nmme]["longitudes"].size
        nyy = self.median_data[nmme]["latitudes"].size
        geotransform = self.get_geotransform(nmme)
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

    def get_num_of_months(self, nmme):
        """Get the number of months for a given metric and NMME model."""
        return self.median_data[nmme]["num_months"]

    def get_units(self, nmme):
        """Get the units for a given metric and NMME model."""
        return self.median_data[nmme]["units"]

    def get_long_name(self, nmme):
        """Get the long_name for a given metric and NMME model."""
        return self.median_data[nmme]["long_name"]

    def get_generating_process(self, nmme):
        """Get the generating process for a given metric and NMME model."""
        return self.filename_elements[nmme]["GP"]

    def get_startdates_by_month(self, nmme):
        """Get the startdates for each month for a given metric and NMME
        model."""
        return self.median_data[nmme]["startdates_by_month"]

    def get_median_values(self, nmme):
        """Get the metric values for a given metric and NMME model."""
        return self.median_data[nmme]["median"]

    def get_filename_elements(self, nmme):
        """Get the filename elements for a given metric and NMME model."""
        return self.filename_elements[nmme]

# Private module methods
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} topdir metric"
    print(txt)
    print("[INFO] where:")
    print("[INFO] topdir: top directory with all S2S metrics netCDF files.")
    print("[INFO] metric: name of metric to process.")

def _read_cmd_args():
    """Read command line arguments."""
    if len(sys.argv) != 3:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)
    topdir = sys.argv[1]
    if not os.path.exists(topdir):
        print(f"[ERR] {topdir} does not exist!")
        sys.exit(1)
    metric = sys.argv[2]
    if metric not in _METRICS:
        print(f"[ERR] Unknown metric {metric} requested!")
        sys.exit(1)
    return topdir, metric

def _find_nmme_in_filename(filename):
    """Find and return NMME model in filename."""
    basename = os.path.basename(filename)
    for nmme_model in _NMME_MODELS:
        upper = nmme_model.upper()
        if upper in basename:
            return upper
    return None

def _get_first_startdate(filename_elements, metricfile):
    """Get start date from name of metric file"""
    try:
        yyyymmdd = filename_elements["DP"].split("_")[0]
    except KeyError:
        print(f"[ERR] Cannot resolve data period from file name {metricfile}")
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
    topdir, metric = _read_cmd_args()
    mgt = _MetricGeoTiff(topdir, metric)
    mgt.save_nmme_metric_filenames()
    mgt.save_filename_elements()
    mgt.calc_medians()
    metadata = {}
    # Write single GeoTIFF file for each month for each NMME model.
    # Only one metric is processed.
    for nmme in _NMME_MODELS:
        if nmme not in mgt.nmme_metric_files:
            continue
        num_months = mgt.get_num_of_months(nmme)
        metadata["varname"] = metric
        metadata["units"] = mgt.get_units(nmme)
        metadata["long_name"] = mgt.get_long_name(nmme)
        metadata["generating_process"] = \
            mgt.get_generating_process(nmme)
        startdates_by_month = \
            mgt.get_startdates_by_month(nmme)
        median_values = mgt.get_median_values(nmme)
        for imonth in range(0, num_months):
            metadata["forecast_month"] = f"{imonth + 1}"
            metadata["valid_year_and_month"] = \
                f"{startdates_by_month[imonth].year:4d}" + \
                f"-{startdates_by_month[imonth].month:2d}"
            var2d = median_values[imonth, ::-1, :] # Flip y-axis
            geotiff_filename = \
                mgt.make_geotiff_filename(metric, nmme, imonth)
            output_raster = mgt.create_output_raster(nmme,
                                                     geotiff_filename)
            output_raster.SetMetadata(metadata)
            output_raster.GetRasterBand(1).SetNoDataValue(-9999)
            output_raster.GetRasterBand(1).WriteArray(var2d)
            output_raster.GetRasterBand(1).SetMetadata(metadata)
            output_raster.FlushCache() # Write to disk
            del output_raster

# Invoke driver
if __name__ == "__main__":
    _driver()
