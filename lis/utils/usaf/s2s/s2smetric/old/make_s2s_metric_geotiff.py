#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
#------------------------------------------------------------------------------
#
# SCRIPT: make_s2s_metric_geotiff.py
#
# PURPOSE: Generates GeoTIFF file from S2S metrics netCDF file.
#
# REQUIREMENTS as of 4 Oct 2021:
# * Python 3.8 or higher
# * UNIDATA NetCDF4 Python library
# * GDAL Python library (bundled in osgeo package)
#
# REVISION HISTORY:
# * 4 Oct 2021: Eric Kemp (SSAI), first version.
# * 7 Oct 2021: Eric Kemp (SSAI), second attempt.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import sys

# Third party modules
# NOTE: pylint cannot see the Dataset class in netCDF4 since that latter is not
# written in Python. We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
from osgeo import gdal, osr

# Private constants
_VARNAMES = ["RootZone_SM_ANOM", "RootZone_SM_SANOM",
             "Streamflow_ANOM", "Streamflow_SANOM",
             "Surface_SM_ANOM", "Surface_SM_SANOM"]

# Private class
class _MetricGeoTiff:
    """Class for building GeoTIFFs from S2S Metric netCDF file."""

    def __init__(self, metricfile):
        """Constructor"""
        self.metricfile = metricfile

        # Parse the name of the metric file, and save elements.  Assumes
        # filename obeys Air Force Weather file naming convention.
        self.metric_filename_elements = {}
        element_list = self.metricfile.split("_")
        for element in element_list:
            key = element.split(".")[0]
            value = element.split(".")[1]
            self.metric_filename_elements[key] = value

        # Find number of months in metrics file, and save latitudes and
        # longitudes.
        rootgrp = nc4_dataset(self.metricfile, 'r', format="NETCDF4_CLASSIC")
        num_months = rootgrp.dimensions["lead"].size
        self.latitudes = rootgrp.variables["latitude"][:]
        self.longitudes = rootgrp.variables["longitude"][:]
        rootgrp.close()

        # Set start and end of each month of data in the metric file.
        self.startdates_by_month = []
        self.enddates_by_month = []
        for imonth in range(0, num_months):
            if imonth == 0:
                self.startdates_by_month.append(self.get_first_startdate())
                newdate2 = _set_newdate(self.get_first_startdate())
                self.enddates_by_month.append(newdate2)
            else:
                newdate1 = _set_newdate(self.startdates_by_month[imonth-1])
                newdate2 = _set_newdate(self.enddates_by_month[imonth-1])
                self.startdates_by_month.append(newdate1)
                self.enddates_by_month.append(newdate2)

    def get_first_startdate(self):
        """Get start date from name of metric file."""
        try:
            yyyymmdd = self.metric_filename_elements["DP"].split("_")[0]
        except KeyError:
            print("[ERR] Cannot resolve data period from file name %s" \
                  %(self.metricfile))
            sys.exit(1)
        return datetime.datetime(year=int(yyyymmdd[0:4]),
                                 month=int(yyyymmdd[4:6]),
                                 day=int(yyyymmdd[6:8]))

    def get_geotransform(self):
        """Set affine transformation from image coordinate space to
        georeferenced space. See
        https://gdal.org/tutorials/geotransforms_tut.html"""
        xmin = self.longitudes.min()
        xmax = self.longitudes.max()
        ymin = self.latitudes.min()
        ymax = self.latitudes.max()
        nxx = self.longitudes.size
        nyy = self.latitudes.size
        xres = (xmax - xmin) / float(nxx-1)
        yres = (ymax - ymin) / float(nyy-1)
        # Below is based on gdal.org/tutorials/geotransforms_tut.html
        # xmin is x-coordinate of upper-left corner of upper-left pixel
        # ymax in y-coordinate of upper-left corner of upper-left pixel
        # Third variable is row rotation, set to zero for north-up image
        # Fourth variable is column rotation, set to zero
        # Sixth variable is n-s pixel resolution (negative for north-up image)
        return (xmin, xres, 0, ymax, 0, -1*yres)

    def make_geotiff_filename(self, varname, i_ens_member, imonth):
        """Make name of new geotiff file."""
        startdate = self.startdates_by_month[imonth]
        enddate = self.enddates_by_month[imonth]
        filename = "PS.%s" %(self.metric_filename_elements["PS"])
        filename += "_SC.%s" %(self.metric_filename_elements["SC"])
        filename += "_DI.%s" %(self.metric_filename_elements["DI"])
        filename += "_GP.%s-%2.2d" \
            %(self.metric_filename_elements["GP"], i_ens_member)
        filename += "_GR.%s" %(self.metric_filename_elements["GR"])
        filename += "_AR.%s" %(self.metric_filename_elements["AR"])
        filename += "_PA.LIS-S2S-%s" %(varname.replace("_","-").upper())
        filename += "_DP.%4.4d%2.2d%2.2d-%4.4d%2.2d%2.2d" %(startdate.year,
                                                            startdate.month,
                                                            startdate.day,
                                                            enddate.year,
                                                            enddate.month,
                                                            enddate.day)
        filename += "_TP.%s" %(self.metric_filename_elements["TP"])
        filename += "_DF.TIF"
        return filename

    def get_variable(self, varname):
        """Retrieve variable and select metadata from metric file."""
        rootgrp = nc4_dataset(self.metricfile, 'r', format="NETCDF4_CLASSIC")
        var = rootgrp.variables[varname][:,:,:,:]
        units = rootgrp.variables[varname].units
        long_name = rootgrp.variables[varname].long_name
        rootgrp.close()
        return var, units, long_name

    def get_num_of_months(self):
        """Fetch number of months in LIS output."""
        rootgrp = nc4_dataset(self.metricfile, 'r', format="NETCDF4_CLASSIC")
        num_months = rootgrp.dimensions["lead"].size
        rootgrp.close()
        return num_months

    def get_ensemble_size(self):
        """Fetch number of ensemble members in LIS output."""
        rootgrp = nc4_dataset(self.metricfile, 'r', format="NETCDF4_CLASSIC")
        num_ens = rootgrp.dimensions["ens"].size
        rootgrp.close()
        return num_ens

    def create_output_raster(self, outfile):
        """Create the output raster file (the GeoTIFF), including map
        projection"""
        nxx = self.longitudes.size
        nyy = self.latitudes.size
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

# Private module methods
def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s metricfile" %(sys.argv[0])
    print(txt)
    print("[INFO] where:")
    print("[INFO] metricfile: netcdf file with metric data")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct
    if len(sys.argv) != 2:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if metric file can be opened.
    metricfile = sys.argv[1]
    rootgrp = nc4_dataset(metricfile, mode="r", format="NETCDF4_CLASSIC")
    rootgrp.close()

    return metricfile

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
    metricfile = _read_cmd_args()
    mgt = _MetricGeoTiff(metricfile)
    metadata = {}
    metadata["generating_process"] = mgt.metric_filename_elements["GP"]
    num_months = mgt.get_num_of_months()
    ensemble_size = mgt.get_ensemble_size()
    # Write single GeoTIFF file for each month and ensemble member of each
    # metric.
    for varname in _VARNAMES:
        var, units, long_name = mgt.get_variable(varname)
        metadata["varname"] = varname
        metadata["units"] = units
        metadata["long_name"] = long_name
        for i_ens_member in range(0, ensemble_size):
            metadata["ensemble_member"] = "%s" %(i_ens_member)
            for imonth in range(0, num_months):
                metadata["forecast_month"] = "%s" %(imonth + 1)
                metadata["valid_year_and_month"] = "%4.4d-%2.2d" \
                    %(mgt.startdates_by_month[imonth].year,
                      mgt.startdates_by_month[imonth].month)
                var2d = var[i_ens_member, imonth, :, :]
                geotiff_filename = \
                    mgt.make_geotiff_filename(varname, i_ens_member, imonth)
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
