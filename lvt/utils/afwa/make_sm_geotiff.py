#!/usr/bin/env python3

"""
#------------------------------------------------------------------------------
#
# SCRIPT: make_sm_geotiff.py
#
# PURPOSE:  Extracts soil moisture climatology, soil moisture anomaly, and
# latitude/longitude information from netCDF files, and converts to GeoTIFF.
#
# REQUIREMENTS as of 25 June 2021:
# * Python 3.8 or higher.
# * UNIDATA NetCDF4 Python library (for reading netCDF4 files)
# * GDAL Python library (bundled in osgeo package).
#
# REVISION HISTORY:
# 25 June 2021: Eric Kemp (SSAI), first version, based on code provided by
#               Sujay Kumar (NASA GSFC).
# 28 June 2021: Eric Kemp (SSAI), add support for each soil layer.
# 29 June 2021: Eric Kemp (SSAI), added processing for climatologies for all
#               months.  Also tweaked pylint checking for netCDF4 module.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import sys

# Third party modules
# NOTE: pylint cannot see the Dataset class in netCDF4 since the latter is
# not written in Python.  We therefore disable a check for this line to
# avoid a known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
from osgeo import gdal, osr

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s ldtfile tsfile finalfile" %(sys.argv[0])
    txt += " anomaly_gt_prefix climo_gt_prefix"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]   ldtfile: LDT parameter file with full lat/lon data")
    print("[INFO]   tsfile: LVT 'TS' soil moisture anomaly file")
    print("[INFO]   finalfile: LVT 'FINAL' soil moisture anomaly file")
    print("[INFO]   anomaly_gt_prefix: prefix for new anomaly GeoTIFF files")
    print("[INFO]   climo_gt_prefix: prefix for new climatology GeoTIFF files")

def _read_cmd_args():
    """Read command line arguments."""
    # Check if argument count is correct
    if len(sys.argv) != 6:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if LDT parameter file can be opened.
    _ldtfile = sys.argv[1]
    ncid_ldt = nc4_dataset(_ldtfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_ldt.close()

    # Check of LVT TS anomaly file can be opened
    _tsfile = sys.argv[2]
    ncid_lvt = nc4_dataset(_tsfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_lvt.close()

    # Check of LVT FINAL anomaly file can be opened
    _finalfile = sys.argv[3]
    ncid_lvt = nc4_dataset(_finalfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_lvt.close()

    _outfile_anomaly_prefix = sys.argv[4]
    _outfile_climo_prefix = sys.argv[5]

    return _ldtfile, _tsfile, _finalfile, \
        _outfile_anomaly_prefix, _outfile_climo_prefix

def _make_geotransform(lon, lat, nxx, nyy):
    """Set affine transformation from image coordinate space to georeferenced
    space.  See https://gdal.org/tutorials/geotransforms_tut.html"""
    xmin = lon.min()
    xmax = lon.max()
    ymin = lat.min()
    ymax = lat.max()
    xres = (xmax - xmin) / float(nxx)
    yres = (ymax - ymin) / float(nyy)
    # Sujay's original code
    #geotransform = (xmax, xres, 0, ymin, 0, -yres)
    # Eric's code...Based on gdal.org/tutorials/geotransforms_tut.html
    # xmin is x-coordinate of upper-left corner of upper-left pixel
    # ymax is y-coordinate of upper-left corner of upper-left pixel
    # Third variable is row rotation, set to zero for north-up image
    # Fourth variable is column rotation, set to zero
    # Sixth variable is n-s pixel resolution (negative for north-up image)
    _geotransform = (xmin, xres, 0, ymax, 0, -1*yres)
    #print(_geotransform)
    return _geotransform

def _create_output_raster(outfile, nxx, nyy, _geotransform, var1):
    """Create the output raster file (the GeoTIFF), including map projection"""
    _output_raster = gdal.GetDriverByName('GTiff').Create(outfile,
                                                          nxx, nyy, 1,
                                                          gdal.GDT_Float32)

    _output_raster.SetGeoTransform(_geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326) # Corresponds to WGS 84
    _output_raster.GetRasterBand(1).SetNoDataValue(-9999)
    _output_raster.SetProjection(srs.ExportToWkt())
    _output_raster.GetRasterBand(1).WriteArray(var1)
    return _output_raster

_MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
           "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

# Main driver
if __name__ == "__main__":

    # Get the file names for this invocation.
    ldtfile, tsfile, finalfile, anomaly_gt_prefix, climo_gt_prefix = \
        _read_cmd_args()

    # First, fetch latitude/longitudes.  This is pulled from the LDT parameter
    # file, since LVT output has data voids over water.
    ncid = nc4_dataset(ldtfile, 'r', format='NETCDF4')
    longitudes = ncid.variables["lon"][:,:]
    latitudes = ncid.variables["lat"][:,:]
    ncid.close()

    # Next, fetch the soil moisture anomalies from the LVT 'TS' file.
    ncid = nc4_dataset(tsfile, 'r', format='NETCDF4')
    for i in range(0, 4): # Loop across four LSM layers
        sm_anomalies = ncid.variables["SoilMoist"][i,:,:]
        nrows, ncols = sm_anomalies.shape

        # Write soil moisture anomalies to GeoTIFF
        sm1 = sm_anomalies[::-1, :]
        geotransform = _make_geotransform(longitudes, latitudes, ncols, nrows)
        outfile_anomaly = "%s.layer%d.tif" %(anomaly_gt_prefix,
                                             i)
        output_raster = _create_output_raster(outfile_anomaly,
                                              ncols, nrows, geotransform,
                                              sm1)
        output_raster.FlushCache() # Write to disk
        del output_raster
    ncid.close()

    # Next, fetch the monthly soil moisture climatologies from the LVT 'FINAL'
    # file.
    ncid = nc4_dataset(finalfile, 'r', format='NETCDF4')
    for imonth in range(0, 12):
        month = _MONTHS[imonth]
        climo_name = "SoilMoist_%s_climo" %(month)
        for i in range(0, 4): # Loop across four LSM layers
            sm_climo = ncid.variables[climo_name][i,:,:]
            nrows, ncols = sm_climo.shape

            # Write soil moisture climatology to GeoTIFF
            sm1 = sm_climo[::-1, :]
            geotransform = _make_geotransform(longitudes, latitudes,
                                              ncols, nrows)
            outfile_climo = "%s.%s.layer%d.tif" %(climo_gt_prefix,
                                                  month,
                                                  i)
            output_raster = _create_output_raster(outfile_climo,
                                                  ncols, nrows, geotransform,
                                                  sm1)
            output_raster.FlushCache() # Write to disk
            del output_raster
    ncid.close()
