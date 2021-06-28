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
#
#------------------------------------------------------------------------------
"""

# Standard modules
import sys

# Third party modules
import netCDF4 as nc4
from osgeo import gdal, osr

def _usage():
    """Print command line usage."""
    print("[INFO] Usage: %s ldtfile lvtfile anomaly_gt_prefix climo_gt_prefix"
          % (sys.argv[0]))
    print("[INFO]   where:")
    print("[INFO]    ldtfile: LDT parameter file with full lat/lon data")
    print("[INFO]    lvtfile: LVT 'TS' soil moisture anomaly file")
    print("[INFO]    anomaly_gt_prefix: prefix for new anomaly GeoTIFF file")
    print("[INFO]    climo_gt_prefix: prefix for new climatology GeoTIFF file")

def _read_cmd_args():
    """Read command line arguments."""
    # Check if argument count is correct
    if len(sys.argv) != 5:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if LDT parameter file can be opened.
    # NOTE:  Pylint complains about netCDF4 not having Dataset as a
    # member.  But this is false -- Dataset is a class, and we use
    # the constructor below.  Since this appears to be a bug in pylint, we
    # disable the no-member check in this section.
    # pylint: disable=no-member
    _ldtfile = sys.argv[1]
    ncid_ldt = nc4.Dataset(_ldtfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_ldt.close()

    # Check of LVT anomaly file can be opened
    _lvtfile = sys.argv[2]
    ncid_lvt = nc4.Dataset(_lvtfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_lvt.close()
    # pylint: enable=no-member

    _outfile_anomaly_prefix = sys.argv[3]
    _outfile_climo_prefix = sys.argv[4]

    return _ldtfile, _lvtfile, _outfile_anomaly_prefix, _outfile_climo_prefix

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
    print(_geotransform)
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

# Main driver
if __name__ == "__main__":

    # Get the file names for this invocation.
    ldtfile, lvtfile, outfile_anomaly_prefix, outfile_climo_prefix = \
        _read_cmd_args()

    # Fetch data from LVT output file
    # NOTE:  Pylint complains about netCDF4 not having Dataset as a
    # member.  But this is false -- Dataset is a class, and we use
    # the constructor below.  Since this appears to be a bug in pylint, we
    # disable the no-member check in this function.
    # pylint: disable=no-member
    ncid = nc4.Dataset(ldtfile, 'r', format='NETCDF4')
    longitudes = ncid.variables["lon"][:,:]
    latitudes = ncid.variables["lat"][:,:]
    ncid.close()

    for i in range(0, 4): # Loop across four LSM layers
        ncid = nc4.Dataset(lvtfile, 'r', format='NETCDF4')
        sm_anomalies = ncid.variables["SoilMoist"][i,:,:]
        nrows, ncols = sm_anomalies.shape
        # pylint: enable=no-member

        # Write soil moisture anomalies to GeoTIFF
        sm1 = sm_anomalies[::-1, :]
        geotransform = _make_geotransform(longitudes, latitudes, ncols, nrows)
        outfile_anomaly = "%s.layer%d.tif" %(outfile_anomaly_prefix,
                                             i)
        output_raster = _create_output_raster(outfile_anomaly,
                                              ncols, nrows, geotransform,
                                              sm1)
        output_raster.FlushCache() # Write to disk
        del output_raster

