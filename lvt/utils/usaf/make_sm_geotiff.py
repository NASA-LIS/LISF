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
# 09 July 2021: Eric Kemp (SSAI), added metadata to GeoTIFF files describing
#               raster fields.  Also changed numbering of soil layers from
#               0-3 to 1-4.
# 06 Dec 2022:  Eric Kemp (SSAI), updates to improve pylint score.
# 16 May 2023:  Eric Kemp (SSAI), updates to use 557 WW file convention for
#               output.
# 02 Jun 2023:  Eric Kemp (SSAI), further updates for 557 WW file
#               convention for output.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import sys

# Third party modules
# NOTE: pylint cannot see the Dataset class in netCDF4 since the latter is
# not written in Python.  We therefore disable a check for this line to
# avoid a known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
from osgeo import gdal, osr

_MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
           "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

_SOIL_LAYERS = {
    "NOAH" :    ["0-0.1 m", "0.1-0.4 m",  "0.4-1.0 m",  "1.0-2.0 m"],
    "NOAHMP" : ["0-0.1 m", "0.1-0.4 m",  "0.4-1.0 m",  "1.0-2.0 m"],
    "JULES" :   ["0-0.1 m", "0.1-0.35 m", "0.35-1.0 m", "1.0-3.0 m"],
}

_557WW_SOIL_LAYERS = {
    "NOAH" : ["D0CM-D10CM", "D10CM-D40CM", "D40CM-D100CM",
              "D100CM-D200CM"],
    "NOAHMP" : ["D0CM-D10CM", "D10CM-D40CM", "D40CM-D100CM",
                "D100CM-D200CM"],
    "JULES" : ["D0CM-D10CM", "D10CM-D35CM", "D35CM-D100CM",
                "D100CM-D300CM"],
}

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} ldtfile tsfile finalfile"
    txt += " LSM yyyymmddhh"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]   ldtfile: LDT parameter file with full lat/lon data")
    print("[INFO]   tsfile: LVT 'TS' soil moisture anomaly file")
    print("[INFO]   finalfile: LVT 'FINAL' soil moisture anomaly file")
    print("[INFO]   LSM: land surface model")
    print("[INFO]   yyyymmddhh: Valid date and time (UTC)")

def _read_cmd_args():
    """Read command line arguments."""
    # Check if argument count is correct
    if len(sys.argv) != 6:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if LDT parameter file can be opened.
    ldtfile = sys.argv[1]
    ncid_ldt = nc4_dataset(ldtfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_ldt.close()

    # Check of LVT TS anomaly file can be opened
    tsfile = sys.argv[2]
    ncid_lvt = nc4_dataset(tsfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_lvt.close()

    # Check of LVT FINAL anomaly file can be opened
    finalfile = sys.argv[3]
    ncid_lvt = nc4_dataset(finalfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_lvt.close()

    lsm = sys.argv[4]
    yyyymmddhh = sys.argv[5]

    if lsm not in ["NOAH", "NOAHMP", "JULES"]:
        print(f"[ERR] Unknown LSM {lsm}")
        print("Options are NOAH, NOAHMP, JULES")
        sys.exit(1)

    cmd_args = {
        "ldtfile" : ldtfile,
        "tsfile" : tsfile,
        "finalfile" : finalfile,
        "lsm" : lsm,
        "yyyymmddhh" : yyyymmddhh,
    }
    return cmd_args

def _make_geotransform(lon, lat, nxx, nyy):
    """Set affine transformation from image coordinate space to georeferenced
    space.  See https://gdal.org/tutorials/geotransforms_tut.html"""
    xmin = lon.min()
    xmax = lon.max()
    ymin = lat.min()
    ymax = lat.max()
    xres = (xmax - xmin) / float(nxx)
    yres = (ymax - ymin) / float(nyy)
    # Based on gdal.org/tutorials/geotransforms_tut.html
    # xmin is x-coordinate of upper-left corner of upper-left pixel
    # ymax is y-coordinate of upper-left corner of upper-left pixel
    # Third variable is row rotation, set to zero for north-up image
    # Fourth variable is column rotation, set to zero
    # Sixth variable is n-s pixel resolution (negative for north-up image)
    geotransform = (xmin, xres, 0, ymax, 0, -1*yres)
    return geotransform

def _create_output_raster(outfile, nxx, nyy, geotransform, var1):
    """Create the output raster file (the GeoTIFF), including map projection"""
    output_raster = gdal.GetDriverByName('GTiff').Create(outfile,
                                                         nxx, nyy, 1,
                                                         gdal.GDT_Float32)

    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326) # Corresponds to WGS 84
    output_raster.GetRasterBand(1).SetNoDataValue(-9999)
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(var1)

    return output_raster

def _set_metadata(varname, soil_layer, model, \
                  yyyymmddhh, \
                  climomonth=None):
    """Create metadata dictionary for output to GeoTIFF file"""
    validdt = datetime.datetime(year=int(yyyymmddhh[0:4]),
                                month=int(yyyymmddhh[4:6]),
                                day=int(yyyymmddhh[6:8]),
                                hour=int(yyyymmddhh[8:10]))
    metadata = { 'varname' : f'{varname}',
                 'units' : 'm3/m3',
                 'soil_layer' : f'{soil_layer}',
                 'land_surface_model' : f'{model}' }
    if climomonth is None:
        time_string = f"Valid {validdt.hour:02}Z {validdt.day} "
        time_string += f"{_MONTHS[validdt.month-1]} {validdt.year:04}"
        metadata["valid_date_time"] = time_string
    else:
        time_string = f"Updated {validdt.hour:02}Z {validdt.day} "
        time_string += f"{_MONTHS[validdt.month-1]} {validdt.year:04}"
        metadata["update_date_time"] = time_string
        metadata["climo_month"] = climomonth

    return metadata

def _make_outfile_anomaly(lsm, i, yyyymmddhh):
    """Create anomaly filename"""
    filename = "PS.557WW_SC.U_DI.C"
    filename += f"_GP.LIS-{lsm}"
    filename += "_GR.C0P09DEG_AR.GLOBAL"
    filename += f"_LY.{_557WW_SOIL_LAYERS[lsm][i]}"
    filename += f"_PA.SM-ANOMALY"
    filename += f"_DD.{yyyymmddhh[0:8]}"
    filename += f"_DT.{yyyymmddhh[8:10]}00_DF.TIF"
    return filename

def _proc_sm_anomalies(cmd_args, longitudes, latitudes):
    """Process soil moisture anomalies"""
    # Next, fetch the soil moisture anomalies from the LVT 'TS' file.
    ncid = nc4_dataset(cmd_args["tsfile"], 'r', format='NETCDF4')
    for i in range(0, 4): # Loop across four LSM layers
        sm_anomalies = ncid.variables["SoilMoist"][i,:,:]
        nrows, ncols = sm_anomalies.shape

        soil_layer = _SOIL_LAYERS[cmd_args["lsm"]][i]

        # Write soil moisture anomalies to GeoTIFF
        sm1 = sm_anomalies[::-1, :]
        geotransform = _make_geotransform(longitudes, latitudes, ncols, nrows)
        outfile_anomaly = \
            _make_outfile_anomaly(cmd_args["lsm"], i,
                                  cmd_args["yyyymmddhh"])
        varname = "Soil Moisture Anomaly"
        output_raster = _create_output_raster(outfile_anomaly,
                                              ncols, nrows, geotransform,
                                              sm1)
        metadata = _set_metadata(varname=varname,
                                 soil_layer=soil_layer,
                                 model=cmd_args["lsm"],
                                 yyyymmddhh=cmd_args["yyyymmddhh"])
        output_raster.GetRasterBand(1).SetMetadata(metadata)
        output_raster.FlushCache() # Write to disk
        del output_raster
    ncid.close()

def _make_outfile_climo(lsm, i, month, yyyymmddhh):
    """Create climatology filename"""
    filename = "PS.557WW_SC.U_DI.C_DC.CLIMO"
    filename += f"_GP.LIS-{lsm}"
    filename += "_GR.C0P09DEG_AR.GLOBAL"
    filename += f"_LY.{_557WW_SOIL_LAYERS[lsm][i]}"
    filename += f"_PA.SM-{month}"
    filename += f"_DP.20080101-{yyyymmddhh[0:8]}"
    filename += f"_DF.TIF"
    return filename

def _proc_sm_climo(cmd_args, longitudes, latitudes):
    """Process soil moisture climatology data"""
    ncid = nc4_dataset(cmd_args["finalfile"], 'r', format='NETCDF4')
    for imonth in range(0, 12):
        month = _MONTHS[imonth]
        climo_name = f"SoilMoist_{month}_climo"
        for i in range(0, 4): # Loop across four LSM layers
            nrows, ncols = ncid.variables[climo_name][i,:,:].shape

            # Write soil moisture climatology to GeoTIFF
            geotransform = _make_geotransform(longitudes, latitudes,
                                              ncols, nrows)
            outfile_climo = \
                _make_outfile_climo(cmd_args["lsm"], i,
                                    month, cmd_args["yyyymmddhh"])
            varname = "Climatological Soil Moisture"
            output_raster = \
                _create_output_raster(outfile_climo,
                                      ncols, nrows, geotransform,
                                      ncid.variables[climo_name][i,::-1,:])
            metadata = \
                _set_metadata(varname=varname,
                              soil_layer=_SOIL_LAYERS[cmd_args["lsm"]][i],
                              model=cmd_args["lsm"],
                              yyyymmddhh=cmd_args["yyyymmddhh"],
                              climomonth=month)
            output_raster.GetRasterBand(1).SetMetadata(metadata)

            output_raster.FlushCache() # Write to disk
            del output_raster
    ncid.close()

def _main():
    """Main driver"""
    # Get the file names for this invocation.
    cmd_args = _read_cmd_args()

    # First, fetch latitude/longitudes.  This is pulled from the LDT parameter
    # file, since LVT output has data voids over water.
    ncid = nc4_dataset(cmd_args["ldtfile"], 'r', format='NETCDF4')
    longitudes = ncid.variables["lon"][:,:]
    latitudes = ncid.variables["lat"][:,:]
    ncid.close()

    # Next, fetch the soil moisture anomalies from the LVT 'TS' file.
    _proc_sm_anomalies(cmd_args, longitudes, latitudes)

    # Next, fetch the monthly soil moisture climatologies from the LVT 'FINAL'
    # file.
    _proc_sm_climo(cmd_args, longitudes, latitudes)

# Main driver
if __name__ == "__main__":
    _main()
