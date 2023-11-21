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
# SCRIPT: daily_s2spost_nc.py
#
# PURPOSE:  Merges daily netCDF output from LIS-NoahMP and LIS-HYMAP2 into
# single, CF-compliant netCDF4 file for distribution.  Rewritten to use
# NetCDF4 Python library instead of NCO software, to reduce runtime.
#
# REQUIREMENTS as of 02 Jun 2023:
# * Python 3.9 or higher.
# * UNIDATA NetCDF4 Python library
# * Numpy Python library
#
# REFERENCES:
# https://cfconventions.org for specifications of NetCDF Climate and Forecast
#   (CF) Metadata Conventions.
# https://unidata.ucar.edu/software/udunits for documentation on UDUNITS2
#   library, which CF is generally consistent with for unit specifications.
#
# REVISION HISTORY:
# 09 Sep 2022: Eric Kemp/SSAI, first version.
# 14 Nov 2022: K. Arsenault/SAIC, removed fields for FOC.
# 02 Jun 2023: K. Arsenault, updated the s2spost filenaming conventions
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
import copy
import datetime
import os
import sys

# Third-party libraries
# NOTE: pylint cannot see the Dataset class in netCDF4 since the latter is not
# written in Python.  We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module, too-many-branches, too-many-statements, too-many-locals
from netCDF4 import Dataset as nc4_dataset
import numpy as np
import yaml

# pylint: enable=no-name-in-module

_cell_methods = {
    "AvgSurfT_tavg" : "time: mean area: point where land",
    "BaseflowStor_tavg" : "time: mean area: point where land",
    "Elevation_inst" : "area: point where land",
    "Evap_tavg" : "time: mean area: point where land",
    "FloodedArea_tavg" : "time: mean area: point where land",
    "FloodedFrac_tavg" : "time: mean area: point where land",
    "FloodStor_tavg" : "time: mean area: point where land",
    "GWS_tavg" : "time: mean area: point where land",
    "Greenness_inst" : "area: point where land",
    "LWdown_f_tavg" : "time: mean",
    "Psurf_f_tavg" : "time: mean",
    "Qair_f_tavg" : "time: mean",
    "Qs_acc" : "time: sum area: point where land",
    "Qsb_acc" : "time: sum area: point where land",
    "RelSMC_tavg" : "time: mean area: point where land",
    "RiverStor_tavg" : "time: mean area: point where land",
    "Snowcover_tavg" : "time: mean area: point where land",
    "SnowDepth_tavg" : "time: mean area: point where land",
    "SoilMoist_tavg" : "time: mean area: point where land",
    "SoilTemp_tavg" : "time: mean area: point where land",
    "Streamflow_tavg" : "time: mean area: point where land",
    "SWdown_f_tavg" : "time: mean",
    "SWE_tavg" : "time: mean area: point where land",
    "SWS_tavg" : "time: mean area: point where land",
    "Tair_f_max" : "time: maximum",
    "Tair_f_min" : "time: minimum",
    "Tair_f_tavg" : "time: mean",
    "TotalPrecip_acc" : "time: sum",
    "TWS_tavg" : "time: mean area: point where land",
    "Wind_f_tavg" : "time: mean",
}

_new_standard_names = {
    "AvgSurfT_tavg" : "surface_temperature",
    "Elevation_inst" : "height_above_mean_sea_level",
    "Evap_tavg" : "water_evapotranspiration_flux",
    "LANDMASK" : "land_binary_mask",
    "Psurf_f_tavg" : "surface_air_pressure",
    "Qg_tavg" : "downward_heat_flux_at_ground_level_in_soil",
    "SnowDepth_tavg" : "surface_snow_thickness",
    "Soiltype_inst" : "soil_type",
    "Streamflow_tavg" : "water_volume_transport_in_river_channel",
    "TotalPrecip_acc" : "precipitation_amount",
}

_remove_standard_names_list = ["BaseflowStor_tavg", "RelSMC_tavg",
                               "FloodedArea_tavg", "FloodedFrac_tavg",
                               "FloodStor_tavg", "Greenness_inst",
                               "GWS_tavg", "Landcover_inst",
                               "Landmask_inst",
                               "RiverDepth_tavg",
                               "RiverStor_tavg",
                               "SoilMoist_tavg", "Soiltype_inst",
                               "SWS_tavg", "TWS_tavg"]

_new_units = {
    "ensemble" : "1",
    "FloodedFrac_tavg" : "1",
    "Greenness_inst" : "1",
    "LANDMASK" : "1",
    "Landmask_inst" : "1",
    "RelSMC_tavg" : "1",
    "Qair_f_tavg" : "1",
    "SoilMoist_tavg" : "1",
}

# Private methods.
def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} configfile noahmp_file hymap2_file"
    txt += " output_dir fcst_date YYYYMMDDHH model_forcing"
    print(txt)
    print("[INFO] where:")
    print("[INFO] configfile: Path to s2spost config file")
    print("[INFO] noahmp_file: LIS-NoahMP netCDF file (2d ensemble gridspace)")
    txt = "[INFO] hymap2_file: LIS-HYMAP2 netCDF file (2d ensemble gridspace)"
    print(txt)
    print("[INFO] output_dir: Directory to write merged output")
    print("[INFO] fcst_yyyymmdd: Initial forecast date of daily files")
    print("[INFO] YYYYMMDDHH: Valid lead year,month,day,hour of data (in UTC)")
    print("[INFO] model_forcing: ID for atmospheric forcing for LIS")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 8:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if input files exist.
    configfile = sys.argv[1]
    if not os.path.exists(configfile):
        print(f"[ERR] {configfile} does not exist!")
        sys.exit(1)

    noahmp_file = sys.argv[2]
    if not os.path.exists(noahmp_file):
        print(f"[ERR] {noahmp_file} does not exist!")
        sys.exit(1)

    hymap2_file = sys.argv[3]
    if not os.path.exists(hymap2_file):
        print(f"[ERR] {hymap2_file} does not exist!")
        sys.exit(1)

    # Create output directory if it doesn't exist.
    output_dir = sys.argv[4]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get valid initial forecast date.
    fcst_yyyymmdd = sys.argv[5]

    if len(fcst_yyyymmdd) != 8:
        print("[ERR] Invalid length for fcst_yyyymmdd, must be 8 characters!")
        sys.exit(1)
    fcstyear = int(fcst_yyyymmdd[0:4])
    fcstmonth = int(fcst_yyyymmdd[4:6])
    fcstday = int(fcst_yyyymmdd[6:8])

    try:
        fcst_date = datetime.datetime(fcstyear, fcstmonth, fcstday)
    except ValueError:
        print("[ERR] Invalid fcst_yyyymmdd passed to script!")
        sys.exit(1)

    # Get valid lead date and time of data.
    yyyymmddhh = sys.argv[6]

    if len(yyyymmddhh) != 10:
        print("[ERR] Invalid length for YYYYMMDDHH, must be 10 characters!")
        sys.exit(1)
    year = int(yyyymmddhh[0:4])
    month = int(yyyymmddhh[4:6])
    day = int(yyyymmddhh[6:8])
    hour = int(yyyymmddhh[8:10])

    try:
        curdt = datetime.datetime(year, month, day, hour)
    except ValueError:
        print("[ERR] Invalid YYYYMMDDHH passed to script!")
        sys.exit(1)

    # Get ID of model forcing
    model_forcing = sys.argv[7].upper()

    return configfile, noahmp_file, hymap2_file, output_dir, fcst_date, \
        curdt, model_forcing

def _read_config(configfile):
    """Read from s2spost config file."""
    config = configparser.ConfigParser()
    config.read(configfile)
    return config

def _create_final_filename(output_dir, fcst_date, curdt, model_forcing, domain):
    """Create final filename, following 557 convention."""
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing}"
    name += "_GR.C0P25DEG"
    if domain == 'AFRICOM':
        name += "_AR.AFRICA"
    if domain == 'GLOBAL':
        name += "_AR.GLOBAL"

    name += "_PA.ALL"
    name += f"_DD.{fcst_date.year:04d}{fcst_date.month:02d}{fcst_date.day:02d}"
    name += f"_DT.{curdt.hour:02d}00"
    name += f"_FD.{curdt.year:04d}{curdt.month:02d}{curdt.day:02d}"
    name += f"_DT.{curdt.hour:02d}00"
    name += "_DF.NC"
    if len(os.path.basename(name)) > 128:
        print("[ERR] Output file name is too long!")
        print(f"[ERR] {os.path.basename(name)} exceeds 128 characters!")
        sys.exit(1)
    return name

def _merge_files(ldtfile, noahmp_file, hymap2_file, merge_file, fcst_date):
    """Copy LDT, NoahMP and HYMAP2 fields into same file."""

    src1 = nc4_dataset(noahmp_file, "r")
    src2 = nc4_dataset(hymap2_file, "r")
    src3 = nc4_dataset(ldtfile, "r")
    dst = nc4_dataset(merge_file, "w", format="NETCDF4", )

    # Define all dimensions, variables, and attributes from src1
    dimension_dict = {
        "east_west" : "lon",
        "north_south" : "lat",
        "SoilMoist_profiles" : "soil_layer",
        "SoilTemp_profiles" : "soil_layer",
        "RelSMC_profiles" : "soil_layer",
    }
    for dimname in src1.dimensions:
        # Soil dimensions will be replaced by soil_layer, so just
        # write it once.
        if dimname in ["SoilTemp_profiles", \
                       "RelSMC_profiles"]:
            continue
        dimname1 = dimname
        if dimname in dimension_dict:
            dimname1 = dimension_dict[dimname]
        dst.createDimension(dimname1, src1.dimensions[dimname].size)
    attrs = copy.deepcopy(src1.__dict__)
    attrs["Conventions"] = "CF-1.8"
    del attrs["conventions"]
    del attrs["missing_value"]
    del attrs["SOIL_LAYER_THICKNESSES"]
    attrs["source"] = "Noah-MP.4.0.1+template open water+HYMAP2"

    dst.setncatts(attrs)

    for name, variable in src1.variables.items():

        # Special handling for certain variables to match CF convention
        if name == "lat":
            dst.createVariable(name, variable.datatype, ("lat"), zlib=True,
                               complevel=6, shuffle=True)
        elif name == "lon":
            dst.createVariable(name, variable.datatype, ("lon"),zlib=True,
                               complevel=6, shuffle=True)
        elif name in ["SoilMoist_tavg", "SoilTemp_tavg",
                      "RelSMC_tavg"]:
            dst.createVariable(name, variable.datatype, \
                               ("ensemble", "time", "soil_layer","lat", "lon"),zlib=True,
                               complevel=6, shuffle=True)
        else:
            # Need to account for new CF dimension names
            dimensions = []
            for dimension in variable.dimensions:
                if dimension in dimension_dict:
                    dimensions.append(dimension_dict[dimension])
                else:
                    dimensions.append(dimension)
            if len(dimensions) == 3:
                dimensions = ['ensemble', 'time', 'lat', 'lon']
            if name == "Landcover_inst" or name == "Soiltype_inst" or name == "Elevation_inst":
                dimensions = ['lat', 'lon']
            var = dst.createVariable(name, variable.datatype,
                                     dimensions,zlib=True, complevel=6, shuffle=True )
            if name == "Landcover_inst":
                var.flag_values = [np.float32(i) for i in range(1,22)]
            elif name == "Soiltype_inst":
                var.flag_values = [np.float32(i) for i in range(1,17)]
        # Extra CF attributes
        attrs = copy.deepcopy(src1[name].__dict__)
        if name == "time":
            attrs["calendar"] = "standard"
            attrs["axis"] = "T"
            attrs["bounds"] = "time_bnds"
        elif name == "lat":
            attrs["axis"] = "Y"
        elif name == "lon":
            attrs["axis"] = "X"
        elif name == "ensemble":
            attrs["axis"] = "E"
        elif name in ["SoilMoist_tavg"]:
            attrs["long_name"] = "volumetric soil moisture content"
        elif name == "Soiltype_inst":
            attrs["flag_meanings"] = \
                "sand loamy_sand sandy_loam silt_loam silt loam" + \
                " sandy_clay_loam silty_clay_loam clay_loam sandy_clay" + \
                " silty_clay clay organic_material water bedrock" + \
                " other+land-ice"
            attrs["valid_range"] = [1., 16.]
            del attrs["units"]
        elif name == "Landcover_inst":
            attrs["flag_meanings"] = \
                "evergreen_needleleaf_forest evergreen_broadleaf_forest" + \
                " deciduous_needleleaf_forest deciduous_broadleaf_forest" + \
                " mixed_forests closed_shrublands open_shrublands" + \
                " woody_savannas savannas grasslands permanent_wetlands" + \
                " croplands urban_and_built-up" + \
                " cropland+natural_vegetation_mosaic snow_and_ice" + \
                " barren_or_sparsely_vegetated water wooded_tundra" + \
                " mixed_tundra barren_tundra water"
            attrs["valid_range"] = [1., 21.]
            del attrs["units"]
        if name in _cell_methods:
            attrs["cell_methods"] = _cell_methods[name]
        if name in _new_standard_names:
            attrs["standard_name"] = _new_standard_names[name]
        if name in _remove_standard_names_list:
            del attrs["standard_name"]
        if name in _new_units:
            attrs["units"] = _new_units[name]
        dst[name].setncatts(attrs)

    # Add select variables and attributes from src2
    # In future, may want to blend RunoffStor and BaseflowStor with above (KRA)
    src2_excludes = ["lat", "lon", "time", "ensemble", "RunoffStor_tavg",
                     "BaseflowStor_tavg"]
    for name, variable in src2.variables.items():
        if name in src2_excludes:
            continue
        dimensions = []
        for dimension in variable.dimensions:
            if dimension in dimension_dict:
                dimensions.append(dimension_dict[dimension])
            else:
                dimensions.append(dimension)
        if len(dimensions) == 3:
                dimensions = ['ensemble', 'time', 'lat', 'lon']
        dst.createVariable(name, variable.datatype,
                           dimensions,zlib=True, complevel=6, shuffle=True )
        # Extra CF attributes
        attrs = copy.deepcopy(src2[name].__dict__)
        if name in _cell_methods:
            attrs["cell_methods"] = _cell_methods[name]
        if name in _new_standard_names:
            attrs["standard_name"] = _new_standard_names[name]
        if name in _remove_standard_names_list:
            del attrs["standard_name"]
        if name in _new_units:
            attrs["units"] = _new_units[name]
        dst[name].setncatts(attrs)

    # Add LANDMASK variable and attributes from src3
    dimensions = []
    for dimension in src3["LANDMASK"].dimensions:
        if dimension in dimension_dict:
            dimensions.append(dimension_dict[dimension])
        else:
            dimensions.append(dimension)
    dst.createVariable("LANDMASK", src3["LANDMASK"].datatype,
                       dimensions,zlib=True, complevel=6, shuffle=True )
    attrs = copy.deepcopy(src3["LANDMASK"].__dict__)
    attrs["flag_values"] = [np.float32(i) for i in range(0,2)]
    attrs["flag_meanings"] = "water land"
    del attrs["standard_name"]
    attrs["long_name"] = "land mask from LDT"
    dst["LANDMASK"].setncatts(attrs)

    # Add time_bnds variable
    dst.createDimension("nv", 2)
    dst.createVariable("time_bnds", "f4", ("time", "nv"))

    # Add soil_layer and soil_layer_thickness variables
    dst.createVariable("soil_layer", "i4", ("soil_layer"))
    attrs = {
        "long_name" : "soil layer level",
        "axis" : "Z",
        "positive" : "down",
    }
    dst["soil_layer"].setncatts(attrs)
    dst.createVariable("soil_layer_thickness", "f4", ("soil_layer"))
    attrs = {
        "long_name" : "soil layer thicknesses",
        "units" : "m",
    }
    dst["soil_layer_thickness"].setncatts(attrs)

    # add atime forecast_reference_time
    dst.createVariable("atime", "f8")
    attrs = {"standard_name": "forecast_reference_time",
             "units": "hours since " + fcst_date.strftime("%Y-%m-%d") + " 00:00"}
    dst["atime"].setncatts(attrs)

    # Write data from src1
    for name, variable in src1.variables.items():
        # Special handling for lat and lon, which should be 1d arrays
        # for lat/lon projection in CF convention
        if name == "lat":
            dst[name][:] = src1[name][:,0]
        elif name == "lon":
            dst[name][:] = src1[name][0,:]
        # Special handling for soil variables. CF convention requires
        # switching the soil_layer and ensemble dimension.
        elif name in ["SoilMoist_tavg", "SoilTemp_tavg",
                      "RelSMC_tavg"]:
            _ns = dst.dimensions["soil_layer"].size
            _es = dst.dimensions["ensemble"].size
            for e in range(0, _es):
                for i in range(0, _ns):
                    dst[name][e,0,i,:,:] = src1[name][i,e,:,:]
        #elif len(variable.dimensions) == 4:
        #    dst[name][:,0,:,:,:] = src1[name][:,:,:,:]
        elif len(variable.dimensions) == 3:
            if name == "Landcover_inst" or name == "Soiltype_inst" or name == "Elevation_inst":
                dst[name][:,:] = src1[name][0,:,:]
            else:
                dst[name][:,0,:,:] = src1[name][:,:,:]
        elif len(variable.dimensions) == 2:
            dst[name][:,:] = src1[name][:,:]
        elif len(variable.dimensions) == 1:
            dst[name][:] = src1[name][:]
        if name == 'TotalPrecip_acc':
            # apply LANDMASK
            prec = np.array(src1[name][:])
            mask = np.array(src3["LANDMASK"][:])
            nens = prec.shape[0]
            for i in range (0, nens):
                prec[i,:,:] = np.where(mask ==1, prec[i,:,:],-9999.)
            dst[name][:,0,:,:] = prec

    # Write data from src2
    for name, variable in src2.variables.items():
        if name in src2_excludes:
            continue
        if len(variable.dimensions) == 4:
            dst[name][:,0,:,:,:] = src2[name][:,:,:,:]
        elif len(variable.dimensions) == 3:
            dst[name][:,0,:,:] = src2[name][:,:,:]
        elif len(variable.dimensions) == 2:
            dst[name][:,:] = src2[name][:,:]
        elif len(variable.dimensions) == 1:
            dst[name][:] = src2[name][:]

    # Write data from src3
    dst["LANDMASK"][:,:] = src3["LANDMASK"][:,:]

    # Write time_bnds
    dst["time_bnds"][0,0] = -1440.
    dst["time_bnds"][0,1] = 0.

    # writ atime
    dst["atime"][()] = 0.

    # Write soil layer data
    dst["soil_layer"][0] = 1
    dst["soil_layer"][1] = 2
    dst["soil_layer"][2] = 3
    dst["soil_layer"][3] = 4
    dst["soil_layer_thickness"][0] = 0.1
    dst["soil_layer_thickness"][1] = 0.3
    dst["soil_layer_thickness"][2] = 0.6
    dst["soil_layer_thickness"][3] = 1.0
    del dst['time'].time_increment

    src1.close()
    src2.close()
    src3.close()
    dst.close()
    sys.exit()
# Test driver
if __name__ == "__main__":

    # Get the file and directory names
    _configfile, _noahmp_file, _hymap2_file, _output_dir, _fcst_date, _curdt, _model_forcing \
        = _read_cmd_args()
    # load config file
    with open(_configfile, 'r', encoding="utf-8") as file:
        _config = yaml.safe_load(file)

    _final_file = _create_final_filename(_output_dir, _fcst_date, _curdt, _model_forcing,
                                        _config["EXP"]["DOMAIN"])

    _ldtfile = _config['SETUP']['supplementarydir'] + '/lis_darun/' + \
        _config["SETUP"]["ldtinputfile"]
    _merge_files(_ldtfile, _noahmp_file, _hymap2_file, _final_file, _fcst_date)
