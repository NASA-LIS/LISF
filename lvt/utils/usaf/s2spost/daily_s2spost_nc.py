#!/usr/bin/env python3

"""
#------------------------------------------------------------------------------
#
# SCRIPT: daily_s2spost_nc.py
#
# PURPOSE:  Merges daily netCDF output from LIS-NoahMP and LIS-HYMAP2 into
# single, CF-compliant netCDF4 file for distribution.
#
# REQUIREMENTS as of 13 Sep 2021:
# * Python 3.8 or higher.
# * netCDF Operator (NCO) binaries version 5.0.1 or later.
#
# REFERENCES:
# https://cfconventions.org for specifications of NetCDF Climate and Forecast
#   (CF) Metadata Conventions.
# https://nco.sourceforge.net for NCO utilities documentation and source code.
# https://unidata.ucar.edu/software/udunits for documentation on UDUNITS2
#   library, which CF is generally consistent with for unit specifications.
#
# REVISION HISTORY:
# 13 Sep 2021: Eric Kemp (SSAI), first version.
# 14 Sep 2021: Eric Kemp (SSAI), refactored to please cfchecker and pylint.
# 15 Sep 2021: Eric Kemp (SSAI), added reading landmask from LDT param file.
# 16 Sep 2021: Eric Kemp (SSAI), renamed.
# 17 Sep 2021: Eric Kemp (SSAI), swapped order of soil_layer and ensemble
#   dimensions; added ID of model forcing for LIS to filenames.
# 23 Sep 2021: Eric Kemp (SSAI), added SmLiqFrac_inst.
# 27 Oct 2021: Eric Kemp/SSAI, addressed pylint string objections.
# 29 Oct 2021: Eric Kemp/SSAI, added config file.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
import datetime
import os
import subprocess
import sys

# Path to NCO binaries.  Hardwired here due to Air Force security requirements.
# This is intended as an internal constant, hence the name is prefixed with
# "_".
#_NCO_DIR = "/usr/local/other/nco/5.0.1/bin" # On Discover

# NOTE: This script is supposed to be entirely self-contained, and not to be
# used as an imported module.  So all the function names are prefixed with "_"
# to indicate "for internal use."

def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} configfile noahmp_file hymap2_file"
    txt += " output_dir YYYYMMDDHH model_forcing"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]  configfile: Path to s2spost config file")
    print("[INFO]  noahmp_file: LIS-Noah netCDF file (2d ensemble gridspace)")
    txt = "[INFO]  hymap2_file: LIS-HYMAP2 netCDF file (2d ensemble gridspace)"
    print(txt)
    print("[INFO]  output_dir: Directory to write merged output")
    print("[INFO]  YYYYMMDDHH is valid year,month,day,hour of data (in UTC)")
    print("[INFO]  model_forcing: ID for atmospheric forcing for LIS")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 7:
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

    # Get valid date and time of data.
    yyyymmddhh = sys.argv[5]

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
    model_forcing = sys.argv[6]

    return configfile, noahmp_file, hymap2_file, output_dir, curdt, \
        model_forcing

def _read_config(configfile):
    """Read from s2spost config file."""
    config = configparser.ConfigParser()
    config.read(configfile)
    return config

def _check_nco_binaries(ncodir):
    """Check to see if necessary NCO binaries are available."""
    nco_bins = ["ncap2", "ncatted", "ncks", "ncpdq", "ncrename"]
    for nco_bin in nco_bins:
        path = f"{ncodir}/{nco_bin}"
        if not os.path.exists(path):
            print(f"[ERR] Cannot find {path} for converting LIS netCDF4 data!")
            print("[ERR] Make sure NCO package is installed on the system!")
            print("[ERR] And update nco_dir in config file if necessary!")
            sys.exit(1)

def _run_cmd(cmd, error_msg):
    """Handle running shell command and checking error."""
    #print(cmd)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print(error_msg)
        sys.exit(1)

def _create_final_filename(output_dir, curdt, model_forcing):
    """Create final filename, following 557 convention."""
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing}"
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S"
    name += f"_DD.{curdt.year:04d}{curdt.month:02d}{curdt.day:02d}"
    name += f"_DT.{curdt.hour:02d}00"
    name += "_DF.NC"
    if len(os.path.basename(name)) > 128:
        print("[ERR] Output file name is too long!")
        print(f"[ERR] {os.path.basename(name)} exceeds 128 characters!")
        sys.exit(1)
    return name

def _merge_files(ncodir, ldtfile, noahmp_file, hymap2_file, merge_file):
    """Copy LDT, NoahMP and HYMAP2 fields into same file."""

    # Copy NoahMP fields.
    cmd = f"{ncodir}/ncks {noahmp_file} -6 {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncks!")

    # Exclude coordinate variables from HYMAP2 file, since these are already
    # copied from NoahMP.  Also exclude unnecessary HYMAP2 storage fields.
    cmd = f"{ncodir}/ncks -A -C"
    cmd += " -x -v lat,lon,time,ensemble,RunoffStor_tavg,BaseflowStor_tavg"
    cmd += f" {hymap2_file} -6 {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncks!")

    # Only copy the LANDMASK field from the LDT parameter file.  This version
    # actually discriminates land and water.
    cmd = f"{ncodir}/ncks -A -C -v LANDMASK {ldtfile} -6 {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncks!")

def _add_time_attrs(ncodir, merge_file):
    """Add calendar, axis, and bounds attributes to time."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a calendar,time,c,c,'standard'"
    cmd += " -a axis,time,c,c,'T'"
    cmd += " -a bounds,time,c,c,'time_bnds'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _add_time_bnds(ncodir, merge_file):
    """Add a time_bnds variable and appropriate dimension."""
    cmd = f"{ncodir}/ncap2"
    cmd +=""" -s 'defdim("nv",2)'"""
    cmd += " -s 'time_bnds[$time,$nv]={-1440f, 0f}'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncap2!")

def _add_soil_layer_data(ncodir, merge_file):
    """Add a soil_layer dimension and variable, plus soil layer thickness."""
    cmd = f"{ncodir}/ncap2"
    cmd +=""" -s 'defdim("soil_layer",$RelSMC_profiles.size)'"""
    cmd += " -s 'soil_layer[$soil_layer]={1,2,3,4}'"
    cmd += """ -s 'soil_layer@long_name="soil layer level"'"""
    cmd += """ -s 'soil_layer@axis="Z"'"""
    cmd += """ -s 'soil_layer@positive="down"'"""
    cmd += " -s 'soil_layer_thickness[$soil_layer]={0.1f, 0.3f, 0.6f, 1.0f}'"
    cmd += """ -s 'soil_layer_thickness@long_name="soil layer thicknesses"'"""
    cmd += """ -s 'soil_layer_thickness@units="m"'"""
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncap2!")

def _rename_lat_lon_dims(ncodir, merge_file):
    """Rename the spatial dimensions to match CF convention."""
    cmd = f"{ncodir}/ncrename -O"
    cmd += " -d east_west,lon"
    cmd += " -d north_south,lat"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncrename!")

def _rename_lat_lon_fields(ncodir, merge_file):
    """Rename lat, lon fields, so we can later create 1d lat, lon arrays."""
    cmd = f"{ncodir}/ncrename -O"
    cmd += " -v lon,old_lon"
    cmd += " -v lat,old_lat"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncrename!")

def _add_axis_lat_lon(ncodir, merge_file):
    """Add axis attributes to to old_lat and old_lon."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a axis,old_lat,c,c,'Y'"
    cmd += " -a axis,old_lon,c,c,'X'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _create_1d_lat_lon_fields(ncodir, merge_file):
    """Create 1d lat and lon coordinate arrays."""
    cmd = f"{ncodir}/ncap2 -v"
    cmd += " -s 'lon=old_lon(1,:)'"
    cmd += " -s 'lat=old_lat(:,1)'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncap2!")

def _add_cell_methods_1(cmd):
    """First set of cell_method additions.  Split this way to appease
    pylint."""
    cmd += \
        " -a cell_methods,Albedo_tavg,c,c,'time: mean area: point where land'"
    cmd += " -a cell_methods,AvgSurfT_inst,c,c,'area: point where land'"
    cmd += \
       " -a cell_methods,AvgSurfT_tavg,c,c,'time: mean area: point where land'"
#    cmd += " -a cell_methods,BaseflowStor_tavg,c,c,"
#    cmd +=        "'time: mean area: point where land'"
    cmd += " -a cell_methods,CanopInt_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,Elevation_inst,c,c,'area: point where land'"
    cmd += \
        " -a cell_methods,Evap_tavg,c,c,'time: mean area: point where land'"
    cmd += " -a cell_methods,FloodedArea_tavg,c,c,"
    cmd +=        "'time: mean area: point where land'"
    cmd += " -a cell_methods,FloodedFrac_tavg,c,c,"
    cmd +=        "'time: mean area: point where land'"
    cmd += " -a cell_methods,FloodStor_tavg,c,c,"
    cmd +=        "'time: mean area: point where land'"
    cmd += " -a cell_methods,GWS_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,Greenness_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,LWdown_f_tavg,c,c,'time: mean'"
    cmd += " -a cell_methods,Psurf_f_tavg,c,c,'time: mean'"
    cmd += " -a cell_methods,Qair_f_tavg,c,c,'time: mean'"
    cmd += " -a cell_methods,Qg_tavg,c,c,'time: mean area: point where land'"
    cmd += " -a cell_methods,Qh_tavg,c,c,'time: mean area: point where land'"
    cmd += " -a cell_methods,Qle_tavg,c,c,'time: mean area: point where land'"
    cmd += " -a cell_methods,Qs_acc,c,c,'time: sum area: point where land'"
    cmd += " -a cell_methods,Qsb_acc,c,c,'time: sum area: point where land'"
    cmd += " -a cell_methods,RelSMC_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,RHMin_inst,c,c,'area: point where land'"
    return cmd

def _add_cell_methods_2(cmd):
    """Second set of cell_method additions.  Split this way to appease
    pylint."""
    cmd += " -a cell_methods,RiverDepth_tavg,c,c,"
    cmd +=        "'time: mean area: point where land'"
    cmd += " -a cell_methods,RiverFlowVelocity_tavg,c,c,"
    cmd +=        "'time: mean area: point where land'"
    cmd += " -a cell_methods,RiverStor_tavg,c,c,"
    cmd +=        "'time: mean area: point where land'"
#    cmd += " -a cell_methods,RunoffStor_tavg,c,c,"
#    cmd +=        "'time: mean area: point where land'"
    cmd += " -a cell_methods,SmLiqFrac_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,Snowcover_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,SnowDepth_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,SoilMoist_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,SoilMoist_tavg,c,c,"
    cmd +=       "'time: mean area: point where land'"
    cmd += " -a cell_methods,SoilTemp_inst,c,c,'area: point where land'"
    cmd += \
       " -a cell_methods,SoilTemp_tavg,c,c,'time: mean area: point where land'"
    cmd += " -a cell_methods,Streamflow_tavg,c,c,"
    cmd +=       "'time: mean area: point where land'"
    cmd += \
       " -a cell_methods,SurfElev_tavg,c,c,'time: mean area: point where land'"
    cmd += " -a cell_methods,SWdown_f_tavg,c,c,'time: mean'"
    cmd += " -a cell_methods,SWE_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,SWS_tavg,c,c,'time: mean area: point where land'"
    cmd += " -a cell_methods,Tair_f_max,c,c,'time: maximum'"
    cmd += " -a cell_methods,Tair_f_min,c,c,'time: minimum'"
    cmd += " -a cell_methods,Tair_f_tavg,c,c,'time: mean'"
    cmd += " -a cell_methods,TotalPrecip_acc,c,c,'time: sum'"
    cmd += " -a cell_methods,TWS_inst,c,c,'area: point where land'"
    cmd += " -a cell_methods,Wind_f_tavg,c,c,'time: mean'"
    return cmd

def _add_cell_methods(ncodir, merge_file):
    """Add the cell_methods attribute to select variables."""
    cmd = f"{ncodir}/ncatted"
    cmd = _add_cell_methods_1(cmd)
    cmd = _add_cell_methods_2(cmd)
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _modify_standard_names(ncodir, merge_file):
    """Modify standard_name attributes for certain variables."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a standard_name,AvgSurfT_inst,m,c,'surface_temperature'"
    cmd += " -a standard_name,AvgSurfT_tavg,m,c,'surface_temperature'"
    cmd += " -a standard_name,CanopInt_inst,m,c,'canopy_water_amount'"
    cmd += " -a standard_name,Elevation_inst,m,c,'height_above_mean_sea_level'"
    cmd += " -a standard_name,Evap_tavg,m,c,'water_evapotranspiration_flux'"
    cmd += " -a standard_name,LANDMASK,m,c,'land_binary_mask'"
    cmd += " -a standard_name,Psurf_f_inst,m,c,'surface_air_pressure'"
    cmd += " -a standard_name,Psurf_f_tavg,m,c,'surface_air_pressure'"
    cmd += " -a standard_name,Qg_tavg,m,c,"
    cmd +=      "'downward_heat_flux_at_ground_level_in_soil'"
    cmd += " -a standard_name,SnowDepth_inst,m,c,'surface_snow_thickness'"
    cmd += " -a standard_name,Soiltype_inst,m,c,'soil_type'"
    cmd += " -a standard_name,Streamflow_tavg,m,c,"
    cmd +=      "'water_volume_transport_in_river_channel'"
    cmd += " -a standard_name,TotalPrecip_acc,m,c,'precipitation_amount'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _remove_standard_names(ncodir, merge_file):
    """Remove standard_name attributes for non-standard variables."""
    cmd = f"{ncodir}/ncatted"
#    cmd += " -a standard_name,BaseflowStor_tavg,d,,"
    cmd += " -a standard_name,RelSMC_inst,d,,"
    cmd += " -a standard_name,FloodedArea_tavg,d,,"
    cmd += " -a standard_name,FloodedFrac_tavg,d,,"
    cmd += " -a standard_name,FloodStor_tavg,d,,"
    cmd += " -a standard_name,Greenness_inst,d,,"
    cmd += " -a standard_name,GWS_inst,d,,"
    cmd += " -a standard_name,Landcover_inst,d,,"
    cmd += " -a standard_name,Landmask_inst,d,,"
    cmd += " -a standard_name,RelSMC_inst,d,,"
    cmd += " -a standard_name,RHMin_inst,d,,"
    cmd += " -a standard_name,RiverDepth_tavg,d,,"
    cmd += " -a standard_name,RiverFlowVelocity_tavg,d,,"
    cmd += " -a standard_name,RiverStor_tavg,d,,"
#    cmd += " -a standard_name,RunoffStor_tavg,d,,"
    cmd += " -a standard_name,SmLiqFrac_inst,d,,"
    cmd += " -a standard_name,SoilMoist_inst,d,,"
    cmd += " -a standard_name,SoilMoist_tavg,d,,"
    cmd += " -a standard_name,SurfElev_tavg,d,,"
    cmd += " -a standard_name,SWS_tavg,d,,"
    cmd += " -a standard_name,TWS_inst,d,,"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _modify_units(ncodir, merge_file):
    """Modify units for certain dimensionless variables."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a units,ensemble,m,c,'1'"
    cmd += " -a units,FloodedFrac_tavg,m,c,'1'"
    cmd += " -a units,Greenness_inst,m,c,'1'"
    cmd += " -a units,LANDMASK,m,c,'1'"
    cmd += " -a units,Landmask_inst,m,c,'1'"
    cmd += " -a units,RelSMC_inst,m,c,'1'"
    cmd += " -a units,Qair_f_inst,m,c,'1'"
    cmd += " -a units,Qair_f_tavg,m,c,'1'"
    cmd += " -a units,SmLiqFrac_inst,m,c,'1'"
    cmd += " -a units,SoilMoist_inst,m,c,'1'"
    cmd += " -a units,SoilMoist_tavg,m,c,'1'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _modify_soilmoist_long_name(ncodir, merge_file):
    """Modify long_name of SoilMoist fields to emphasize volumetric units."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a long_name,SoilMoist_inst,m,c,"
    cmd +=    "'volumetric soil moisture content'"
    cmd += " -a long_name,SoilMoist_tavg,m,c,"
    cmd +=    "'volumetric soil moisture content'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_landmask_attrs(ncodir, merge_file):
    """Add metadata defining the landmask."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a flag_values,LANDMASK,c,f,"
    cmd +=      "'0, 1'"
    cmd += " -a flag_meanings,LANDMASK,c,c,"
    cmd +=    "'water land'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_soiltype_attrs(ncodir, merge_file):
    """Add and correct metadata defining the soil types."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a flag_values,Soiltype_inst,c,f,"
    cmd +=      "'1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16'"
    cmd += " -a flag_meanings,Soiltype_inst,c,c,"
    cmd +=    "'sand loamy_sand sandy_loam silt_loam silt loam sandy_clay_loam"
    cmd +=    " silty_clay_loam clay_loam sandy_clay silty_clay clay"
    cmd +=    " organic_material water bedrock other+land-ice'"
    cmd += " -a valid_range,Soiltype_inst,c,f,"
    cmd +=      "'1, 16'"
    cmd += " -a units,Soiltype_inst,d,,"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_landcover_attrs(ncodir, merge_file):
    """Add metadata defining the land cover categories."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a flag_values,Landcover_inst,c,f,"
    cmd +=      "'1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,"
    cmd +=      "12, 13, 14, 15, 16, 17, 18, 19, 20, 21'"
    cmd += " -a flag_meanings,Landcover_inst,c,c,"
    cmd +=   "'evergreen_needleleaf_forest evergreen_broadleaf_forest"
    cmd +=   " deciduous_needleleaf_forest deciduous_broadleaf_forest"
    cmd +=   " mixed_forests closed_shrublands open_shrublands"
    cmd +=   " woody_savannas savannas grasslands permanent_wetlands"
    cmd +=   " croplands urban_and_built-up cropland+natural_vegetation_mosaic"
    cmd +=   " snow_and_ice barren_or_sparsely_vegetated water wooded_tundra"
    cmd +=   " mixed_tundra barren_tundra water'"
    cmd += " -a valid_range,Landcover_inst,c,i,"
    cmd +=      "'1, 21'"
    cmd += " -a units,Landcover_inst,d,,"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _update_global_attrs(ncodir, merge_file):
    """Update global attributes."""
    cmd = f"{ncodir}/ncatted"
    cmd += " -a Conventions,global,c,c,'CF-1.8'"
    cmd += " -a conventions,global,d,,"
    cmd += " -a missing_value,global,d,,"
    cmd += " -a SOIL_LAYER_THICKNESSES,global,d,,"
    cmd += " -a source,global,m,c,"
    cmd += "'Noah-MP.4.0.1+template open water+HYMAP2'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncatted!")

def _rename_soil_moist_temp_fields(ncodir, merge_file):
    """Rename the soil moisture and temperature variables so we can create
       new copies with common vertical dimension."""
    cmd = f"{ncodir}/ncrename -O"
    cmd += " -v RelSMC_inst,old_RelSMC_inst"
    cmd += " -v SmLiqFrac_inst,old_SmLiqFrac_inst"
    cmd += " -v SoilMoist_inst,old_SoilMoist_inst"
    cmd += " -v SoilMoist_tavg,old_SoilMoist_tavg"
    cmd += " -v SoilTemp_inst,old_SoilTemp_inst"
    cmd += " -v SoilTemp_tavg,old_SoilTemp_tavg"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncrename!")

def _create_new_soil_moist_temp_fields(ncodir, merge_file):
    """Create new copies of soil fields with new vertical coordinate."""
    cmd = f"{ncodir}/ncap2 -v"
    cmd += " -s 'RelSMC_inst[$soil_layer,$ensemble,$lat,$lon]="
    cmd += "old_RelSMC_inst(:,:,:,:)'"
    cmd += " -s 'SmLiqFrac_inst[$soil_layer,$ensemble,$lat,$lon]="
    cmd += "old_SmLiqFrac_inst(:,:,:,:)'"
    cmd += " -s 'SoilMoist_inst[$soil_layer,$ensemble,$lat,$lon]="
    cmd += "old_SoilMoist_inst(:,:,:,:)'"
    cmd += " -s 'SoilMoist_tavg[$soil_layer,$ensemble,$lat,$lon]="
    cmd += "old_SoilMoist_tavg(:,:,:,:)'"
    cmd += " -s 'SoilTemp_inst[$soil_layer,$ensemble,$lat,$lon]="
    cmd += "old_SoilTemp_inst(:,:,:,:)'"
    cmd += " -s 'SoilTemp_tavg[$soil_layer,$ensemble,$lat,$lon]="
    cmd += "old_SoilTemp_tavg(:,:,:,:)'"
    cmd += f" {merge_file}"
    _run_cmd(cmd, "[ERR] Problem with ncap2!")

def _reverse_soil_layer_ensemble_dims(ncodir, merge_file):
    """Reverse the order of the soil_layer and ensemble dimensions."""
    tmp_file = f"{merge_file}.new"
    cmd = f"{ncodir}/ncpdq"
    cmd += " -a ensemble,soil_layer"
    cmd += f" {merge_file} {tmp_file}"
    _run_cmd(cmd, "[ERR] Problem with ncpdq!")
    os.rename(tmp_file, merge_file)

def _copy_to_final_file(ncodir, merge_file, final_file):
    """Copy data to new netCDF4 file, excluding "old" variables.  This should
    also screen out unneeded dimensions."""
    cmd = f"{ncodir}/ncks"
    cmd += " -x -v old_RelSMC_inst,old_SmLiqFrac_inst,"
    cmd += "old_SoilMoist_inst,old_SoilMoist_tavg,"
    cmd += "old_SoilTemp_inst,old_SoilTemp_tavg,old_lat,old_lon"
    cmd += f" {merge_file} -7 -L 1 {final_file}"
    _run_cmd(cmd, "[ERR] Problem with ncks!")

def _driver():
    """Main driver."""

    # Get the file and directory names
    configfile, noahmp_file, hymap2_file, output_dir, curdt, model_forcing \
        = _read_cmd_args()
    merge_file = f"{output_dir}/merge_tmp.nc"
    final_file = _create_final_filename(output_dir, curdt, model_forcing)
    config = _read_config(configfile)
    ncodir = config["s2spost"]["nco_dir"]
    ldtfile = config["s2spost"]["ldt_file"]

    # Make sure we can find the required NCO binaries.
    _check_nco_binaries(ncodir)

    # Merge the LDT, NoahMP, and HYMAP2 fields into the same file.
    _merge_files(ncodir, ldtfile, noahmp_file, hymap2_file, merge_file)

    # Incremental updates to the coordinate variables, including adding
    # soil_layer as a new vertical coordinate, and renaming dimensions.
    _add_time_attrs(ncodir, merge_file)
    _add_time_bnds(ncodir, merge_file)
    _add_soil_layer_data(ncodir, merge_file)
    _rename_lat_lon_dims(ncodir, merge_file)
    _rename_lat_lon_fields(ncodir, merge_file)
    _add_axis_lat_lon(ncodir, merge_file)
    _create_1d_lat_lon_fields(ncodir, merge_file)

    # Additional updates to metadata.
    _add_cell_methods(ncodir, merge_file)
    _modify_standard_names(ncodir, merge_file)
    _remove_standard_names(ncodir, merge_file)
    _modify_units(ncodir, merge_file)
    _modify_soilmoist_long_name(ncodir, merge_file)
    _update_landmask_attrs(ncodir, merge_file)
    _update_soiltype_attrs(ncodir, merge_file)
    _update_landcover_attrs(ncodir, merge_file)
    _update_global_attrs(ncodir, merge_file)

    # Create new soil moisture and temperature variables with a common
    # vertical coordinate.  This is a two-step process.
    _rename_soil_moist_temp_fields(ncodir, merge_file)
    _create_new_soil_moist_temp_fields(ncodir, merge_file)

    # Reverse order of soil_layer and ensemble dimensions, as recommended by
    # CF convention.
    _reverse_soil_layer_ensemble_dims(ncodir, merge_file)

    # Copy most data to final file, filtering out redundant older fields.
    _copy_to_final_file(ncodir, merge_file, final_file)

    # Delete temporary merge file.
    os.unlink(merge_file)

# Invoke the main driver
if __name__ == "__main__":
    _driver()
