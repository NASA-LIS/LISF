#!/usr/bin/env python3

"""
#------------------------------------------------------------------------------
#
# SCRIPT: monthly_s2spost_nc.py
#
# PURPOSE: Read daily S2S CF-convention netCDF files, calculate monthly
# averages and accumulations, and write to new CF-convention netCDF file.
#
# REQUIREMENTS as of 16 Sep 2021:
# * Python 3.8 or higher
# * UNIDATA NetCDF4 Python library
#
# REVISION HISTORY:
# 16 Sep 2021: Eric Kemp (SSAI), first version.
# 17 Sep 2021: Eric Kemp (SSAI), renamed script, tweaked variable list and
#   attributes; added ID for model forcing for LIS to filenames.
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import sys
import time

# Third-party libraries
# NOTE: pylint cannot see the Dataset class in netCDF4 since the latter is not
# written in Python.  We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module

# Lists of variables to process.  These are intended as internal constants,
# hence the names are prefixed with "_".

_VAR_ACC_LIST = ["Qs_acc", "Qsb_acc", "TotalPrecip_acc"]

_VAR_TAVG_LAND_LIST = ["Qle_tavg", "Qh_tavg", "Qg_tavg", "Evap_tavg",
                       "AvgSurfT_tavg", "Albedo_tavg",
                       "SoilMoist_tavg", "SoilTemp_tavg",
                       "Streamflow_tavg", "FloodedFrac_tavg", "SurfElev_tavg",
                       "SWS_tavg", "RiverStor_tavg", "RiverDepth_tavg",
                       "RiverFlowVelocity_tavg", "FloodStor_tavg",
                       "FloodedArea_tavg"]

_VAR_TAVG_F_LIST = ["Wind_f_tavg", "Tair_f_tavg",
                    "Qair_f_tavg", "Psurf_f_tavg",
                    "SWdown_f_tavg", "LWdown_f_tavg"]


_VAR_TAVG_TWSGWS_LIST = ["TWS_inst", "GWS_inst"]

_VAR_TAIR_MAX_LIST = ["Tair_f_max"]
_VAR_TAIR_MIN_LIST = ["Tair_f_min"]

_VAR_INST_LIST = ["AvgSurfT_inst", "SWE_inst", "SnowDepth_inst",
                  "SoilMoist_inst", "SoilTemp_inst",
                  "SmLiqFrac_inst",
                  "CanopInt_inst",
                  "Snowcover_inst", "Wind_f_inst", "Tair_f_inst",
                  "Qair_f_inst", "Psurf_f_inst",
                  "SWdown_f_inst", "LWdown_f_inst",
                  "Greenness_inst", "RelSMC_inst"]

_CONST_LIST = ["lat", "lon", "ensemble", "soil_layer",
               "soil_layer_thickness",
               "Landmask_inst", "LANDMASK", "Landcover_inst",
               "Soiltype_inst", "Elevation_inst"]

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s input_dir output_dir start_yyyymmdd end_yyyymmdd"
    txt += " model_forcing"
    print(txt)
    print("[INFO] where:")
    print("[INFO]  input_dir: directory with daily S2S files in CF convention")
    print("[INFO]  output_dir: directory for output file")
    print("[INFO]  start_yyyymmdd: Starting date/time of daily files")
    print("[INFO]  end_yyyymmdd: Ending date/time of daily files")
    print("[INFO]  model_forcing: ID for atmospheric forcing for LIS")

def _proc_date(yyyymmdd):
    """Convert YYYYMMDD string to Python date object."""
    if len(yyyymmdd) != 8:
        print("[ERR] Invalid length for YYYYMMDD, must be 8 characters!")
        sys.exit(1)
    year = int(yyyymmdd[0:4])
    month = int(yyyymmdd[4:6])
    day = int(yyyymmdd[6:8])
    try:
        dateobj = datetime.date(year, month, day)
    except ValueError:
        print("[ERR] Invalid YYYYMMDD passed to script!")
        sys.exit(1)
    return dateobj

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 6:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if input directory exists.
    input_dir = sys.argv[1]
    if not os.path.exists(input_dir):
        print("[ERR] Directory %s does not exist!")
        sys.exit(1)

    # Create output directory if it doesn't exist.
    output_dir = sys.argv[2]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get valid starting and ending dates of data.
    start_yyyymmdd = sys.argv[3]
    startdate = _proc_date(start_yyyymmdd)
    end_yyyymmdd = sys.argv[4]
    enddate = _proc_date(end_yyyymmdd)
    if startdate > enddate:
        print("[ERR] Start date is after end date!")
        sys.exit(1)

    # Get ID for model forcing for LIS
    model_forcing = sys.argv[5]

    return input_dir, output_dir, startdate, enddate, model_forcing

def _check_filename_size(name):
    """Make sure filename does not exceed 128 characters, per Air Force
    requirement."""
    if len(os.path.basename(name)) > 128:
        print("[ERR] Output file name is too long!")
        print("[ERR] %s exceeds 128 characters!" %(os.path.basename(name)))
        sys.exit(1)

def _create_daily_s2s_filename(input_dir, curdate, model_forcing):
    """Create path to daily S2S netCDF file."""
    name = "%s" %(input_dir)
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += "_GP.LIS-S2S-%s" %(model_forcing)
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S"
    name += "_DD.%4.4d%2.2d%2.2d" %(curdate.year, curdate.month, curdate.day)
    name += "_DT.0000"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

def _create_monthly_s2s_filename(output_dir, startdate, enddate,
                                 model_forcing):
    """Create path to monthly S2S netCDF file."""
    name = "%s" %(output_dir)
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += "_GP.LIS-S2S-%s" %(model_forcing)
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S"
    name += "_DP.%4.4d%2.2d%2.2d-%4.4d%2.2d%2.2d" \
        %(startdate.year, startdate.month, startdate.day,
          enddate.year, enddate.month, enddate.day)
    name += "_TP.0000-0000"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

def _copy_dims_gattrs(ncid_in, ncid_out):
    """Copy dimensions and global attributes from one netCDF file to
    another."""
    for dimname in ncid_in.dimensions:
        ncid_out.createDimension(dimname, ncid_in.dimensions[dimname].size)
    for gattrname in ncid_in.__dict__:
        ncid_out.setncattr(gattrname, ncid_in.__dict__[gattrname])

def _create_firstguess_monthly_file(infile, outfile):
    """Read daily S2S file, and most fields copy to monthly S2S file.
    This allows us to cleanly copy dimensions and all attributes.  The
    numerical values of the arrays in the monthly S2S
    file will be replaced later in the script."""

    if not os.path.exists(infile):
        print("[ERR] %s does not exist!" %(infile))
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    # Create monthly file, copying dimensions and global attributes.
    # NOTE:  We clean up the global attributes later in the script.
    ncid_out = nc4_dataset(outfile, "w", format='NETCDF4_CLASSIC')
    _copy_dims_gattrs(ncid_in, ncid_out)

    # Copy the fields.
    varnames = _CONST_LIST + _VAR_INST_LIST + _VAR_ACC_LIST + \
        _VAR_TAVG_LAND_LIST + _VAR_TAVG_F_LIST + _VAR_TAVG_TWSGWS_LIST + \
        _VAR_TAIR_MAX_LIST + _VAR_TAIR_MIN_LIST
    for varname in varnames:
        var_in = ncid_in.variables[varname]
        if "missing_value" in var_in.__dict__:
            var_out = \
                ncid_out.createVariable(varname, var_in.datatype,
                                        dimensions=var_in.dimensions,
                                        zlib=True,
                                        complevel=1,
                                        shuffle=True,
                                        fill_value=var_in.missing_value)
        else:
            var_out = \
                ncid_out.createVariable(varname, var_in.datatype,
                                        dimensions=var_in.dimensions,
                                        zlib=True,
                                        complevel=1,
                                        shuffle=True)

        for attrname in var_in.__dict__:
            if attrname == "_FillValue":
                continue
            var_out.setncattr(attrname, var_in.__dict__[attrname])
        if len(var_out.shape) == 4:
            var_out[:,:,:,:] = var_in[:,:,:,:]
        elif len(var_out.shape) == 3:
            var_out[:,:,:] = var_in[:,:,:]
        elif len(var_out.shape) == 2:
            var_out[:,:] = var_in[:,:]
        elif len(var_out.shape) == 1:
            var_out[:] = var_in[:]

    ncid_out.close()
    ncid_in.close()

def _read_second_daily_file(infile):
    """Read the second daily S2S file and copy the acc and tavg fields in
    appropriate dictionaries. We use the second file to start, since acc
    and tavg are valid for the prior 24-hr period."""

    if not os.path.exists(infile):
        print("[ERR] %s does not exist!" %(infile))
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    accs = {}
    tavgs = {}
    tavgs["counter"] = 1

    # Copy the values of the fields we will average or accumulate
    varnames = _VAR_TAVG_LAND_LIST + _VAR_TAVG_F_LIST + \
        _VAR_TAVG_TWSGWS_LIST + _VAR_TAIR_MAX_LIST + _VAR_TAIR_MIN_LIST
    for varname in varnames:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            tavgs[varname] = var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            tavgs[varname] = var_in[:,:,:]
        elif len(var_in.shape) == 2:
            tavgs[varname] = var_in[:,:]
    for varname in _VAR_ACC_LIST:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            accs[varname] = var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            accs[varname] = var_in[:,:,:]
        elif len(var_in.shape) == 2:
            accs[varname] = var_in[:,:]

    ncid_in.close()
    return accs, tavgs

def _read_next_daily_file(infile, accs, tavgs):
    """Read next daily S2S file and copy the required variable values to
    appropriate dictionaries."""

    if not os.path.exists(infile):
        print("[ERR] %s does not exist!" %(infile))
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    tavgs["counter"] += 1

    # Add the values of the fields we will average or accumulate
    varnames = _VAR_TAVG_LAND_LIST + _VAR_TAVG_F_LIST + \
        _VAR_TAVG_TWSGWS_LIST + _VAR_TAIR_MAX_LIST + _VAR_TAIR_MIN_LIST
    for varname in varnames:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            tavgs[varname][:,:,:,:] += var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            tavgs[varname][:,:,:] += var_in[:,:,:]
        elif len(var_in.shape) == 2:
            tavgs[varname][:,:] += var_in[:,:]
    for varname in _VAR_ACC_LIST:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            accs[varname][:,:,:,:] = var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            accs[varname][:,:,:] = var_in[:,:,:]
        elif len(var_in.shape) == 2:
            accs[varname][:,:] = var_in[:,:]

    ncid_in.close()
    return accs, tavgs

def _finalize_tavgs(tavgs):
    """Finalize averages by dividing by sample size."""
    count = tavgs["counter"]
    del tavgs["counter"]
    for varname in tavgs:
        if varname == "counter":
            continue
        if len(tavgs[varname].shape) == 4:
            tavgs[varname][:,:,:,:] /= count
        elif len(tavgs[varname].shape) == 3:
            tavgs[varname][:,:,:] /= count
        elif len(tavgs[varname].shape) == 2:
            tavgs[varname][:,:] /= count
    return tavgs

def _update_monthly_s2s_values(outfile, accs, tavgs):
    """Update the values in the monthly S2S file."""
    if not os.path.exists(outfile):
        print("[ERR] %s does not exist!" %(outfile))
        sys.exit(1)
    ncid = nc4_dataset(outfile, 'a', format='NETCDF4_CLASSIC')
    for dictionary in [accs, tavgs]:
        for varname in dictionary:
            var = ncid.variables[varname]
            if len(var.shape) == 4:
                var[:,:,:,:] = dictionary[varname][:,:,:,:]
            elif len(var.shape) == 3:
                var[:,:,:] = dictionary[varname][:,:,:]
            elif len(var.shape) == 2:
                var[:,:] = dictionary[varname][:,:]
    ncid.close()

def _add_time_data(infile, outfile, startdate, enddate):
    """Add time information to outfile, matching CF convention.  This
    requires pulling more data from the last daily file."""

    if not os.path.exists(infile):
        print("[ERR] %s does not exist!" %(infile))
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')
    if not os.path.exists(outfile):
        print("[ERR] %s does not exist!" %(outfile))
        sys.exit(1)
    ncid_out = nc4_dataset(outfile, 'a', format='NETCDF4_CLASSIC')

    # Copy the time array from the last daily file.
    var_in = ncid_in.variables["time"]
    var_out = ncid_out.createVariable("time", var_in.datatype,
                                      dimensions=var_in.dimensions,
                                      zlib=True,
                                      complevel=1,
                                      shuffle=True)
    for attrname in var_in.__dict__:
        if attrname == "_FillValue":
            continue
        var_out.setncattr(attrname, var_in.__dict__[attrname])
    var_out[:] = var_in[:]

    # Copy the time_bnds array from the last daily file.  But, we will change
    # the value to span one month of data.
    var_in = ncid_in.variables["time_bnds"]
    var_out = ncid_out.createVariable("time_bnds", var_in.datatype,
                                      dimensions=var_in.dimensions,
                                      zlib=True,
                                      complevel=1,
                                      shuffle=True)
    var_out[:,:] = var_in[:,:]
    # Count number of minutes between start and end dates.
    var_out[0,0] = ((enddate - startdate).days)  * (-24*60) # Days to minutes
    var_out[0,1] = 0

    ncid_in.close()
    ncid_out.close()

def _update_cell_methods(outfile):
    """Update cell_method attributes for select variables."""

    if not os.path.exists(outfile):
        print("[ERR] %s does not exist!" %(outfile))
        sys.exit(1)
    ncid = nc4_dataset(outfile, 'a', format='NETCDF4_CLASSIC')

    # Elaborate on monthly calculations of most variables
    varnames = _VAR_ACC_LIST + _VAR_TAVG_LAND_LIST + _VAR_TAVG_F_LIST + \
        _VAR_TAVG_TWSGWS_LIST + _VAR_TAIR_MAX_LIST + _VAR_TAIR_MIN_LIST

    for varname in varnames:

        var = ncid.variables[varname]

        # Special handling for TWS_inst and GWS_inst -- we want the monthly
        # means.
        if varname in _VAR_TAVG_TWSGWS_LIST:
            var.cell_methods = \
                "time: mean (interval: 1 day) area: point where land"

        # Special handling for Tair_f_max and Tair_f_min:  We want the monthly
        # averages of the daily maxs and mins
        elif varname in _VAR_TAIR_MAX_LIST:
            var.cell_methods = \
                "time: mean (interval: 1 day comment: daily maxima)"
        elif varname in _VAR_TAIR_MIN_LIST:
            var.cell_methods = \
                "time: mean (interval: 1 day comment: daily minima)"

        # Clarify monthly mean of daily means over land
        elif varname in _VAR_TAVG_LAND_LIST:
            var.cell_methods = \
                "time: mean (interval: 1 day comment: daily means)" + \
                " area: point where land"

        # Clarify monthly mean of daily means everywhere
        elif varname in _VAR_TAVG_F_LIST:
            var.cell_methods = \
                "time: mean (interval: 1 day comment: daily means)"

        # Clarify monthly accumulations
        elif varname in _VAR_ACC_LIST:
            var.cell_methods = \
                "time: sum (interval: 1 day comment: daily sums)"

    ncid.close()

def _cleanup_global_attrs(outfile):
    """Clean-up global attributes."""
    if not os.path.exists(outfile):
        print("[ERR] %s does not exist!" %(outfile))
        sys.exit(1)
    ncid = nc4_dataset(outfile, 'a', format='NETCDF4_CLASSIC')
    ncid.history = "created on date: %s" %(time.ctime())
    del ncid.NCO
    del ncid.history_of_appended_files

def _driver():
    """Main driver."""

    # Get the directories and dates
    input_dir, output_dir, startdate, enddate, model_forcing = _read_cmd_args()

    # Loop through dates
    curdate = startdate
    seconddate = startdate + datetime.timedelta(days=1)
    delta = datetime.timedelta(days=1)
    while curdate <= enddate:
        infile = _create_daily_s2s_filename(input_dir, curdate, model_forcing)
        #print("[INFO] Reading %s" %(infile))
        if curdate == startdate:
            tmp_outfile = "%s/tmp_monthly.nc" %(output_dir)
            _create_firstguess_monthly_file(infile, tmp_outfile)
        elif curdate == seconddate:
            accs, tavgs = _read_second_daily_file(infile)
        else:
            accs, tavgs = _read_next_daily_file(infile, accs, tavgs)
        curdate += delta

    # Finalize averages (dividing by number of days).  Then write
    # tavgs and accs to file
    tavgs = _finalize_tavgs(tavgs)
    _update_monthly_s2s_values(tmp_outfile, accs, tavgs)
    del accs
    del tavgs

    # Clean up a few details.
    infile = _create_daily_s2s_filename(input_dir, enddate, model_forcing)
    _add_time_data(infile, tmp_outfile, startdate, enddate)
    _update_cell_methods(tmp_outfile)
    _cleanup_global_attrs(tmp_outfile)

    # Rename the output file
    outfile = _create_monthly_s2s_filename(output_dir, startdate,
                                           enddate, model_forcing)
    os.rename(tmp_outfile, outfile)

# Invoke driver
if __name__ == "__main__":
    _driver()
