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
# 27 Oct 2021: Eric Kemp/SSAI, addressed pylint string format complaints.
# 29 Oct 2021: Eric Kemp/SSAI, added config file.
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
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

# Private methods.
def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} configfile input_dir output_dir"
    txt += " start_yyyymmdd end_yyyymmdd model_forcing"
    print(txt)
    print("[INFO] where:")
    print("[INFO]  configfile: path to s2spost config file")
    print("[INFO]  input_dir: directory with daily S2S files in CF convention")
    print("[INFO]  output_dir: directory for output file")
    print("[INFO]  start_yyyymmdd: Starting date/time of daily files")
    print("[INFO]  end_yyyymmdd: Ending date/time of daily files")
    print("[INFO]  model_forcing: ID for atmospheric forcing for LIS")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 7:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if config file exists.
    configfile = sys.argv[1]
    if not os.path.exists(configfile):
        print(f"[ERR] Directory {configfile} does not exist!")
        sys.exit(1)

    # Check if input directory exists.
    input_dir = sys.argv[2]
    if not os.path.exists(input_dir):
        print(f"[ERR] Directory {input_dir} does not exist!")
        sys.exit(1)

    # Create output directory if it doesn't exist.
    output_dir = sys.argv[3]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get valid starting and ending dates of data.
    start_yyyymmdd = sys.argv[4]
    startdate = _proc_date(start_yyyymmdd)
    end_yyyymmdd = sys.argv[5]
    enddate = _proc_date(end_yyyymmdd)
    if startdate > enddate:
        print("[ERR] Start date is after end date!")
        sys.exit(1)

    # Get ID for model forcing for LIS
    model_forcing = sys.argv[6]

    return configfile, input_dir, output_dir, startdate, enddate, model_forcing

def _read_config(configfile):
    """Read from s2spost config file."""
    config = configparser.ConfigParser()
    config.read(configfile)
    return config

def _make_varlists(config):
    """Build lists of variables."""
    varlists = {}
    varlists["var_acc_list"] = config["s2spost"]["var_acc_list"].split()
    varlists["var_tavg_land_list"] = \
        config["s2spost"]["var_tavg_land_list"].split()
    varlists["var_tavg_f_list"] = \
        config["s2spost"]["var_tavg_f_list"].split()
    varlists["var_tavg_twsgws_list"] = \
        config["s2spost"]["var_tavg_twsgws_list"].split()
    varlists["var_tair_max_list"] = \
        config["s2spost"]["var_tair_max_list"].split()
    varlists["var_tair_min_list"] = \
        config["s2spost"]["var_tair_min_list"].split()
    varlists["var_inst_list"] = config["s2spost"]["var_inst_list"].split()
    varlists["const_list"] = config["s2spost"]["const_list"].split()
    return varlists

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

def _check_filename_size(name):
    """Make sure filename does not exceed 128 characters, per Air Force
    requirement."""
    if len(os.path.basename(name)) > 128:
        print("[ERR] Output file name is too long!")
        print(f"[ERR] {os.path.basename(name)} exceeds 128 characters!")
        sys.exit(1)

def _create_daily_s2s_filename(input_dir, curdate, model_forcing):
    """Create path to daily S2S netCDF file."""
    name = f"{input_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing}"
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S"
    name += f"_DD.{curdate.year:04d}{curdate.month:02d}{curdate.day:02d}"
    name += "_DT.0000"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

def _create_monthly_s2s_filename(output_dir, startdate, enddate,
                                 model_forcing):
    """Create path to monthly S2S netCDF file."""
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing}"
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S"
    name += f"_DP.{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
    name += f"-{enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
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

def _create_firstguess_monthly_file(varlists, infile, outfile):
    """Read daily S2S file, and most fields copy to monthly S2S file.
    This allows us to cleanly copy dimensions and all attributes.  The
    numerical values of the arrays in the monthly S2S
    file will be replaced later in the script."""

    if not os.path.exists(infile):
        print(f"[ERR] {infile} does not exist!")
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    # Create monthly file, copying dimensions and global attributes.
    # NOTE:  We clean up the global attributes later in the script.
    ncid_out = nc4_dataset(outfile, "w", format='NETCDF4_CLASSIC')
    _copy_dims_gattrs(ncid_in, ncid_out)

    # Copy the fields.
    varnames = []
    for listname in ["const_list", "var_inst_list", "var_acc_list",
                     "var_tavg_land_list", "var_tavg_f_list",
                     "var_tavg_twsgws_list", "var_tair_max_list",
                     "var_tair_min_list"]:
        varnames += varlists[listname]
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

def _read_second_daily_file(varlists, infile):
    """Read the second daily S2S file and copy the acc and tavg fields in
    appropriate dictionaries. We use the second file to start, since acc
    and tavg are valid for the prior 24-hr period."""

    if not os.path.exists(infile):
        print(f"[ERR] {infile} does not exist!")
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    accs = {}
    tavgs = {}
    tavgs["counter"] = 1

    # Copy the values of the fields we will average or accumulate
    varnames = []
    for listname in ["var_tavg_land_list", "var_tavg_f_list",
                     "var_tavg_twsgws_list", "var_tair_max_list",
                     "var_tair_min_list"]:
        varnames += varlists[listname]
    for varname in varnames:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            tavgs[varname] = var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            tavgs[varname] = var_in[:,:,:]
        elif len(var_in.shape) == 2:
            tavgs[varname] = var_in[:,:]
    for varname in varlists["var_acc_list"]:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            accs[varname] = var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            accs[varname] = var_in[:,:,:]
        elif len(var_in.shape) == 2:
            accs[varname] = var_in[:,:]

    ncid_in.close()
    return accs, tavgs

def _read_next_daily_file(varlists, infile, accs, tavgs):
    """Read next daily S2S file and copy the required variable values to
    appropriate dictionaries."""

    if not os.path.exists(infile):
        print(f"[ERR] {infile} does not exist!")
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    tavgs["counter"] += 1

    # Add the values of the fields we will average or accumulate
    varnames = []
    for listname in ["var_tavg_land_list", "var_tavg_f_list",
                     "var_tavg_twsgws_list", "var_tair_max_list",
                     "var_tair_min_list"]:
        varnames += varlists[listname]
    for varname in varnames:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            tavgs[varname][:,:,:,:] += var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            tavgs[varname][:,:,:] += var_in[:,:,:]
        elif len(var_in.shape) == 2:
            tavgs[varname][:,:] += var_in[:,:]
    for varname in varlists["var_acc_list"]:
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
        print(f"[ERR] {outfile} does not exist!")
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
        print(f"[ERR] {infile} does not exist!")
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')
    if not os.path.exists(outfile):
        print(f"[ERR] {outfile} does not exist!")
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

def _update_cell_methods(varlists, outfile):
    """Update cell_method attributes for select variables."""

    if not os.path.exists(outfile):
        print(f"[ERR] {outfile} does not exist!")
        sys.exit(1)
    ncid = nc4_dataset(outfile, 'a', format='NETCDF4_CLASSIC')

    # Elaborate on monthly calculations of most variables
    varnames = []
    for listname in ["var_acc_list", "var_tavg_land_list", "var_tavg_f_list",
                     "var_tavg_twsgws_list", "var_tair_max_list",
                     "var_tair_min_list"]:
        varnames += varlists[listname]
    for varname in varnames:

        var = ncid.variables[varname]

        # Special handling for TWS_inst and GWS_inst -- we want the monthly
        # means.
        if varname in varlists["var_tavg_twsgws_list"]:
            var.cell_methods = \
                "time: mean (interval: 1 day) area: point where land"

        # Special handling for Tair_f_max and Tair_f_min:  We want the monthly
        # averages of the daily maxs and mins
        elif varname in varlists["var_tair_max_list"]:
            var.cell_methods = \
                "time: mean (interval: 1 day comment: daily maxima)"
        elif varname in varlists["var_tair_min_list"]:
            var.cell_methods = \
                "time: mean (interval: 1 day comment: daily minima)"

        # Clarify monthly mean of daily means over land
        elif varname in varlists["var_tavg_land_list"]:
            var.cell_methods = \
                "time: mean (interval: 1 day comment: daily means)" + \
                " area: point where land"

        # Clarify monthly mean of daily means everywhere
        elif varname in varlists["var_tavg_f_list"]:
            var.cell_methods = \
                "time: mean (interval: 1 day comment: daily means)"

        # Clarify monthly accumulations
        elif varname in varlists["var_acc_list"]:
            var.cell_methods = \
                "time: sum (interval: 1 day comment: daily sums)"

    ncid.close()

def _cleanup_global_attrs(outfile):
    """Clean-up global attributes."""
    if not os.path.exists(outfile):
        print(f"[ERR] {outfile} does not exist!")
        sys.exit(1)
    ncid = nc4_dataset(outfile, 'a', format='NETCDF4_CLASSIC')
    ncid.history = f"created on date: {time.ctime()}"
    del ncid.NCO
    del ncid.history_of_appended_files

def _driver():
    """Main driver."""

    # Get the directories and dates
    configfile, input_dir, output_dir, startdate, enddate, model_forcing \
        = _read_cmd_args()
    config = _read_config(configfile)
    varlists = _make_varlists(config)

    # Loop through dates
    curdate = startdate
    seconddate = startdate + datetime.timedelta(days=1)
    while curdate <= enddate:
        infile = _create_daily_s2s_filename(input_dir, curdate, model_forcing)
        #print("[INFO] Reading %s" %(infile))
        if curdate == startdate:
            tmp_outfile = f"{output_dir}/tmp_monthly.nc"
            _create_firstguess_monthly_file(varlists, infile, tmp_outfile)
        elif curdate == seconddate:
            accs, tavgs = _read_second_daily_file(varlists, infile)
        else:
            accs, tavgs = _read_next_daily_file(varlists, infile, accs, tavgs)
        curdate += datetime.timedelta(days=1)

    # Finalize averages (dividing by number of days).  Then write
    # tavgs and accs to file
    tavgs = _finalize_tavgs(tavgs)
    _update_monthly_s2s_values(tmp_outfile, accs, tavgs)
    del accs
    del tavgs

    # Clean up a few details.
    infile = _create_daily_s2s_filename(input_dir, enddate, model_forcing)
    _add_time_data(infile, tmp_outfile, startdate, enddate)
    _update_cell_methods(varlists, tmp_outfile)
    _cleanup_global_attrs(tmp_outfile)

    # Rename the output file
    outfile = _create_monthly_s2s_filename(output_dir, startdate,
                                           enddate, model_forcing)
    os.rename(tmp_outfile, outfile)

# Invoke driver
if __name__ == "__main__":
    _driver()
