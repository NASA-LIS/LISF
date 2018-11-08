#!/usr/bin/env python
#------------------------------------------------------------------------------
#
# SCRIPT: convert_nc2surf.py
# 
# PURPOSE:  Reads netCDF4 files from LDT and LVT containing JULES-related
# variables on GALWEM global grid, and convert to SURF format for input
# by GALWEM.
#
# REQUIREMENTS:  
# * Python 2.6 or 2.7
# * UKMO MULE Python library (for writing SURF files)
# * UNIDATA NetCDF4 Python library (for reading netCDF4 files)
# * NumPy Python library (for array objects)
#
# REVISION HISTORY:
# 26 Jan 2018:  Eric Kemp (SSAI), first version.
# 01 Feb 2018:  Eric Kemp (SSAI), now reads command line args.
# 21 Mar 2018:  Eric Kemp (SSAI), updated netCDF variable names to remain
#               consistent with LVT.
# 12 Apr 2018:  Eric Kemp (SSAI), bug fix for specifying soil layer
#               thicknesses.
# 13 Apr 2018:  Eric Kemp (SSAI), added unit conversion for soil moisture.
# 16 Apr 2018:  Eric Kemp (SSAI), fixed LBPLEV setting for snow and surface
#               temperature fields.
# 18 Apr 2018:  Eric Kemp (SSAI), disabled LEVEL_DEPENDENT_CONSTANTS section
#               for most SURF files.
#
#------------------------------------------------------------------------------
# Standard modules
import datetime
import decimal
import os
import sys

#------------------------------------------------------------------------------
# Non-standard modules
import mule
import netCDF4 
import numpy as np

#------------------------------------------------------------------------------
# Version of UM model.  Used below when assigning metadata to SURF files.
_MODEL_VERSION = 1006

#------------------------------------------------------------------------------
# SURF Variable codes, based on STASHmaster_A file used with GALWEM.
# Each key is a combination of the variable name and the source program that
# produced the netCDF file.
varIds = {
    # Surface temperature after timestep (1 level)...From LVT
    "water_temp:LVT" : {  
        "LBFC" : 16,       # PPF
        "ITEM_CODE" : 24,  # Item
        "LBTYP" : 18,      # CFFF
        "LBVC" : 129,      # Value at surface
        "LBPLEV" : 0,      # STASH pseudo dimension
        },
    # SST anomaly (1 level)...From LVT
    "SSTanom:LVT" : { # FIXME:  Not supported yet
        "LBFC" : 0,        # PPF
        "ITEM_CODE" : 39,  # Item 
        "LBTYP" : 0,       # CFFF
        "LBVC" : 129,      # Value at surface
        "LBPLEV" : 0,      # STASH pseudo dimension
        },
    # Frac of sea ice in sea after tstep (1 level)...From LVT
    "aice:LVT" : { 
        "LBFC" :      37, # PPF
        "ITEM_CODE" : 31, # Item
        "LBTYP" : 134,    # CFFF
        "LBVC" : 129,     # Value at surface
        "LBPLEV" : 0,     # STASH pseudo dimension
        },
    # Sea ice depth (mean over ice) (1 level)...From LVT
    "hi:LVT" : { 
        "LBFC" : 687,     # PPF 
        "ITEM_CODE" : 32, # Item 
        "LBTYP" : 0,      # CFFF
        "LBVC" : 129,     # Value at surface
        "LBPLEV" : 0,     # STASH pseudo dimension
        },
    # Snow amount over land aft tstp (1 level)...From LVT
    "SWE_inst:LVT" : { 
        "LBFC" :      93, # PPF
        "ITEM_CODE" : 23, # Item
        "LBTYP" : 121,    # CFFF
        "LBVC" : 129,     # Value at surface
        "LBPLEV" : 1,     # STASH pseudo dimension
        },
    # Snow depth on ground on tiles (1 level)...From LVT
    "SnowDepth_inst:LVT" : {
        "LBFC" : 0,        # PPF
        "ITEM_CODE" : 376, # Item
        "LBTYP" : 0,       # CFFF
        "LBVC" : 129,      # Value at surface
        "LBPLEV" : 1,      # STASH pseudo dimension
        },
    # Soil moisture content in a layer (4 levels)...From LVT
    "SoilMoist_inst:LVT" : {
        "LBFC" : 122,      # PPF
        "ITEM_CODE" : 9,   # Item
        "LBTYP" : 191,     # CFFF
        "LBVC" : 6,        # Based on UKMO sample
        "LBPLEV" : 0,      # STASH pseudo dimension
        },
    # Deep soil temp after timestep (4 levels)...From LVT
    "SoilTemp_inst:LVT" : {
        "LBFC" :      23, # PPF   
        "ITEM_CODE" : 20, # Item
        "LBTYP" : 190,    # CFFF
        "LBVC" : 6,       # Based on UKMO sample
        "LBPLEV" : 0,     # STASH pseudo dimension
        },
    # Surface temperature on tiles (1 level)...From LVT
    "AvgSurfT_inst:LVT" : {
        "LBFC" : 1510,     # PPF 
        "ITEM_CODE" : 233, # Item
        "LBTYP" : 0,       # CFFF
        "LBVC" : 129,      # Value at surface
        "LBPLEV" : 1,      # STASH pseudo dimension
        },
    # Vol SMC at Wilting after timestep...From LDT
    "JULES_SM_WILT:LDT" : {
        "LBFC" : 329,     # PPF
        "ITEM_CODE" : 40, # Item 
        "LBTYP" : 491,    # CFFF
        "LBVC" : 129,     # Value at surface
        "LBPLEV" : 0,     # STASH pseudo dimension
        },
    # Vol SMC at Crit Pt after timestep...From LDT
    "JULES_SM_CRIT:LDT" : {
        "LBFC" : 330,     # PPF
        "ITEM_CODE" : 41, # Item
        "LBTYP" : 492,    # CFFF
        "LBVC" : 129,     # Value at surface
        "LBPLEV" : 0,     # STASH pseudo dimension
        },
    }

#------------------------------------------------------------------------------
# Class for generating the SURF files.
class nc2surf(object):

    #--------------------------------------------------------------------------
    # Initialize object
    def __init__(self,ncfile_lvt,ncfile_ldt):
        # Open LDT netCDF file
        self.ncid_ldt = netCDF4.Dataset(ncfile_ldt,mode='r', 
                                        format='NETCDF4_CLASSIC')
        
        # Open LVT netCDF file
        self.ncid_lvt = netCDF4.Dataset(ncfile_lvt,mode='r',
                                        format='NETCDF4_CLASSIC')
        # Fetch valid time
        time = self.ncid_lvt.variables["time"]
        yyyymmdd = time.__dict__["begin_date"]
        self.year = int(yyyymmdd[0:4])
        self.month = int(yyyymmdd[4:6])
        self.day = int(yyyymmdd[6:])
        hhmmss = time.__dict__["begin_time"]
        self.hour = int(hhmmss[0:2])
        self.minute = int(hhmmss[2:4])
        self.second = int(hhmmss[4:])
        del time
        # Fetch constants
        self.ncols = len(self.ncid_lvt.dimensions['east_west'])
        self.nrows = len(self.ncid_lvt.dimensions['north_south'])
        self.nlevs = self.ncid_lvt.__dict__['NUM_SOIL_LAYERS']
        self.dx = self.ncid_lvt.__dict__['DX']
        self.dy = self.ncid_lvt.__dict__['DY']
        self.start_lat = self.ncid_lvt.variables['latitude'][0,0]

        # Convert soil layer depths to meters.  Use fixed precision.
        slt_cm = \
            self.ncid_lvt.__dict__['SOIL_LAYER_THICKNESSES']               
        cm2m = decimal.Decimal('0.01')
        slt_m = \
            [cm2m*decimal.Decimal("%s" %x) for x in slt_cm]
        self.soil_layer_thicknesses = slt_m

        # Find soil layer depths
        self.soil_layer_depths = self.soil_layer_thicknesses[:]
        for i in range(1,len(self.soil_layer_thicknesses)):
            self.soil_layer_depths[i] += self.soil_layer_depths[i-1]

        # Determine first non-negative longitude and index
        lons = self.ncid_lvt.variables['longitude'][0,:]
        self.start_lon = None
        for i in range(0,len(lons)):
            if not lons[i] < 0:
                self.start_lon = lons[i]
                break
        if self.start_lon == None:
            print "ERROR finding starting longitude for LIS grid!"
            sys.exit(1)
        self.i_pm = i # First index at or past prime meridian
        del lons

    #--------------------------------------------------------------------------
    # Internal method for creating SURF object from template.  Refer to
    # Unified Model Documentation Paper F03 for meaning of metadata codes.
    def _create_surf_template(self,file_type,varfields):

        # For SST field, currently only 1 field will be written
        # (Surface Temperature After Timestep).  Longer term, we will
        # add SST Anomaly.
        _glu_sst_template = {
            'fixed_length_header' : {
                'sub_model' : 1,
                'vert_coord_type' : 5,
                'horiz_grid_type' : 0,
                'dataset_type' : 4,
                'run_identifier' : 0,
                'calendar' :   1,
                'grid_staggering' :   6,
                'time_type' :   0,
                'model_version' :   _MODEL_VERSION,
                },
            'integer_constants' : {
                'num_times' :   1,
                'num_levels' :   1, 
                'num_field_types' :   1, # Change to 2 when we add SST anom
                },
            'real_constants' : {
                'north_pole_lat' :   90.0,
                'north_pole_lon' :   0.0,
                },
            'level_dependent_constants' : {
                'dims' : (1,None),  # MULE wants this
                },
            }

        _glu_ice_template = {
            "fixed_length_header" : {
                'sub_model' :   1,
                'vert_coord_type' :   5,
                'horiz_grid_type' :   0,
                'dataset_type' :   4,
                'run_identifier' :   0,
                'calendar' :   1,
                'grid_staggering' :   6,
                'time_type' :   0,
                'model_version' :   _MODEL_VERSION,
                },
            'integer_constants' : {
                'num_times' :   1,
                'num_levels' :   1, 
                'num_field_types' :   2, # Two variables in ice file
                },
            'real_constants' : {
                'north_pole_lat' :   90.0,
                'north_pole_lon' :   0.0,
                },
            'level_dependent_constants' : {
                'dims' : (1,None),  # MULE wants this
                },
            }

        _glu_snow_template = {
            'fixed_length_header' : {
                'sub_model' :   1,
                'vert_coord_type' :   1,
                'horiz_grid_type' :   0,
                'dataset_type' :   4,
                'run_identifier' :   0,
                'calendar' :   1,
                'grid_staggering' :   6,
                'time_type' :   0,
                'model_version' :  _MODEL_VERSION,
                },
            'integer_constants' : {
                'num_times' :   1,
                'num_levels' :   1,  
                'num_field_types' :  2, # Two variables in snow file
                },
            'real_constants' : {
                'north_pole_lat' :   90.0,
                'north_pole_lon' :   0.0,
                },
            'level_dependent_constants' : {
                'dims' : (1,None),  # MULE wants this
                },
            }

        _glu_smc_template = {
            'fixed_length_header' : {
                'sub_model' :   1,
                'vert_coord_type' :   1,
                'horiz_grid_type' :   0,
                'dataset_type' :   4,
                'run_identifier' :   0,
                'calendar'   : 1,
                'grid_staggering' :   6,
                'time_type' : 0,
                'model_version' : _MODEL_VERSION,
                },
            'integer_constants' : {
                'num_times' :   1,
#                'num_levels' :   1,  # In sample _glu_smc file
                'num_levels' :   4,   # MULE wants this instead
                'num_field_types' :   5, # Five variables in smc file
                },
            'real_constants' : {
                'north_pole_lat' :   90.0,
                'north_pole_lon' :   0.0,
                },
            'level_dependent_constants' : {
                'dims' : (4,None),  # MULE wants this for 4 soil layers
                },
            }

        _templates = {
            '_glu_sst' : _glu_sst_template,
            '_glu_ice' : _glu_ice_template,
            '_glu_snow' : _glu_snow_template,
            '_glu_smc' : _glu_smc_template,
            }

        # Select appropriate template
        self.template = _templates[file_type]

        # Customize appropriate settings
        self.template["fixed_length_header"]["t1_year"] = self.year
        self.template["fixed_length_header"]["t1_month"] = self.month
        self.template["fixed_length_header"]["t1_day"] = self.day
        self.template["fixed_length_header"]["t1_hour"] = self.hour
        self.template["fixed_length_header"]["t1_minute"] = self.minute
        self.template["fixed_length_header"]["t1_second"] = self.second

        self.template["fixed_length_header"]["t2_year"] = self.year
        self.template["fixed_length_header"]["t2_month"] = self.month
        self.template["fixed_length_header"]["t2_day"] = self.day
        self.template["fixed_length_header"]["t2_hour"] = self.hour
        self.template["fixed_length_header"]["t2_minute"] = self.minute
        self.template["fixed_length_header"]["t2_second"] = self.second

        self.template["fixed_length_header"]["t3_year"] = self.year
        self.template["fixed_length_header"]["t3_month"] = self.month
        self.template["fixed_length_header"]["t3_day"] = self.day
        self.template["fixed_length_header"]["t3_hour"] = self.hour
        self.template["fixed_length_header"]["t3_minute"] = self.minute
        self.template["fixed_length_header"]["t3_second"] = self.second

        self.template["integer_constants"]["num_cols"] = self.ncols
        self.template["integer_constants"]["num_rows"] = self.nrows
        self.template["integer_constants"]["num_field_types"] = len(varfields)

        self.template["real_constants"]["col_spacing"] = self.dx
        self.template["real_constants"]["row_spacing"] = self.dy
        self.template["real_constants"]["start_lat"] = self.start_lat
        self.template["real_constants"]["start_lon"] = self.start_lon

        if file_type == "_glu_smc":
            self.template["level_dependent_constants"]["soil_thickness"] = \
                self.soil_layer_thicknesses[:]

    #--------------------------------------------------------------------------
    # Internal method to create and attach Field object to SURF object
    # Refer to Unified Model Documentation Paper F03 for meaning of metadata.
    def _add_field(self,key,var2d,lblev,ilev,nlev,var2d_provider,surf):
        # Determine how many fields are already in the SURF object
        num_fields = len(surf.fields)
        # Create Field3 object and populate records
        field = mule.Field3.empty()
        field.lbyr = self.year
        field.lbmon = self.month
        field.lbdat = self.day
        field.lbhr = self.hour
        field.lbmin = self.minute
        field.lbsec = self.second
        field.lbyrd = self.year
        field.lbmond = self.month
        field.lbdatd = self.day
        field.lbhrd = self.hour
        field.lbmind = self.minute
        field.lbsecd = self.second
        field.lbtim = 1 # Dates use Gregorian calendar
        field.lbft = 0 # No difference between valid time and data time
        #field.lblrec = foo # Let MULE determine this
        field.lbcode = 1 # Lat/lon grid
        field.lbhem = 0 # Global field
        # Note:  Set lbrow and lbnpt to zero for land-packed field
        field.lbrow = self.nrows
        field.lbnpt = self.ncols
        field.lbext = 0 # No extra data
        field.lbpack = 0 # No packing
        field.lbrel = 2
        field.lbfc = varIds[key]["LBFC"]
        field.lbcfc = 0 # Always 0 for UM
        field.lbproc = 0 # No processing -- raw field
        field.lbvc = varIds[key]["LBVC"]
        field.lbrvc = 0
        field.lbexp = 0
        #field.lbegin = 0 # For AncilFile???
        #field.lbnrec = foo # Disk length...Let MULE handle this
        #field.lbproj = 900 # Not populated for UM post 8.4
        field.lbproj = mule._INTEGER_MDI
        field.lbtyp = varIds[key]["LBTYP"]
        field.lblev = lblev
        field.lbrsvd1 = 0
        field.lbrsvd2 = 0
        field.lbrsvd3 = 0
        field.lbrsvd4 = 0
        field.lbsrce = 10000*_MODEL_VERSION + 1111 
        # FIXME...Set to 3 for landsea mask
        field.lbuser1 = 1 # Real field
        field.lbuser2 = (num_fields*self.nrows*self.ncols) + 1
        field.lbuser3 = 0 # No rim or halo sizes
        field.lbuser4 = varIds[key]["ITEM_CODE"]
        #field.lbuser5 = 0
        field.lbuser5 = varIds[key]["LBPLEV"] # "STASH Psuedo dimension"
        field.lbuser6 = 0 # Free space for users...Let MULE handle this
        field.lbuser7 = 1 # Atmosphere
        field.bdatum = 0 # Datum value constant subtracted from each field 
                         # value
        if nlev == 1:
            field.blev = 0 # Surface AGL
        else:
            #field.blev = self.soil_layer_depths[ilev]
            field.blev = self.soil_layer_thicknesses[ilev]
        field.bplat = 90
        field.bplon = 0
        # Grid settings
        field.bgor = 0
        field.bzy = self.start_lat - self.dy
        field.bdy = self.dy
        field.bzx = self.start_lon - self.dx
        field.bdx = self.dx
        field.bmdi = mule._REAL_MDI
        field.bmks = 1.0
        field.raw[1] = self.year
        field.set_data_provider(var2d_provider)
        surf.fields.append(field)
        return surf

    #--------------------------------------------------------------------------
    # Method for creating SURF object with fields
    def create_surf_file(self,file_type,varlist,surffile):

        # Create template for SURF file
        self._create_surf_template(file_type,varlist)

        # Create SURF object
        surf = mule.AncilFile.from_template(self.template)

        # The UM RECON preprocessor does not work if the SURF file (other
        # than _glu_smc) has a LEVEL_DEPENDENT_CONSTANTS section.  However,
        # the MULE library executes validation code that requires the 
        # LEVEL_DEPENDENT_CONSTANTS section to exist.  The current 
        # workaround is to disable the validation for this specific SURF 
        # file, and remove the troublesome section.
        if file_type in ["_glu_snow","_glu_ice","_glu_sst"]:
            def dummy_validate(*args, **kwargs):
                pass
            surf.validate = dummy_validate
            surf.level_dependent_constants = None

        # Create Field3 object for each variable
        for key in varlist:

            # See if the varname and source are recognized
            if key not in varIds.keys():
                print "WARN, %s not recognized!" %(key)
                continue

            # See if the source is recognized
            infile_type = key.split(":")[1]
            if infile_type not in ["LDT","LVT"]:
                print "ERROR, invalid infile type %s" %(infile_type)
                print "Found in %s" %(key)
                print "Internal error, aborting..."
                sys.exit(1)
            
            # Trim the varname to exclude the source
            varid = key.split(":")[0]

            # Attempt to retrieve the variable from the appropriate
            # netCDF file.
            if infile_type == "LDT":
                try:
                    var = self.ncid_ldt.variables[varid]
                except:
                    print "WARN, %s not available in LDT file!" %(varid)
                    continue
            elif infile_type == "LVT":
                try:
                    var = self.ncid_lvt.variables[varid]
                except:
                    print "WARN, %s not available in LVT file!" %(varid)
                    continue

            # Save the "missing data" value for this netCDF variable
            if infile_type == "LVT":
                fillValue = var._FillValue
            elif infile_type == "LDT":
                fillValue = var.missing_value

            # At this point we have a reference to the variable.  Copy to
            # a NumPy array, and record the number of vertical levels.
            if var.ndim == 2:
                var = var[:,:]
                nlev = 1
            elif var.ndim == 3:
                var = var[:,:,:]
                nlev = var.shape[0]
            else:
                print "ERROR, unsupported array with ",ndim,' dimensions!'
                sys.exit(1)

            # Loop through each level.
            for ilev in range(0,nlev):

                # In 2D case, work with the whole array.
                if var.ndim == 2:
                    var2d = var[:,:]
                    lblev = 9999 # Indicates surface level
                # In 3D case, pull out the current vertical level as
                # a 2D array.
                else:
                    var2d = var[ilev,:,:]
                    lblev = ilev+1 # Use 1-based indexing

                # MULE doesn't like masked arrays, so pull the raw
                # data out in this case.
                if type(var2d) == np.ma.core.MaskedArray:
                    var2d = var2d.data

                # Update the missing value to match that used in SURF
                var2d = np.where(var2d == fillValue, mule._REAL_MDI, var2d)

                # EMK...For SoilMoist, convert from m3/m3 to kg m-2.
                if varid == "SoilMoist_inst":
                    soil_layer_thicknesses = \
                        self.ncid_lvt.getncattr("SOIL_LAYER_THICKNESSES")
                    dzsoil = soil_layer_thicknesses[ilev]*0.01 # cm to m
                    var2d = np.where(var2d == mule._REAL_MDI,
                                     mule._REAL_MDI,
                                     var2d*1000*dzsoil)

                # Rotate the field to match the 0:360 longitudinal convention
                # used by GALWEM.  Then create a "provider" of the data.
                var2d_for_surf = np.roll(var2d,self.i_pm,axis=1)
                var2d_provider = mule.ArrayDataProvider(var2d_for_surf)

                # Now add the field to the SURF object.
                print "var: %s, ilev: %s" %(key,ilev)
                surf = self._add_field(key,var2d_for_surf,lblev,ilev, \
                                           nlev,var2d_provider,surf)
                
        # All fields have been added to the SURF object.  Write to file.
        surf.to_file(surffile)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Print command line usage.
def usage():
    print "Usage: %s yyyymmddhh lvt_nc ldt_nc" %(sys.argv[0])
    print "   where:"
    print "           yyyymmddhh is valid year/month/day/hour in UTC"
    print "           lvt_nc is name of LVT netCDF file to convert to SURF"
    print "           ldt_nc is name of LDT netCDF file with JULES ancillaries"

#------------------------------------------------------------------------------
# Read command line arguments
def read_cmd_args():

    # Check if argument count is correct
    if len(sys.argv) != 4:
        print "[ERR] Invalid number of command line arguments!"
        usage()
        sys.exit(1)

    # Convert yyyymmddhh argument to a datetime object
    yyyymmddhh = sys.argv[1]
    try:
        year   = int(yyyymmddhh[0:4])
        month  = int(yyyymmddhh[4:6])
        day    = int(yyyymmddhh[6:8])
        hour   = int(yyyymmddhh[8:10])
        validdt = datetime.datetime(year,month,day,hour)
    except:
        print "[ERR] Cannot process yyyymmddhh argument!"
        usage()
        sys.exit(1)

    # See if lvt_nc file can be opened
    lvt_nc = sys.argv[2]
    try:
        ncid_lvt = netCDF4.Dataset(lvt_nc,mode='r',
                                   format='NETCDF4_CLASSIC')
    except:
        print "[ERR] Cannot open %s with netCDF4 library!" %(lvt_nc)
        usage()
        sys.exit(1)
    ncid_lvt.close()

    # See if ldt_nc file can be opened
    ldt_nc = sys.argv[3]
    try:
        ncid_ldt = netCDF4.Dataset(ldt_nc,mode='r',
                                   format='NETCDF4_CLASSIC')
    except:
        print "[ERR] Cannot open %s with netCDF4 library!" %(ldt_nc)
        usage()
        sys.exit(1)
    ncid_ldt.close()
    
    return validdt,lvt_nc,ldt_nc

#------------------------------------------------------------------------------
# Main Driver.

if __name__ == "__main__":

    # Process command line arguments
    validdt,lvt_nc,ldt_nc = read_cmd_args()

    # Create SURF object
    surf = nc2surf(lvt_nc,ldt_nc)

    # Generate glu_sst SURF file
    file_type = "_glu_sst"
    # FIXME...Add support for SST anomaly
    varfields = ["water_temp:LVT"]
    surffile = "%4.4d%2.2d%2.2dT%2.2d00Z%s" \
        %(validdt.year,
          validdt.month,
          validdt.day,
          validdt.hour,
          file_type)
    surf.create_surf_file(file_type,varfields,surffile)

    # Generate glu_ice SURF file
    file_type = "_glu_ice"
    varfields = ["aice:LVT","hi:LVT"]
    surffile = "%4.4d%2.2d%2.2dT%2.2d00Z%s" \
        %(validdt.year,
          validdt.month,
          validdt.day,
          validdt.hour,
          file_type)
    surf.create_surf_file(file_type,varfields,surffile)

    # Generate glu_snow SURF file
    file_type = "_glu_snow"
    varfields = ["SWE_inst:LVT","SnowDepth_inst:LVT"]
    surffile = "%4.4d%2.2d%2.2dT%2.2d00Z%s" \
        %(validdt.year,
          validdt.month,
          validdt.day,
          validdt.hour,
          file_type)
    surf.create_surf_file(file_type,varfields,surffile)

    # Generate glu_smc SURF file
    file_type = "_glu_smc"
    varfields = ["SoilMoist_inst:LVT","SoilTemp_inst:LVT", "AvgSurfT_inst:LVT",
                 "JULES_SM_WILT:LDT","JULES_SM_CRIT:LDT"]
    surffile = "%4.4d%2.2d%2.2dT%2.2d00Z%s" \
        %(validdt.year,
          validdt.month,
          validdt.day,
          validdt.hour,
          file_type)
    surf.create_surf_file(file_type,varfields,surffile)

    print "convert_nc2surf.py completed!"
    sys.exit(0)


