#!/usr/bin/env python
"""
#------------------------------------------------------------------------------
#
# SCRIPT: convert_nc2surf.py
#
# PURPOSE:  Reads netCDF4 files from LDT and LVT containing JULES-related
# variables on GALWEM global grid, and convert to SURF format for input
# by GALWEM.
#
# REQUIREMENTS as of 26 March 2020:
# * Python 3.6
# * UKMO MULE Python library 2020.01.1
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
# 19 Jul 2019:  Eric Kemp (SSAI), adjust soil moisture so it is no less than
#               0.1 * the wilting point.  Also, reset any negative snow fields
#               to zero.
# 31 Mar 2020:  Jim Geiger, Eric Kemp:  Upgraded to Python 3.6.  Also,
#               upgraded to MULE 2020.01.1.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import copy
import datetime
import decimal
import sys

#------------------------------------------------------------------------------
# Non-standard modules
import numpy as np
import netCDF4
import mule

#------------------------------------------------------------------------------
# Version of UM model.  Used below when assigning metadata to SURF files.
_MODEL_VERSION = 1006

#------------------------------------------------------------------------------
# SURF Variable codes, based on STASHmaster_A file used with GALWEM.
# Each key is a combination of the variable name and the source program that
# produced the netCDF file.
_VARIDS = {
    # Surface temperature after timestep (1 level)...From LVT
    "water_temp:LVT": {
        "LBFC": 16,       # PPF
        "ITEM_CODE": 24,  # Item
        "LBTYP": 18,      # CFFF
        "LBVC": 129,      # Value at surface
        "LBPLEV": 0,      # STASH pseudo dimension
    },
    # SST anomaly (1 level)...From LVT
    "SSTanom:LVT": {  # FUTURE:  Not supported yet
        "LBFC": 0,        # PPF
        "ITEM_CODE": 39,  # Item
        "LBTYP": 0,       # CFFF
        "LBVC": 129,      # Value at surface
        "LBPLEV": 0,      # STASH pseudo dimension
    },
    # Frac of sea ice in sea after tstep (1 level)...From LVT
    "aice:LVT": {
        "LBFC":      37,  # PPF
        "ITEM_CODE": 31,  # Item
        "LBTYP": 134,    # CFFF
        "LBVC": 129,     # Value at surface
        "LBPLEV": 0,     # STASH pseudo dimension
    },
    # Sea ice depth (mean over ice) (1 level)...From LVT
    "hi:LVT": {
        "LBFC": 687,     # PPF
        "ITEM_CODE": 32,  # Item
        "LBTYP": 0,      # CFFF
        "LBVC": 129,     # Value at surface
        "LBPLEV": 0,     # STASH pseudo dimension
    },
    # Snow amount over land aft tstp (1 level)...From LVT
    "SWE_inst:LVT": {
        "LBFC":      93,  # PPF
        "ITEM_CODE": 23,  # Item
        "LBTYP": 121,    # CFFF
        "LBVC": 129,     # Value at surface
        "LBPLEV": 1,     # STASH pseudo dimension
    },
    # Snow depth on ground on tiles (1 level)...From LVT
    "SnowDepth_inst:LVT": {
        "LBFC": 0,        # PPF
        "ITEM_CODE": 376,  # Item
        "LBTYP": 0,       # CFFF
        "LBVC": 129,      # Value at surface
        "LBPLEV": 1,      # STASH pseudo dimension
    },
    # Soil moisture content in a layer (4 levels)...From LVT
    "SoilMoist_inst:LVT": {
        "LBFC": 122,      # PPF
        "ITEM_CODE": 9,   # Item
        "LBTYP": 191,     # CFFF
        "LBVC": 6,        # Based on UKMO sample
        "LBPLEV": 0,      # STASH pseudo dimension
    },
    # Deep soil temp after timestep (4 levels)...From LVT
    "SoilTemp_inst:LVT": {
        "LBFC":      23,  # PPF
        "ITEM_CODE": 20,  # Item
        "LBTYP": 190,    # CFFF
        "LBVC": 6,       # Based on UKMO sample
        "LBPLEV": 0,     # STASH pseudo dimension
    },
    # Surface temperature on tiles (1 level)...From LVT
    "AvgSurfT_inst:LVT": {
        "LBFC": 1510,     # PPF
        "ITEM_CODE": 233,  # Item
        "LBTYP": 0,       # CFFF
        "LBVC": 129,      # Value at surface
        "LBPLEV": 1,      # STASH pseudo dimension
    },
    # Vol SMC at Wilting after timestep...From LDT
    "JULES_SM_WILT:LDT": {
        "LBFC": 329,     # PPF
        "ITEM_CODE": 40,  # Item
        "LBTYP": 491,    # CFFF
        "LBVC": 129,     # Value at surface
        "LBPLEV": 0,     # STASH pseudo dimension
    },
    # Vol SMC at Crit Pt after timestep...From LDT
    "JULES_SM_CRIT:LDT": {
        "LBFC": 330,     # PPF
        "ITEM_CODE": 41,  # Item
        "LBTYP": 492,    # CFFF
        "LBVC": 129,     # Value at surface
        "LBPLEV": 0,     # STASH pseudo dimension
    },
}

# Number of levels per output variable.
_NLEVS = {
    "water_temp" : 1,
    "aice" : 1,
    "hi" : 1,
    "SWE_inst" : 1,
    "SnowDepth_inst" : 1,
    "SoilMoist_inst" : 4,
    "SoilTemp_inst" : 4,
    "AvgSurfT_inst" : 1,
    "JULES_SM_WILT" : 1,
    "JULES_SM_CRIT" : 1,
}

#------------------------------------------------------------------------------
# Grab the missing data values from the MULE library
_INTEGER_MDI = copy.deepcopy(mule.IntegerConstants.MDI)
_REAL_MDI = copy.deepcopy(mule.RealConstants.MDI)

#------------------------------------------------------------------------------
def _get_infile_type(key):
    """Get the infile_type based on the key. Sanity check included."""
    infile_type = key.split(":")[1]
    if infile_type not in ["LDT", "LVT"]:
        print("[ERR] Invalid infile type %s" % (infile_type))
        print("[ERR] Found in %s" % (key))
        print("[ERR] Internal error, aborting...")
        sys.exit(1)
    return infile_type

#------------------------------------------------------------------------------
def _handle_snowfields_var2d(var2d):
    """Make sure snow fields are nonnegative."""
    var2d_tmp = np.where(var2d < 0, 0, var2d)
    var2d = np.where(var2d == _REAL_MDI,
                     _REAL_MDI,
                     var2d_tmp)
    del var2d_tmp
    return var2d

#------------------------------------------------------------------------------
class Nc2Surf:
    """Class for generating the SURF files."""

    #--------------------------------------------------------------------------
    def __init__(self, ncfile_lvt, ncfile_ldt):
        """Initializes object"""

        # NOTE:  Pylint complains about netCDF4 not having Dataset as a
        # member.  But this is false -- Dataset is a class, and we use
        # the constructor below.  Since this appears to be a bug, we
        # disable the no-member check in this function.
        # pylint: disable=no-member

        # Open LDT netCDF file
        self.ncid_ldt = netCDF4.Dataset(ncfile_ldt, mode='r',
                                        format='NETCDF4_CLASSIC')

        # Open LVT netCDF file
        self.ncid_lvt = netCDF4.Dataset(ncfile_lvt, mode='r',
                                        format='NETCDF4_CLASSIC')

        # pylint: enable=no-member

        # Fetch valid time
        time = self.ncid_lvt.variables["time"]
        yyyymmdd = time.__dict__["begin_date"]
        self.time = {}
        self.time["year"] = int(yyyymmdd[0:4])
        self.time["month"] = int(yyyymmdd[4:6])
        self.time["day"] = int(yyyymmdd[6:])
        hhmmss = time.__dict__["begin_time"]
        self.time["hour"] = int(hhmmss[0:2])
        self.time["minute"] = int(hhmmss[2:4])
        self.time["second"] = int(hhmmss[4:])
        del time

        # Fetch grid constants
        self.grid = {}
        self.grid["ncols"] = len(self.ncid_lvt.dimensions['east_west'])
        self.grid["nrows"] = len(self.ncid_lvt.dimensions['north_south'])
        self.grid["nlevs"] = self.ncid_lvt.__dict__['NUM_SOIL_LAYERS']
        self.grid["dx"] = self.ncid_lvt.__dict__['DX']
        self.grid["dy"] = self.ncid_lvt.__dict__['DY']
        self.grid["start_lat"] = self.ncid_lvt.variables['latitude'][0, 0]

        # Convert soil layer depths to meters.  Use fixed precision.
        slt_cm = \
            self.ncid_lvt.__dict__['SOIL_LAYER_THICKNESSES']
        cm2m = decimal.Decimal('0.01')
        slt_m = \
            [cm2m*decimal.Decimal("%s" % x) for x in slt_cm]
        # EMK...Use below for old pre-LIS 7.2 LVT output, which
        # erroneously had output in meters
        # slt_m = \
        #    [decimal.Decimal("%s" %x) for x in slt_cm]
        self.grid["soil_layer_thicknesses"] = slt_m

        # Find soil layer depths
        self.grid["soil_layer_depths"] = self.grid["soil_layer_thicknesses"][:]
        for i in range(1, len(self.grid["soil_layer_thicknesses"])):
            self.grid["soil_layer_depths"][i] += \
                self.grid["soil_layer_depths"][i-1]

        # Determine first non-negative longitude and index
        lons = self.ncid_lvt.variables['longitude'][0, :]
        self.grid["start_lon"] = None
        for i, lon in enumerate(lons):
            if not lon < 0:
                self.grid["start_lon"] = lon
                break
        if self.grid["start_lon"] is None:
            print("[ERR] No starting longitude for LIS grid!")
            sys.exit(1)
        self.grid["i_pm"] = i  # First index at or past prime meridian
        del lons

        self.template = None
        self.lblev = None

    #--------------------------------------------------------------------------
    def _create_surf_template(self, file_type, varfields):
        """
        Internal method for creating SURF object from template. Refer to
        Unified Model Documentation Paper F03 for meaning of metadata codes.
        """

        # For SST field, currently only 1 field will be written
        # (Surface Temperature After Timestep).  Longer term, we will
        # add SST Anomaly.
        _glu_sst_template = {
            'fixed_length_header': {
                'sub_model': 1,
                'vert_coord_type': 5,
                'horiz_grid_type': 0,
                'dataset_type': 4,
                'run_identifier': 0,
                'calendar':   1,
                'grid_staggering':   6,
                'time_type':   0,
                'model_version':   _MODEL_VERSION,
            },
            'integer_constants': {
                'num_times':   1,
                'num_levels':   1,
                'num_field_types':   1,  # Change to 2 when we add SST anom
            },
            'real_constants': {
                'north_pole_lat':   90.0,
                'north_pole_lon':   0.0,
            },
        }

        _glu_ice_template = {
            "fixed_length_header": {
                'sub_model':   1,
                'vert_coord_type':   5,
                'horiz_grid_type':   0,
                'dataset_type':   4,
                'run_identifier':   0,
                'calendar':   1,
                'grid_staggering':   6,
                'time_type':   0,
                'model_version':   _MODEL_VERSION,
            },
            'integer_constants': {
                'num_times':   1,
                'num_levels':   1,
                'num_field_types':   2,  # Two variables in ice file
            },
            'real_constants': {
                'north_pole_lat':   90.0,
                'north_pole_lon':   0.0,
            },
        }

        _glu_snow_template = {
            'fixed_length_header': {
                'sub_model':   1,
                'vert_coord_type':   1,
                'horiz_grid_type':   0,
                'dataset_type':   4,
                'run_identifier':   0,
                'calendar':   1,
                'grid_staggering':   6,
                'time_type':   0,
                'model_version':  _MODEL_VERSION,
            },
            'integer_constants': {
                'num_times':   1,
                'num_levels':   1,
                'num_field_types':  2,  # Two variables in snow file
            },
            'real_constants': {
                'north_pole_lat':   90.0,
                'north_pole_lon':   0.0,
            },
        }

        _glu_smc_template = {
            'fixed_length_header': {
                'sub_model':   1,
                'vert_coord_type':   1,
                'horiz_grid_type':   0,
                'dataset_type':   4,
                'run_identifier':   0,
                'calendar': 1,
                'grid_staggering':   6,
                'time_type': 0,
                'model_version': _MODEL_VERSION,
            },
            'integer_constants': {
                'num_times':   1,
                #'num_levels' :   1,  # In sample _glu_smc file
                'num_levels':   4,   # MULE wants this instead
                'num_field_types':   5,  # Five variables in smc file
            },
            'real_constants': {
                'north_pole_lat':   90.0,
                'north_pole_lon':   0.0,
            },
            'level_dependent_constants': {
                'dims': (4, 1),  # MULE 2020.01.1 wants this.
            },
        }

        _templates = {
            '_glu_sst': _glu_sst_template,
            '_glu_ice': _glu_ice_template,
            '_glu_snow': _glu_snow_template,
            '_glu_smc': _glu_smc_template,
        }

        # Select appropriate template
        self.template = _templates[file_type]

        # Customize appropriate settings
        self.template["fixed_length_header"]["t1_year"] = self.time["year"]
        self.template["fixed_length_header"]["t1_month"] = self.time["month"]
        self.template["fixed_length_header"]["t1_day"] = self.time["day"]
        self.template["fixed_length_header"]["t1_hour"] = self.time["hour"]
        self.template["fixed_length_header"]["t1_minute"] = self.time["minute"]
        self.template["fixed_length_header"]["t1_second"] = self.time["second"]

        self.template["fixed_length_header"]["t2_year"] = self.time["year"]
        self.template["fixed_length_header"]["t2_month"] = self.time["month"]
        self.template["fixed_length_header"]["t2_day"] = self.time["day"]
        self.template["fixed_length_header"]["t2_hour"] = self.time["hour"]
        self.template["fixed_length_header"]["t2_minute"] = self.time["minute"]
        self.template["fixed_length_header"]["t2_second"] = self.time["second"]

        self.template["fixed_length_header"]["t3_year"] = self.time["year"]
        self.template["fixed_length_header"]["t3_month"] = self.time["month"]
        self.template["fixed_length_header"]["t3_day"] = self.time["day"]
        self.template["fixed_length_header"]["t3_hour"] = self.time["hour"]
        self.template["fixed_length_header"]["t3_minute"] = self.time["minute"]
        self.template["fixed_length_header"]["t3_second"] = self.time["second"]

        self.template["integer_constants"]["num_cols"] = self.grid["ncols"]
        self.template["integer_constants"]["num_rows"] = self.grid["nrows"]
        self.template["integer_constants"]["num_field_types"] = len(varfields)

        self.template["real_constants"]["col_spacing"] = self.grid["dx"]
        self.template["real_constants"]["row_spacing"] = self.grid["dy"]
        self.template["real_constants"]["start_lat"] = self.grid["start_lat"]
        self.template["real_constants"]["start_lon"] = self.grid["start_lon"]

    #--------------------------------------------------------------------------
    def _set_field_lb(self, num_fields, key, field):
        """Set the "lb" attributes in the Field object"""
        field.lbyr = self.time["year"]
        field.lbmon = self.time["month"]
        field.lbdat = self.time["day"]
        field.lbhr = self.time["hour"]
        field.lbmin = self.time["minute"]
        field.lbsec = self.time["second"]
        field.lbyrd = self.time["year"]
        field.lbmond = self.time["month"]
        field.lbdatd = self.time["day"]
        field.lbhrd = self.time["hour"]
        field.lbmind = self.time["minute"]
        field.lbsecd = self.time["second"]
        field.lbtim = 1  # Dates use Gregorian calendar
        field.lbft = 0  # No difference between valid time and data time
        # field.lblrec = foo # Let MULE determine this
        field.lbcode = 1  # Lat/lon grid
        field.lbhem = 0  # Global field
        # Note:  Set lbrow and lbnpt to zero for land-packed field
        field.lbrow = self.grid["nrows"]
        field.lbnpt = self.grid["ncols"]
        field.lbext = 0  # No extra data
        field.lbpack = 0  # No packing
        field.lbrel = 3 # For UM 8.1 and upwards
        field.lbfc = _VARIDS[key]["LBFC"]
        field.lbcfc = 0  # Always 0 for UM
        field.lbproc = 0  # No processing -- raw field
        field.lbvc = _VARIDS[key]["LBVC"]
        field.lbrvc = 0
        field.lbexp = 0
        # field.lbegin = 0 # For AncilFile???
        # field.lbnrec = foo # Disk length...Let MULE handle this
        # field.lbproj = 900 # Not populated for UM post 8.4
        field.lbproj = _INTEGER_MDI
        field.lbtyp = _VARIDS[key]["LBTYP"]
        field.lblev = self.lblev
        field.lbrsvd1 = 0
        field.lbrsvd2 = 0
        field.lbrsvd3 = 0
        field.lbrsvd4 = 0
        field.lbsrce = 10000*_MODEL_VERSION + 1111
        # FUTURE...Set to 3 for landsea mask
        field.lbuser1 = 1  # Real field
        field.lbuser2 = (num_fields*self.grid["nrows"]*self.grid["ncols"]) + 1
        field.lbuser3 = 0  # No rim or halo sizes
        field.lbuser4 = _VARIDS[key]["ITEM_CODE"]
        #field.lbuser5 = 0
        field.lbuser5 = _VARIDS[key]["LBPLEV"]  # "STASH Psuedo dimension"
        field.lbuser6 = 0  # Free space for users...Let MULE handle this
        field.lbuser7 = 1  # Atmosphere

    #--------------------------------------------------------------------------
    def _add_field(self, key, ilev, var2d_provider, surf):
        """
        Internal method to create and attach Field object to SURF object.
        Refer to Unified Model Documentation Paper F03 for meaning of metadata.
        """

        nlev = _NLEVS[key.split(":")[0]]

        # Determine how many fields are already in the SURF object
        num_fields = len(surf.fields)

        # Create Field3 object
        field = mule.Field3.empty()

        # Populate the field records starting with "lb"
        self._set_field_lb(num_fields, key, field)

        # Populate remaining records
        field.bdatum = 0  # Datum value constant subtracted from field
        if nlev == 1:
            field.blev = 0  # Surface AGL
        else:
            field.blev = self.grid["soil_layer_thicknesses"][ilev]
        field.bplat = 90
        field.bplon = 0
        # Grid settings
        field.bgor = 0
        field.bzy = self.grid["start_lat"] - self.grid["dy"]
        field.bdy = self.grid["dy"]
        field.bzx = self.grid["start_lon"] - self.grid["dx"]
        field.bdx = self.grid["dx"]
        field.bmdi = _REAL_MDI
        field.bmks = 1.0
        field.raw[1] = self.time["year"]
        field.set_data_provider(var2d_provider)

        # Append to the surf object and return that object
        surf.fields.append(field)
        return surf

    #--------------------------------------------------------------------------
    def _handle_glu_smc(self, file_type, surf):
        """Add soil thickness section to _glu_smc SURF file."""

        # MULE 2020.01.1 does not accept a soil thickness section for
        # Ancils, but it is required by the UM RECON preprocessor when
        # processing soil data.  So, we need to add it here, borrowing
        # from FieldsFile.
        if file_type != "_glu_smc":
            return surf

        # Create empty set of level dependent constants, fill the
        # soil thicknesses, and attach.
        ldc = mule.ff.FF_LevelDependentConstants.empty(4) # 4 soil layers

        # NOTE:  Pylint complains that the FF_LevelDependentConstants instance
        # has no soil_thickness member.  But this demonstrably false.  Since
        # this is a bug, we disable the no-member check here.
        # pylint: disable=no-member
        ldc.soil_thickness[:] = self.grid["soil_layer_thicknesses"][:]
        # pylint: enable=no-member

        surf.level_dependent_constants = ldc
        # This won't pass validation anymore, so we disable it.
        # NOTE: Pylint doesn't like this since *args and **kwargs are not
        # used.  Unfortunately, MULE doesn't have an option to just turn off
        # validation, and arguments must be allowed by the validation function.
        # So, for this particular function, we turn off checking for unused
        # arguments.
        def dummy_validate(*args, **kwargs):
            # pylint: disable=unused-argument
            pass
            # pylint: enable=unused-argument

        surf.validate = dummy_validate
        return surf

    #--------------------------------------------------------------------------
    def _handle_soilmoist_wilt(self):
        # EMK...Pull wilting point data if processing Soil Moisture
        var_wilt = self.ncid_ldt.variables["JULES_SM_WILT"]
        fill_value = var_wilt.missing_value
        var_wilt = var_wilt[:, :]  # Copies to NumPy array
        if isinstance(var_wilt, np.ma.core.MaskedArray):
            var_wilt = var_wilt.data
        # Soil moisture should not be below 0.1 * wilting point
        var_wilt0p1 = np.where(var_wilt == fill_value,
                               _REAL_MDI, 0.1*var_wilt)
        return var_wilt0p1

    #--------------------------------------------------------------------------
    def _get_var_and_fill_value(self, infile_type, varid):
        var = None
        fill_value = None
        if infile_type == "LDT":
            if varid in self.ncid_ldt.variables:
                var = self.ncid_ldt.variables[varid]
            else:
                print("[WARN] %s not available in LDT file!" % (varid))
        elif infile_type == "LVT":
            if varid in self.ncid_lvt.variables:
                var = self.ncid_lvt.variables[varid]
            else:
                print("[WARN] %s not available in LVT file!" % (varid))
        if var is None:
            return var, fill_value # Both are None
        fill_value = var.missing_value
        # At this point we have a reference to the variable.  Copy to
        # a NumPy array, and record the number of vertical levels.
        if var.ndim == 2:
            var = var[:, :]
        elif var.ndim == 3:
            var = var[:, :, :]
        else:
            print("[ERR] Unsupported array with ", var.ndim,
                  ' dimensions!')
            sys.exit(1)

        return var, fill_value

    #--------------------------------------------------------------------------
    def _create_var2d(self, var, ilev, fill_value):
        """Creates a 2d numpy array from the requested variable."""

        # In 2D case, work with the whole array.
        if var.ndim == 2:
            var2d = var[:, :]
            self.lblev = 9999  # Indicates surface level
        # In 3D case, pull out the current vertical level as a 2D array
        else:
            var2d = var[ilev, :, :]
            self.lblev = ilev+1  # Use 1-based indexing

        # MULE doesn't like masked arrays, so pull the raw
        # data out in this case.
        if isinstance(var2d, np.ma.core.MaskedArray):
            var2d = var2d.data

        # Update the missing value to match that used in SURF
        var2d = np.where(var2d == fill_value,
                         _REAL_MDI,
                         var2d)

        return var2d

    #--------------------------------------------------------------------------
    def _handle_soilmoist_var2d(self, ilev, var2d, var_wilt0p1):
        """Adjust SoilMoist to no less than 10% of wilting point, and convert
        from m3/m3 to kg m-2."""
        var2d_tmp = np.where(var2d < var_wilt0p1,
                             var_wilt0p1,
                             var2d)
        var2d = np.where(var2d == _REAL_MDI,
                         _REAL_MDI,
                         var2d_tmp)
        del var2d_tmp
        dzsoil = float(self.grid["soil_layer_thicknesses"][ilev])
        var2d = np.where(var2d == _REAL_MDI,
                         _REAL_MDI,
                         var2d*1000*dzsoil)
        del dzsoil
        return var2d

    #--------------------------------------------------------------------------
    def _create_var2d_provider(self, var2d):
        """Create MULE provider for var2d.  Also, rotate the data to
        UKMO convention."""
        var2d_for_surf = np.roll(var2d, self.grid["i_pm"], axis=1)
        var2d_provider = mule.ArrayDataProvider(var2d_for_surf)
        del var2d_for_surf
        return var2d_provider

    #--------------------------------------------------------------------------
    def create_surf_file(self, file_type, varlist, surffile):
        """Method for creating SURF object with fields"""

        # Create template for SURF file
        self._create_surf_template(file_type, varlist)

        # Create SURF object
        surf = mule.AncilFile.from_template(self.template)

        # Add soil thickness section to _glu_smc SURF file.
        surf = self._handle_glu_smc(file_type, surf)

        # Create Field3 object for each variable
        for key in varlist:

            # See if the varname and source are recognized
            if key not in list(_VARIDS.keys()):
                print("[WARN] %s not recognized!" % (key))
                continue

            # See if the source is recognized
            infile_type = _get_infile_type(key)

            # Trim the varname to exclude the source
            varid = key.split(":")[0]

            # EMK...Pull wilting point data if processing Soil Moisture
            if infile_type == "LVT":
                if varid == "SoilMoist_inst":
                    var_wilt0p1 = self._handle_soilmoist_wilt()

            # Attempt to retrieve the variable from the appropriate
            # netCDF file.
            var, fill_value = self._get_var_and_fill_value(infile_type, varid)
            if var is None:
                continue

            # Loop through each level.
            nlev = _NLEVS[varid]
            for ilev in range(0, nlev):

                var2d = self._create_var2d(var, ilev, fill_value)

                # EMK...For SoilMoist, make sure no less than 0.1*wilting point
                # Wilting point is in m3/m3.  Then convert from m3/m3 to
                # kg m_2.
                if varid == "SoilMoist_inst":
                    var2d = self._handle_soilmoist_var2d(ilev, var2d,
                                                         var_wilt0p1)

                # EMK...For snow fields, make sure nonnegative
                if varid in ["SWE_inst", "SnowDepth_inst"]:
                    var2d = _handle_snowfields_var2d(var2d)

                # Rotate the field to match the 0:360 longitudinal convention
                # used by GALWEM.  Then create a "provider" of the data.
                var2d_provider = self._create_var2d_provider(var2d)

                # Now add the field to the SURF object.
                print("[INFO] Processing %s, ilev: %s" % (key, ilev))

                surf = self._add_field(key, ilev,
                                       var2d_provider, surf)

        # All fields have been added to the SURF object.  Write to file.
        surf.to_file(surffile)

    #--------------------------------------------------------------------------
    def __str__(self):
        return self.__class__.__name__

#------------------------------------------------------------------------------
def usage():
    """Print command line usage."""
    print("[INFO] Usage: %s yyyymmddhh lvt_nc ldt_nc" % (sys.argv[0]))
    print("[INFO]   where:")
    print("[INFO]    yyyymmddhh is valid year/month/day/hour in UTC")
    print("[INFO]    lvt_nc is name of LVT netCDF file to convert to SURF")
    print("[INFO]    ldt_nc is name of LDT netCDF file with JULES ancillaries")

#-----------------------------------------------------------------------------
def read_cmd_args():
    """Read command line arguments"""
    # Check if argument count is correct
    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        usage()
        sys.exit(1)

    # Convert yyyymmddhh argument to a datetime object
    yyyymmddhh = sys.argv[1]
    if len(yyyymmddhh) != 10:
        print("[ERR] Invalid yyyymmddhh argument %s!" %(yyyymmddhh))
        print("[ERR] Check length!")
        usage()
        sys.exit(1)

    print("[INFO] Processing date/time %s" %(yyyymmddhh))
    try:
        year = int(yyyymmddhh[0:4])
        month = int(yyyymmddhh[4:6])
        day = int(yyyymmddhh[6:8])
        hour = int(yyyymmddhh[8:10])
        validdt = datetime.datetime(year, month, day, hour)
    except ValueError:
        print("[ERR] Cannot process date/time %s" %(yyyymmddhh))
        usage()
        sys.exit(1)

    # NOTE:  Pylint complains about netCDF4 not having Dataset as a
    # member.  But this is false -- Dataset is a class, and we use
    # the constructor below.  Since this appears to be a bug, we
    # disable the no-member check in this function.
    # pylint: disable=no-member

    # See if lvt_nc file can be opened
    lvt_nc = sys.argv[2]
    ncid_lvt = netCDF4.Dataset(lvt_nc, mode='r',
                               format='NETCDF4_CLASSIC')
    ncid_lvt.close()

    # See if ldt_nc file can be opened
    ldt_nc = sys.argv[3]
    ncid_ldt = netCDF4.Dataset(ldt_nc, mode='r',
                               format='NETCDF4_CLASSIC')
    ncid_ldt.close()
    # pylint: enable=no-member

    return validdt, lvt_nc, ldt_nc

#------------------------------------------------------------------------------
# Main Driver.
if __name__ == "__main__":

    # Process command line arguments
    VALIDDT, LVT_NC, LDT_NC = read_cmd_args()

    # Create SURF object
    SURF = Nc2Surf(LVT_NC, LDT_NC)

    # Generate glu_sst SURF file
    FILE_TYPE = "_glu_sst"
    # FUTURE...Add support for SST anomaly
    VARFIELDS = ["water_temp:LVT"]
    SURFFILE = "%4.4d%2.2d%2.2dT%2.2d00Z%s" \
        % (VALIDDT.year,
           VALIDDT.month,
           VALIDDT.day,
           VALIDDT.hour,
           FILE_TYPE)
    SURF.create_surf_file(FILE_TYPE, VARFIELDS, SURFFILE)

    # Generate glu_ice SURF file
    FILE_TYPE = "_glu_ice"
    VARFIELDS = ["aice:LVT", "hi:LVT"]
    SURFFILE = "%4.4d%2.2d%2.2dT%2.2d00Z%s" \
        % (VALIDDT.year,
           VALIDDT.month,
           VALIDDT.day,
           VALIDDT.hour,
           FILE_TYPE)
    SURF.create_surf_file(FILE_TYPE, VARFIELDS, SURFFILE)

    # Generate glu_snow SURF file
    FILE_TYPE = "_glu_snow"
    VARFIELDS = ["SWE_inst:LVT", "SnowDepth_inst:LVT"]
    SURFFILE = "%4.4d%2.2d%2.2dT%2.2d00Z%s" \
        % (VALIDDT.year,
           VALIDDT.month,
           VALIDDT.day,
           VALIDDT.hour,
           FILE_TYPE)
    SURF.create_surf_file(FILE_TYPE, VARFIELDS, SURFFILE)

    # Generate glu_smc SURF file
    FILE_TYPE = "_glu_smc"
    VARFIELDS = ["SoilMoist_inst:LVT", "SoilTemp_inst:LVT",
                 "AvgSurfT_inst:LVT",
                 "JULES_SM_WILT:LDT", "JULES_SM_CRIT:LDT"]
    SURFFILE = "%4.4d%2.2d%2.2dT%2.2d00Z%s" \
        % (VALIDDT.year,
           VALIDDT.month,
           VALIDDT.day,
           VALIDDT.hour,
           FILE_TYPE)
    SURF.create_surf_file(FILE_TYPE, VARFIELDS, SURFFILE)

    print("[INFO] convert_nc2surf.py completed!")
    sys.exit(0)
