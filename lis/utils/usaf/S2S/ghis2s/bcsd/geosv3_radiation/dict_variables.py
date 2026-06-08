#!/usr/bin/env python3
#pylint: disable=invalid-name
#pylint: disable=consider-iterating-dictionary
""" Centralized functions to get model variable names used in the SERVIR-S2S BCSD project

Shared dictionary containing all of the attributes of a given 'key' variable name. Also includes
helper functions to quickly extract an attribute. Key variable names are set by me, and should not 
appear in any results, but used as keys to return the information needed.

Example (within Python script):
    from dict_variables import get_all_key_variable_names, get_hydroscs_name

    # Return all key variable names
    key_variable_names = get_all_key_variable_names()
    # Returns ['T', 'Q', 'PS', 'V', 'U', 'LWdown', 'SWdown', 'Pr']

    # Return HydroSFS name, given the key variable name
    var_name_hydrosfs = get_hydrosfs_name('T')
    # Returns 'Tair'

Author: Ryan Zamora
"""

# Labels for the dictionary; if they ever need to be changed, do so here.
# However, in general, it's best to use the helper functions to get attributes
cfsv2_label = 'CFSv2 Name'
cfsv2_var_label = 'CFSv2 Variable Name'
geosv3_label = 'GEOSv3 Name'
hydroscs_label = 'HydroSCS Name'
hydrosfs_label = 'HydroSFS Name'
longname_label = 'Long Name'
units_label = 'Units'
interp_method_label = 'Interpolation Method'
disaggregation_type_label = 'Disaggregation Type'

# Dictionary data
dictionary = {
    'T' : {
        cfsv2_label     : 'tmp2m',
        cfsv2_var_label : 'tmp2m',
        geosv3_label    : 'T2M',
        hydroscs_label  : 'Tair',
        hydrosfs_label  : 'Tair',
        longname_label  : 'Near Surface Air Temperature',
        units_label     : 'K',
        interp_method_label       : 'bilinear',
        disaggregation_type_label : 'TEMP',
    },
    'Q' : {
        cfsv2_label     : 'q2m',
        cfsv2_var_label : 'q2m',
        geosv3_label    : 'Q2M',
        hydroscs_label  : 'Qair',
        hydrosfs_label  : 'Qair',
        longname_label  : 'Near Surface Specific Humidity',
        units_label     : 'kg/kg',
        interp_method_label       : 'conservative',
        disaggregation_type_label : 'TEMP',
    },
    'PS' : {
        cfsv2_label     : 'pressfc',
        cfsv2_var_label : 'pressfc',
        geosv3_label    : 'PS',
        hydroscs_label  : 'Psurf',
        hydrosfs_label  : 'Psurf',
        longname_label  : 'Surface Pressure',
        units_label     : 'Pa',
        interp_method_label       : 'conservative',
        disaggregation_type_label : 'TEMP',
    },
    'V' : {
        cfsv2_label     : 'wnd10m',
        cfsv2_var_label : 'v10',
        geosv3_label    : 'V10M',
        hydroscs_label  : 'Wind_N',
        hydrosfs_label  : 'Wind_N',
        longname_label  : 'Near Surface Northward Wind',
        units_label     : 'm/s',
        interp_method_label       : 'bilinear',
        disaggregation_type_label : 'TEMP',
    },
    'U' : {
        cfsv2_label     : 'wnd10m',
        cfsv2_var_label : 'u10',
        geosv3_label    : 'U10M',
        hydroscs_label  : 'Wind_E',
        hydrosfs_label  : 'Wind_E',
        longname_label  : 'Near Surface Eastward Wind',
        units_label     : 'm/s',
        interp_method_label       : 'bilinear',
        disaggregation_type_label : 'TEMP',
    },
    'LWdown' : {
        cfsv2_label     : 'dlwsfc',
        cfsv2_var_label : 'dlwsfc',
        geosv3_label    : 'LWS',
        hydroscs_label  : 'LWdown',
        hydrosfs_label  : 'LWdown',
        longname_label  : 'Downward Longwave Radiation',
        units_label     : 'W/m^2',
        interp_method_label       : 'conservative',
        disaggregation_type_label : 'TEMP',
    },
    'SWdown' : {
        cfsv2_label     : 'dswsfc',
        cfsv2_var_label : 'dswsfc',
        geosv3_label    : 'SWGDWN',
        hydroscs_label  : 'SWdown',
        hydrosfs_label  : 'SWdown',
        longname_label  : 'Downward Shortwave Radiation',
        units_label     : 'W/m^2',
        interp_method_label       : 'conservative',
        disaggregation_type_label : 'PRCP',
    },
    'Pr' : {
        cfsv2_label     : 'prate',
        cfsv2_var_label : 'prate',
        geosv3_label    : 'PRECTOTCORR',
        hydroscs_label  : 'Rainf',
        hydrosfs_label  : 'Rainf',
        longname_label  : 'Precipitation',
        units_label     : 'kg/m^2/s',
        interp_method_label       : 'conservative',
        disaggregation_type_label : 'PRCP',
    },
}

#
# Functions to return a name, given a key variable name
#
def get_cfsv2_name(key_variable_name: str):
    """
    Helper function to return the CFSv2 variable file name, given the key variable name 
    This is the name used in grib filenames. Differs from get_cfsv2_var_name() only for U & V
    """
    return dictionary[key_variable_name][cfsv2_label]

def get_cfsv2_var_name(key_variable_name: str):
    """
    Helper function to return the CFSv2 variable name, given the key variable name 
    This is the variable name used. Differs from get_cfsv2_name() only for U & V
    """
    return dictionary[key_variable_name][cfsv2_var_label]

def get_geosv3_name(key_variable_name: str):
    """ Helper function to return the GEOSv3 variable name, given the key variable name """
    return dictionary[key_variable_name][geosv3_label]

def get_hydroscs_name(key_variable_name: str):
    """ Helper function to return the HydroSCS variable name, given the key variable name """
    return dictionary[key_variable_name][hydroscs_label]

def get_hydrosfs_name(key_variable_name: str):
    """ Helper function to return the HydroSFS variable name, given the key variable name """
    return dictionary[key_variable_name][hydrosfs_label]

def get_long_name(key_variable_name: str):
    """ Helper function to return the long variable name, given the key variable name """
    return dictionary[key_variable_name][longname_label]

def get_units(key_variable_name: str):
    """ Helper function to return the units, given the key variable name """
    return dictionary[key_variable_name][units_label]

def get_interp_method(key_variable_name: str):
    """ Helper function to return the interpolation method, given the key variable name """
    return dictionary[key_variable_name][interp_method_label]

def get_disaggregation_type(key_variable_name: str):
    """ Helper function to return the disaggregation type, given the key variable name """
    return dictionary[key_variable_name][disaggregation_type_label]

#
# Functions to return all names
#
def get_all_key_variable_names():
    """ Helper function to return all country abbreviations """
    return list(dictionary.keys())

def get_all_cfsv2_names():
    """ Helper function to return all of the CFSv2 file variable names """
    return [get_cfsv2_name(key_variable) for key_variable in dictionary.keys()]

def get_all_cfsv2_var_names():
    """ Helper function to return all of the CFSv2 variable names """
    return [get_cfsv2_var_name(key_variable) for key_variable in dictionary.keys()]

def get_all_geosv3_names():
    """ Helper function to return all of the GEOSv3 variable names """
    return [get_geosv3_name(key_variable) for key_variable in dictionary.keys()]

def get_all_hydroscs_names():
    """ Helper function to return all of the HydroSCS variable names """
    return [get_hydroscs_name(key_variable) for key_variable in dictionary.keys()]

def get_all_hydrosfs_names():
    """ Helper function to return all of the HydroSFS variable names """
    return [get_hydrosfs_name(key_variable) for key_variable in dictionary.keys()]

def get_all_long_names():
    """ Helper function to return all of the long variable names """
    return [get_long_name(key_variable) for key_variable in dictionary.keys()]

def get_all_units():
    """ Helper function to return all of the variable units """
    return [get_units(key_variable) for key_variable in dictionary.keys()]

def get_all_interp_method():
    """ Helper function to return all of the interpolation methods """
    return [get_interp_method(key_variable) for key_variable in dictionary.keys()]

def get_all_disaggregation_type():
    """ Helper function to return all of the disaggregation types """
    return [get_disaggregation_type(key_variable) for key_variable in dictionary.keys()]
#
# Functions to return dictionary for renaming variables
#
def get_rename_dictionary_hydroscs_to_hydrosfs():
    """ Helper function to return renaming dictionary from HydroSCS to HydroSFS """
    return {get_hydroscs_name(key_variable):get_hydrosfs_name(key_variable)
            for key_variable in dictionary.keys()}

def get_rename_dictionary_cfsv2_to_hydrosfs():
    """ Helper function to return renaming dictionary from GEOSv3 to HydroSFS """
    return {get_cfsv2_name(key_variable):get_hydrosfs_name(key_variable)
            for key_variable in dictionary.keys()}

def get_rename_dictionary_geosv3_to_hydrosfs():
    """ Helper function to return renaming dictionary from GEOSv3 to HydroSFS """
    return {get_geosv3_name(key_variable):get_hydrosfs_name(key_variable)
            for key_variable in dictionary.keys()}
