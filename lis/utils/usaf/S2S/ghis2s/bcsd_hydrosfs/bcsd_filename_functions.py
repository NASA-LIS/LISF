#!/usr/bin/env python3
"""
Functions to get file information from the S2S E2ES bcsd_fcst substep

Script
------
bcsd_filename_functions.py

Requirements
------------
Python 3.9

Revision History
----------------
07 Mar 2024: Ryan Zamora, first version

"""
# ----------- #
# - Modules - #
# ----------- #
import calendar
from datetime import datetime
from dateutil.relativedelta import relativedelta
import numpy as np

# ------------- #
# - Variables - #
# ------------- #
str_resolution = '5km'
final_dataset_name = 'HydroSFS'

# ------------------------- #
# - Valid Value Functions - #
# ------------------------- #
def _check_valid_fcst_init_year(value : int):
    """ Checks for valid fcst_init_year value"""
    if not isinstance(value, int):
        raise TypeError('fcst_init_year must be an integer')
    if value not in range(1000,10000):
        raise ValueError('fcst_init_year must be a 4-digit integer')

def _check_valid_fcst_init_month(value : int):
    """Checks for valid fcst_init_month value"""
    if not isinstance(value, int):
        raise TypeError('fcst_init_month must be an integer')
    if value not in range(1,13):
        raise ValueError('fcst_init_month must be an integer with value 1-12')

def _check_valid_ensemble_number(value : int):
    """Checks for valid ensemble_number value"""
    if not isinstance(value, int):
        raise TypeError('ensemble_number must be an integer')
    if value <= 0:
        raise ValueError('ensemble_number must be a positive non-zero integer')

def _check_valid_lead_number(value : int):
    """Checks for valid lead_number value"""
    if not isinstance(value, int):
        raise TypeError('lead_number must be an integer')
    if value < 0:
        raise ValueError('lead_number must be a positive integer')

def _check_valid_dataset_type(value : str):
    """Checks for valid dataset_type value"""
    valid_dataset_types = ['hindcast', 'forecast']
    if not isinstance(value, str):
        raise TypeError('dataset_type must be a string')
    if value not in valid_dataset_types:
        raise ValueError(f"dataset_type must be a string with value: {valid_dataset_types}")


# ------------- #
# - Functions - #
# ------------- #
def generate_file_information(
    file_directory: str, file_name: str, output_type: str = 'file_path') -> str:
    """
    Reusable helper function that returns file information of the given output type
    
    Parameters
    ----------
    file_directory (int) : Directory where the file is located
    file_name      (int) : Name of the file with no path information
    output_type    (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']

    Returns
    -------
    file_information (str) : File information of the given output type

    Examples
    --------

    """
    
    # Valid values for output_type
    valid_output_types = ['file_directory', 'file_name', 'file_path']
    
    # Checks for valid output_type type
    if not isinstance(output_type, str):
        raise TypeError('output_type must be a string')

    # Gather file information of the given output_type
    if output_type == 'file_directory':
        # Directory where the file is located
        file_information = file_directory
    elif output_type == 'file_name':
        # Name of the file with no path information
        file_information = file_name
    elif output_type == 'file_path':
        # Path to the file (including file name)
        file_sep = '/' if file_directory[-1] != '/' else ''
        file_information = f"{file_directory}{file_sep}{file_name}"
    else:
        raise ValueError(f"output_type must be a string with value: {valid_output_types}")
    
    return file_information

def get_forecast_date_strings_from_lead_number(
    fcst_init_year: int, fcst_init_month: int, lead_number: int) -> tuple[str, str]:
    """
    Returns the associated time strings of forecast

    Parameters
    ----------
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    lead_number     (int) : Lead month number

    Returns
    -------
    fcst_year  (str) : Year of the downloaded forecast
    fcst_month (str) : Month of the downloaded forecast

    Examples
    --------
    >>> fcst_year,fcst_month = get_forecast_date_strings_from_lead_number(1991, 5, 2)
    >>> fcst_year
    '1991'
    >>> fcst_month
    '07'
    """

    # Check for valid parameters
    _check_valid_fcst_init_year(fcst_init_year)
    _check_valid_fcst_init_month(fcst_init_month)
    _check_valid_lead_number(lead_number)

    # Generate datetime of forecast
    new_datetime = datetime(fcst_init_year, fcst_init_month,1) + relativedelta(months=lead_number)

    # Formatted forecast year
    fcst_year  = new_datetime.strftime('%Y')

    # Formatted forecast month
    fcst_month = new_datetime.strftime('%m')

    return fcst_year, fcst_month

def get_forecast_date_strings_from_lead_day(
    fcst_init_year: int, fcst_init_month: int, lead_day: int) -> tuple[str, str, str]:
    """
    Returns the associated time strings of forecast

    Parameters
    ----------
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    lead_day        (int) : Lead day number

    Returns
    -------
    fcst_year  (str) : Year of the downloaded forecast
    fcst_month (str) : Month of the downloaded forecast
    fcst_day   (str) : Day of the downloaded forecast

    Examples
    --------
    >>> fcst_year,fcst_month, fcst_day = get_forecast_date_strings_from_lead_number(1991, 5, 2)
    >>> fcst_year
    '1991'
    >>> fcst_month
    '07'
    >>> fcst_day
    '02'
    """
    
    # Generate datetime of forecast
    new_datetime = datetime(fcst_init_year, fcst_init_month,1) + relativedelta(days=lead_day)

    # Formatted forecast year
    fcst_year  = new_datetime.strftime('%Y')

    # Formatted forecast month
    fcst_month = new_datetime.strftime('%m')
    
    # Formatted forecast month
    fcst_day = new_datetime.strftime('%d')

    return fcst_year, fcst_month, fcst_day

# --------------------------------- #
# - DOWNLOADED HYDROSCS FUNCTIONS - #
# --------------------------------- #
def generate_downloaded_hydroscs_filename(
    dir_hydroscs : str, year : int, month : int, day : int, hour : int,
    output_type: str = 'file_path') -> str:
    """
    Returns the associated downloaded HydroSCS file information

    Parameters
    ----------
    year         (int) : Observational year
    month        (int) : Observational month
    day          (int) : Observational day
    hour         (int) : Observational hour    
    dir_hydroscs (str) : Directory of HydroSCS downloads
    output_type  (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']    

    Returns
    -------
    file_information (str) : File information of downloaded dataset for given arguments

    Examples
    --------
    >>> 
    """

    # File directory and name information
    file_pfx = 'HydroSCS_'
    file_sfx = '00.d01.nc'
    
    file_directory = f"{dir_hydroscs}/{year:04d}{month:02d}"
    file_name = f"{file_pfx}{year:04d}{month:02d}{day:02d}{hour:02d}{file_sfx}"
    
    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

# ------------------------------- #
# - RAW MONTHLY FORECAST FUNCTIONS - #
# ------------------------------- #
def generate_raw_monthly_forecast_filename(
    dir_e2es : str, fcst_name : str, dataset_type  : str, fcst_init_year : int, fcst_init_month : int,
    ensemble_number : int, lead_number : int, output_type: str = 'file_path') -> str:
    """
    Returns the associated raw monthly forecast file information

    Parameters
    ----------
    dir_e2es        (str) : Main directory of the E2ES
    dataset_type    (str) : Type of the data (hindcast or forecast)
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    ensemble_number (int) : Ensemble number
    lead_number     (int) : Lead month number
    output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']    

    Returns
    -------
    file_information (str) : File information of raw monthly forecast dataset for given arguments

    Examples
    --------
    >>> x = generate_raw_monthly_forecast_filename('E2ES/', 'hindcast', 1991, 5, 1, 2)
    >>> x
    'E2ES//hindcast/bcsd_fcst/CFSv2_5km/raw/Monthly/may01/1991/ens1/may01.cfsv2.199107.nc'
    >>> y = generate_raw_monthly_forecast_filename('E2ES/', 'forecast', 2023, 9, 5, 0)
    >>> y
    'E2ES//bcsd_fcst/CFSv2_5km/raw/Monthly/sep01/2023/ens5/sep01.cfsv2.202309.nc'
    """

    # Check for valid parameters
    _check_valid_dataset_type(dataset_type)
    _check_valid_fcst_init_year(fcst_init_year)
    _check_valid_fcst_init_month(fcst_init_month)
    _check_valid_ensemble_number(ensemble_number)
    _check_valid_lead_number(lead_number)

    # Directory with dataset_type info
    dir_dataset = f"{dir_e2es}/{dataset_type}"

    # Month name info
    mon01 = calendar.month_abbr[fcst_init_month].lower() + '01'

    # Generate fcst_year and fcst_month
    fcst_year, fcst_month = get_forecast_date_strings_from_lead_number(
        fcst_init_year, fcst_init_month, lead_number)

    # File directory and name information
    file_directory = (
        f"{dir_dataset}/{fcst_name}_{str_resolution}/Monthly/"
        f"{mon01}/{fcst_init_year}/ens{ensemble_number}")
    file_name = f"{mon01}.{fcst_name.lower()}.{fcst_year}{fcst_month}.nc"

    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

# ---------------------------------- #
# - RAW MONTHLY HYDROSCS FUNCTIONS - #
# ---------------------------------- #
def generate_raw_monthly_hydroscs_filename(
    dir_e2es : str, dataset_type  : str, year : int, month : int,
    output_type: str = 'file_path') -> str:
    """
    Returns the associated downloaded HydroSCS file information

    Parameters
    ----------
    dir_e2es     (str) : Directory of HydroSCS downloads
    dataset_type (str) : Type of the data (hindcast or forecast)
    year         (int) : Observational year
    month        (int) : Observational month
    output_type  (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']    

    Returns
    -------
    file_information (str) : File information of downloaded HydroSCS dataset for given arguments

    Examples
    --------
    >>> 
    """

    # Check for valid parameters
    # _check_valid_dataset_type(dataset_type)
    # _check_valid_year(year)
    # _check_valid_month(month)
    # _check_valid_output_type(output_type)

    # Directory with dataset_type info
    dir_dataset = f"{dir_e2es}/{dataset_type}"
    
    # File directory and name information
    file_pfx = 'HydroSCS_'
    file_sfx = '.nc'
    
    # File directory and name information
    file_directory = f"{dir_dataset}/HydroSCS/Monthly"
    file_name = f"{file_pfx}{year:04d}{month:02d}{file_sfx}"
    
    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

def generate_raw_6hourly_forecast_filename(
    dir_e2es : str, fcst_name : str, dataset_type : str, fcst_init_year : int, fcst_init_month : int,
    ensemble_number : int, lead_day : int, output_type: str = 'file_path') -> str:
    """
    Returns the associated raw 6-hourly forecast file information

    Parameters
    ----------
    dir_e2es        (str) : Main directory of the E2ES
    dataset_type    (str) : Type of the data (hindcast or forecast)
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    ensemble_number (int) : Ensemble number
    lead_day        (int) : Lead day number
    output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']    

    Returns
    -------
    file_information (str) : File information of raw 6-hourly forecast dataset for given arguments

    Examples
    --------
    >>> x = generate_raw_6hourly_forecast_filename('~/E2ES', 'hindcast', 1991, 5, 1, 0)
    >>> x
    '~/E2ES/hindcast/bcsd_fcst/CFSv2_5km/raw/6-Hourly/may01/1991/ens1/may01.cfsv2.199105.nc'
    >>> y = generate_raw_6hourly_forecast_filename('~/E2ES', 'forecast', 2023, 9, 5, 0)
    >>> y
    '~/E2ES/bcsd_fcst/CFSv2_5km/raw/6-Hourly/sep01/2023/ens5/sep01.cfsv2.202309.nc'
    """

    # Directory with dataset_type info
    dir_dataset = f"{dir_e2es}/{dataset_type}"

    # Month name info
    mon01 = calendar.month_abbr[fcst_init_month].lower() + '01'

    # Generate fcst_year and fcst_month
    fcst_year, fcst_month, fcst_day = get_forecast_date_strings_from_lead_day(
        fcst_init_year, fcst_init_month, lead_day)

    # File directory and name information
    file_directory = (
        f"{dir_dataset}/{fcst_name}_{str_resolution}/6-Hourly/"
        f"{mon01}/{fcst_init_year}/ens{ensemble_number}")
    file_name = f"{mon01}.{fcst_name.lower()}.{fcst_year}{fcst_month}{fcst_day}.nc"

    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

# ----------------------------------- #
# - RAW CLIMATOLOGY FORECAST FUNCTIONS - #
# ----------------------------------- #
def generate_raw_climatology_forecast_filename(
    dir_e2es : str, fcst_name : str, dataset_type : str, fcst_init_month : int, variable_name : str,
    output_type: str = 'file_path') -> str:
    """
    Returns the associated raw climatology forecast file information

    Parameters
    ----------
    dir_e2es        (str) : Main directory of the E2ES
    dataset_type    (str) : Type of the data (hindcast or forecast)
    fcst_init_month (int) : Initial forecast month
    variable_name   (str) : Variable name (observation)
    output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']

    Returns
    -------
    file_information (str) : File information of raw climatology forecast dataset for given arguments

    Examples
    --------
    >>> x = generate_raw_climatology_forecast_filename('~/E2ES', 'hindcast', 5, 'T2M')
    >>> x
    '~/E2ES/hindcast/bcsd_fcst/CFSv2_5km/raw/Climatology/may01/T2M_fcst_clim.nc'
    >>> y = generate_raw_climatology_forecast_filename('~/E2ES', 'forecast', 9, 'PS')
    >>> y
    '~/E2ES/bcsd_fcst/CFSv2_5km/raw/Climatology/sep01/PS_fcst_clim.nc'
    """

    # Check for valid parameters
    _check_valid_dataset_type(dataset_type)
    _check_valid_fcst_init_month(fcst_init_month)

    # Directory with dataset_type info
    dir_dataset = f"{dir_e2es}/{dataset_type}"

    # Month name info
    mon01 = calendar.month_abbr[fcst_init_month].lower() + '01'

    # File directory and name information
    file_directory = f"{dir_dataset}/{fcst_name}_{str_resolution}/Climatology/{mon01}"
    file_name = f"{variable_name}_fcst_clim.nc"

    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

# -------------------------------------- #
# - RAW CLIMATOLOGY HYDROSCS FUNCTIONS - #
# -------------------------------------- #
def generate_raw_climatology_hydroscs_filename(
    dir_e2es : str, dataset_type : str, var_name: str, output_type: str = 'file_path') -> str:
    """
    Returns the associated downloaded HydroSCS file information

    Parameters
    ----------
    dir_e2es     (str) : Directory of HydroSCS downloads
    dataset_type (str) : Type of the data (hindcast or forecast)
    var_name     (str) : Variable name
    output_type  (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path'] 

    Returns
    -------
    file_information (str) : File information of downloaded HydroSCS dataset for given arguments

    Examples
    --------
    >>> 
    """

    # Directory with dataset_type info
    dir_dataset = f"{dir_e2es}/{dataset_type}"
    
    # File directory and name information
    file_directory = f"{dir_dataset}/HydroSCS/Climatology"
    file_name = f"{var_name}_obs_clim.nc"
    
    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

# -------------------------------- #
# - BCSD MONTHLY FORECAST FUNCTIONS - #
# -------------------------------- #
def generate_bcsd_monthly_forecast_filename(
    dir_e2es : str, fcst_name : str, dataset_type : str, fcst_init_year : int, fcst_init_month : int,
    variable_name : str, output_type: str = 'file_path') -> str:
    """
    Returns the associated bcsd monthly forecast file information

    Parameters
    ----------
    dir_e2es        (str) : Main directory of the E2ES
    dataset_type    (str) : Type of the data (hindcast or forecast)
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    variable_name   (str) : Variable name (observation)
    output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']

    Returns
    -------
    file_information (str) : File information of bcsd monthly forecast dataset for given arguments

    Examples
    --------
    >>> x = generate_bcsd_monthly_forecast_filename('~/E2ES', 'hindcast', 1991, 5, 'T2M')
    >>> x
    '~/E2ES/hindcast/bcsd_fcst/CFSv2_5km/bcsd/Monthly/may01/T2M.CFSv2.may01_1991_1991.nc'
    >>> y = generate_bcsd_monthly_forecast_filename('~/E2ES', 'forecast', 2023, 9, 'PS')
    >>> y
    '~/E2ES/bcsd_fcst/CFSv2_5km/bcsd/Monthly/sep01/PS.CFSv2.sep01_2023_2023.nc'
    """

    # Check for valid parameters
    _check_valid_dataset_type(dataset_type)
    _check_valid_fcst_init_year(fcst_init_year)
    _check_valid_fcst_init_month(fcst_init_month)

    # Directory with dataset_type info
    dir_dataset = f"{dir_e2es}/{dataset_type}"

    # Month name info
    mon01 = calendar.month_abbr[fcst_init_month].lower() + '01'

    # File directory and name information
    file_directory = f"{dir_dataset}/{final_dataset_name}_{fcst_name}/Monthly/{mon01}"
    file_name = f"{variable_name}.{fcst_name}.{mon01}_{fcst_init_year}_{fcst_init_year}.nc"

    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information

# --------------------------------- #
# - BCSD 6-HOURLY FORECAST FUNCTIONS - #
# --------------------------------- #
def generate_bcsd_6hourly_forecast_filename(
    dir_e2es : str, fcst_name : str, dataset_type : str, fcst_init_year : int, fcst_init_month : int,
    ensemble_number : int, lead_day, variable_name : str,
    output_type: str = 'file_path') -> str:
    """
    Returns the associated bcsd 6-hourly forecast file information

    Parameters
    ----------
    dir_e2es        (str) : Main directory of the E2ES
    dataset_type    (str) : Type of the data (hindcast or forecast)
    fcst_init_year  (int) : Initial forecast year
    fcst_init_month (int) : Initial forecast month
    ensemble_number (int) : Ensemble number
    lead_number     (int) : Lead month number
    variable_name   (str) : Variable name (observation)
    output_type     (str) : Type of output wanted ['file_directory', 'file_name', or 'file_path']

    Returns
    -------
    file_information (str) : File information of bcsd 6-hourly forecast dataset for given arguments

    Examples
    --------
    >>> x = generate_bcsd_6hourly_forecast_filename('~/E2ES', 'hindcast', 1991, 5, 1, 2, 'T2M')
    >>> x
    '~/E2ES/hindcast/bcsd_fcst/CFSv2_5km/bcsd/6-Hourly/may01/1991/ens1/T2M.199107.nc4'
    >>> y = generate_bcsd_6hourly_forecast_filename('~/E2ES', 'forecast', 2023, 9, 5, 0, 'PS')
    >>> y
    '~/E2ES/bcsd_fcst/CFSv2_5km/bcsd/6-Hourly/sep01/2023/ens5/PS.202309.nc4'
    """

    # Check for valid parameters
    _check_valid_dataset_type(dataset_type)
    _check_valid_fcst_init_year(fcst_init_year)
    _check_valid_fcst_init_month(fcst_init_month)
    _check_valid_ensemble_number(ensemble_number)

    # Directory with dataset_type info
    dir_dataset = f"{dir_e2es}/{dataset_type}"

    # Month name info
    mon01 = calendar.month_abbr[fcst_init_month].lower() + '01'

    # Generate fcst_year and fcst_month
    fcst_year, fcst_month, fcst_day = get_forecast_date_strings_from_lead_day(
        fcst_init_year, fcst_init_month, lead_day)

    # File directory and name information
    file_directory = (
        f"{dir_dataset}/{final_dataset_name}_{fcst_name}/6-Hourly/"
        f"{mon01}/{fcst_init_year}/ens{ensemble_number}")
    file_name = f"{final_dataset_name}_{fcst_name}_{variable_name}_{fcst_year}{fcst_month}{fcst_day}.nc4"
    
    # Generate file information based on output_type
    file_information = generate_file_information(file_directory, file_name, output_type)
    
    return file_information
