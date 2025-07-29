#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: metricslib.py
#
# PURPOSE: Contains library metrhods used by multiple s2smetrics scripts.
# Based on All_functions.py by Shrad Shukla/UCSB.
#
# REVISION HISTORY:
# 25 Oct 2021: Eric Kemp/SSAI, first version
# 02 Jun 2023: K. Arsenault + S. Mahanama, updated 557 WW file conventions.
#
#------------------------------------------------------------------------------
"""
# Standard modules
import sys
import os
import glob
import xarray as xr
import numpy as np

def _get_snow_cover(sel_cim_data):
    """Return snow cover [%]."""
    return sel_cim_data.Snowcover_tavg

def _get_snow_depth(sel_cim_data):
    """Return snow depth [m]."""
    return sel_cim_data.SnowDepth_tavg

def _get_swe(sel_cim_data):
    """Return Snow Water Equivalent [kg/m^2]."""
    return sel_cim_data.SWE_tavg

def _get_surface_st(sel_cim_data):
    """Return top layer soil temperature."""
    return sel_cim_data.SoilTemp_tavg.isel(soil_layer=0)

def _get_surface_sm(sel_cim_data):
    """Return surface soil moisture."""
    return sel_cim_data.SoilMoist_tavg.isel(soil_layer=0)

def _get_tws(sel_cim_data):
    """Return terrestrial water storage."""
    return sel_cim_data.TWS_tavg

def _get_precip(sel_cim_data):
    """Return total precipitation."""
    return sel_cim_data.TotalPrecip_acc

def _get_air_t(sel_cim_data):
    """Return air temperature."""
    return sel_cim_data.Tair_f_tavg

def _get_et(sel_cim_data):
    """Return evapotranspiration."""
    return sel_cim_data.Evap_tavg

def _get_streamflow(sel_cim_data):
    """Return streamflow."""
    return sel_cim_data.Streamflow_tavg

def sel_var(sel_cim_data, var_name, model):
    """Selects climatology for the given variable."""
    if var_name == "RZSM":
        if model == "CLSM":
            # for clsm the layer-2 is rootzone soil moisture
            var_sel_clim_data = sel_cim_data.SoilMoist_tavg.isel(soil_layer=1)
        elif model in ('NOAHMP', 'NoahMP'):
            term1 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=0) * 0.1
            term2 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=1) * 0.3
            term3 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=2) * 0.6
            var_sel_clim_data = term1 + term2 + term3
        else:
            print(f"[ERR] Unknown model {model}")
            sys.exit(1)

    elif var_name == "TOP40ST":
        if model in ('NOAHMP', 'NoahMP'):
            term1 = sel_cim_data.SoilTemp_tavg.isel(soil_layer=0) * 0.1
            term2 = sel_cim_data.SoilTemp_tavg.isel(soil_layer=1) * 0.3
            var_sel_clim_data = term1 + term2
        else:
            print(f"[ERR] Unknown model {model}")
            sys.exit(1)

    elif var_name == "TOP40SM":
        if model in ('NOAHMP', 'NoahMP'):
            term1 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=0) * 0.1
            term2 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=1) * 0.3
            var_sel_clim_data = term1 + term2
        else:
            print(f"[ERR] Unknown model {model}")
            sys.exit(1)

    elif var_name == 'Total-SM':
        if model == 'CLSM':
            # for clsm the total soil moisture is in the third layer
            var_sel_clim_data = sel_cim_data.SoilMoist_tavg.isel(soil_layer=2)
        elif model in ("NOAHMP", "NoahMP"):
            term1 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=0) * 0.05
            term2 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=1) * 0.15
            term3 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=2) * 0.3
            term4 = sel_cim_data.SoilMoist_tavg.isel(soil_layer=3) * 0.5
            var_sel_clim_data = term1 + term2 + term3 + term4
        else:
            print(f"[ERR] Unknown model {model}")
            sys.exit(1)

    elif var_name == 'Total-Runoff':
        ## Adding total surface runoff with sub-surface runoff
        var_sel_clim_data = sel_cim_data.Qs_tavg + sel_cim_data.Qsb_tavg

    else:
        # Pylint complains about too many elifs in original code. So, for
        # the trivial copies, we combine the remaining elifs into a
        # python dictionary of functions, and call the appropriate function
        # based on the variable name.
        selections = {
            "SFCSM" : _get_surface_sm,
            "TWS" : _get_tws,
            "Precip" : _get_precip,
            "AirT" : _get_air_t,
            "ET" : _get_et,
            "Streamflow" : _get_streamflow,
            "SFCST": _get_surface_st,
            "SWE": _get_swe,
            "SnowDepth": _get_snow_depth,
            "SnowCover": _get_snow_cover
        }
        try:
            var_sel_clim_data = selections[var_name](sel_cim_data)
        except KeyError:
            print(f"[ERR] Unknown var_name: {var_name}")
            sys.exit(1)

    return var_sel_clim_data

def compute_anomaly (target_fcst_data, fcst_clim):
    ''' compute anomaly'''
    out_var = target_fcst_data - fcst_clim
    return out_var

def compute_sanomaly (target_fcst_data, fcst_clim, fcst_std):
    ''' computes standerdized anomaly'''
    out_var = np.ones(target_fcst_data.shape)*-9999.
    fcst_clim_bc = np.broadcast_to(fcst_clim, target_fcst_data.shape)
    fcst_std_bc = np.broadcast_to(fcst_std, target_fcst_data.shape)
    valid_std = np.logical_and(fcst_std_bc > 0, ~np.isnan(fcst_std_bc))
    out_var[valid_std] = (target_fcst_data[valid_std] - fcst_clim_bc[valid_std]) / fcst_std_bc[valid_std]

    return out_var

# Units for variable anomalies.  Standardized anomalies will be dimensionless.
UNITS_ANOM = {
    "RZSM" : "m3 m-3",
    "TOP40SM" : "m3 m-3",
    "SFCSM" : "m3 m-3",
    "TWS" : "mm",
    "Precip" : "kg m-2",
    "AirT" : "K",
    "ET" : "kg m-2 s-1",
    "Streamflow" : "m3 s-1",
    "SWE" : "kg m-2",
    "SnowDepth" : "kg m-2",
    "TOP40ST" : "K",
}
UNITS_SANOM = {
    "RZSM" : "1",
    "TOP40SM" : "1",
    "SFCSM" : "1",
    "TWS" : "1",
    "Precip" : "1",
    "AirT" : "1",
    "ET" : "1",
    "Streamflow" : "1",
    "SWE": "1",
    "SnowDepth": "1",
    "TOP40ST" : "1",    
}

LONG_NAMES_SANOM = {
    "RZSM" : "Root zone soil moisture standardized anomaly",
    "TOP40SM" : "Top 0-40 cm soil moisture standardized anomaly",
    "SFCSM" : "Surface soil moisture standardized anomaly",
    "TWS" : "Terrestrial water storage standardized anomaly",
    "Precip" : "Total precipitation amount standardized anomaly",
    "AirT" : "Air temperature standardized anomaly",
    "ET" : "Total evapotranspiration standardized anomaly",
    "Streamflow" : "Streamflow standardized anomaly",
    "SWE": "Snow water equivalent standardized anomaly",
    "SnowDepth": " Snowdepth standardized anomaly",
    "TOP40ST" : "Top 0-40 cm soil temperature standardized anomaly",
}

LONG_NAMES_ANOM = {
    "RZSM" : "Root zone soil moisture anomaly",
    "TOP40SM" : "Top 0-40 cm soil moisture anomaly",
    "SFCSM" : "Surface soil moisture anomaly",
    "TWS" : "Terrestrial water storage anomaly",
    "Precip" : "Total precipitation amount anomaly",
    "AirT" : "Air temperature anomaly",
    "ET" : "Total evapotranspiration anomaly",
    "Streamflow" : "Streamflow anomaly",
    "SWE": "Snow water equivalent anomaly",
    "SnowDepth": "Snowdepth anomaly",
    "TOP40ST" : "Top 0-40 cm soil temperature anomaly",    
}

def merged_metric_filename(output_dir, startdate, enddate,
                                   model_forcing, domain, weekly=False):
    """Create path to merged S2S metric netCDF file."""
    def _check_filename_size(name):
        """Make sure filename does not exceed 128 characters, per Air Force
        requirement."""
        if len(os.path.basename(name)) > 128:
            print("[ERR] Output file name is too long!")
            print(f"[ERR] {os.path.basename(name)} exceeds 128 characters!")
            sys.exit(1)
            
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing.upper()}"
    name += "_GR.C0P25DEG"
    if domain == 'AFRICOM':
        name += "_AR.AFRICA"
    if domain == 'GLOBAL':
        name += "_AR.GLOBAL"
    if weekly:
        name += "_PA.S2SMETRICS_WEEKLY"
    else:
        name += "_PA.S2SMETRICS"
    name += f"_DD.{startdate.year:04d}{startdate.month:02d}01"
    name += f"_FP.{startdate.year:04d}{startdate.month:02d}{startdate.day:02d}"
    name += f"-{enddate.year:04d}{enddate.month:02d}{enddate.day:02d}"
    name += "_DF.NC"
    _check_filename_size(name)
    return name

def get_anom(path, var_name, metric, weekly=False):
    def preproc(ds_):
        ''' select only the 1st ensemble for streamflow '''
        ds_ = ds_.isel(ens=0)
        return ds_

    if weekly:
        regex = f"{path}/PS.*GP.LIS-S2S-*S2SMETRICS_WEEKLY*DF.NC"
    else:
        regex = f"{path}/PS.*GP.LIS-S2S-*S2SMETRICS_DD*DF.NC"

    files = glob.glob(regex)
    if var_name == 'Streamflow':
        anom_ds = xr.open_mfdataset(files, concat_dim='ens', combine='nested', preprocess=preproc)
    else:
        anom_ds = xr.open_mfdataset(files, concat_dim='ens', combine='nested')
    anom_ds = anom_ds.rename({'latitude': 'lat', 'longitude': 'lon'})
    anom = anom_ds[var_name + '_' + metric]
    return anom.to_dataset(name='anom')

    
    

    


