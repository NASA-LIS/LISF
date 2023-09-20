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
import numpy as np

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
    out_var = np.ones(target_fcst_data.shape)*-99.
    if fcst_std > 0.:
        out_var = (target_fcst_data - fcst_clim) / fcst_std

    return out_var
