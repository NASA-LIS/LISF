'''
This script computes median across all model ensembles for each lead month separately and compares against the S2S produced TIF file.
'''

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

import sys
import argparse
from datetime import date
from dateutil.relativedelta import relativedelta
import xarray as xr
import numpy as np
import yaml

def compute_median (anom, lead):
    da_slice=[]
    for da in anom:
        da_slice.append(da.isel(time=lead))
    da_conc = xr.concat(da_slice, dim = 'ens')
    return da_conc.median(dim = 'ens')

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--fcast_month', required=True, help='forecast start month 1-12')
    parser.add_argument('-y', '--fcast_year', required=True, help='forecast start year')
    parser.add_argument('-c', '--config_file', required=True, help='config file')
    args = parser.parse_args()
    fcst_year = int(args.fcast_year)
    fcst_mon = int(args.fcast_month)
    with open(args.config_file, 'r', encoding="utf-8") as file:
        cfg = yaml.safe_load(file)
    sys.path.append(cfg['SETUP']['LISFDIR'] + '/lis/utils/usaf/s2s/')
    from s2s_modules.shared import utils
    
    end_date = date(fcst_year, fcst_mon, 1) + relativedelta(months=int(cfg["EXP"]["lead_months"]))
    s2smdir = './s2smetric/' + '{:04d}{:02d}'.format(fcst_year, fcst_mon) + '/metrics_cf/NOAHMP/'
    d0 = date(fcst_year, fcst_mon, 1)
    
    for var in cfg["POST"]["metric_vars"]:
        print ("Comparing : ", var)
        
        # read NC files
        anoms = []
        sanoms = []
        print (s2smdir + 'PS.557WW_SC.U_DI.C_GP.LIS-S2S-*_GR.C0P25DEG_AR.GLOBAL_PA.S2SMETRICS_DD.{:04d}{:02d}01_FP.{:04d}{:02d}01-{}_DF.NC'.format(fcst_year, fcst_mon, fcst_year, fcst_mon, end_date.strftime("%Y%m%d")))
        for model in cfg["EXP"]["NMME_models"]:
            ncfile = s2smdir + 'PS.557WW_SC.U_DI.C_GP.LIS-S2S-{}_GR.C0P25DEG_AR.GLOBAL_PA.S2SMETRICS_DD.{:04d}{:02d}01_FP.{:04d}{:02d}01-{}_DF.NC'.format(model.upper(), fcst_year, fcst_mon, fcst_year, fcst_mon, end_date.strftime("%Y%m%d")) 
            ncdata = xr.open_dataset(ncfile)
            anoms.append(ncdata[var + '_ANOM'])
            sanoms.append(ncdata[var + '_SANOM'])
            del ncdata
    
        # read TIF files
        for lead in range(0, int(cfg["EXP"]["lead_months"])):        
            d1 = date(fcst_year, fcst_mon, 1)  + relativedelta(months=lead)
            d2 = d1 + relativedelta(months=1)
            tif_anom =  s2smdir + 'PS.557WW_SC.U_DI.C_GP.LIS-S2S-ANOM_GR.C0P25DEG_AR.GLOBAL_PA.{}_DD.{}_FP.{}-{}_DF.TIF'.format(var.upper(),d0.strftime("%Y%m%d"),d1.strftime("%Y%m%d"),d2.strftime("%Y%m%d"))
            tif_sanom =  s2smdir + 'PS.557WW_SC.U_DI.C_GP.LIS-S2S-SANOM_GR.C0P25DEG_AR.GLOBAL_PA.{}_DD.{}_FP.{}-{}_DF.TIF'.format(var.upper(),d0.strftime("%Y%m%d"),d1.strftime("%Y%m%d"),d2.strftime("%Y%m%d"))
            da_anom = utils.tiff_to_da(tif_anom)
            da_sanom = utils.tiff_to_da(tif_sanom)
        
            anom_rev = da_anom.reindex(y=list(reversed(da_anom.y)))
            sanom_rev = da_sanom.reindex(y=list(reversed(da_sanom.y)))

            nc_med = compute_median (anoms, lead)
            c1 = np.allclose(nc_med.values, anom_rev.values, rtol=1.e-05,   equal_nan=True)
            nc_med = compute_median (sanoms, lead)
            c2 = np.allclose(nc_med.values, sanom_rev.values, rtol=1.e-05,   equal_nan=True)
            print('LEAD : ', lead)
            print ('PS.557WW_SC.U_DI.C_GP.LIS-S2S-ANOM_GR.C0P25DEG_AR.GLOBAL_PA.{}_DD.{}_FP.{}-{}_DF.TIF'.format(var.upper(),d0.strftime("%Y%m%d"),d1.strftime("%Y%m%d"),d2.strftime("%Y%m%d")), '  : ', c1)
            print ('PS.557WW_SC.U_DI.C_GP.LIS-S2S-SANOM_GR.C0P25DEG_AR.GLOBAL_PA.{}_DD.{}_FP.{}-{}_DF.TIF'.format(var.upper(),d0.strftime("%Y%m%d"),d1.strftime("%Y%m%d"),d2.strftime("%Y%m%d")), '  : ', c2)
 


