import os
import sys
import  xarray as xr
import xesmf as xe
import argparse
import yaml
from ghis2s.shared.utils import get_domain_info
from ghis2s.bcsd.bcsd_library import convert_forecast_data_to_netcdf as cfdn

CFSV2_file = '/discover/nobackup/projects/lis/MET_FORCING/CFSv2//Oper_TS/2024/20241217/tmp2m.01.2024121706.daily.grb2'

if __name__ == "__main__":
    # Parse command arguements
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    args = parser.parse_args()

    with open(args.config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    weightdir = config['SETUP']['supplementarydir'] + '/bcsd_fcst/'
    lats, lons = get_domain_info(args.config_file, coord=True)
    resol = round((lats[1] - lats[0])*100)
    cfsv2 = cfdn.wgrib2_to_netcdf(CFSV2_file)
    
    # bilinear
    weight_file = weightdir + f'CFSv2_{resol}km_bilinear.nc'
    print(weight_file)
    ds_out = xr.Dataset(
            {
                "lat": (["lat"], lats),
                "lon": (["lon"], lons),
            }
        )
    run_xe= xe.Regridder(cfsv2, ds_out, "bilinear", periodic=True, 
                         filename=weight_file)
    # conservative
    weight_file = weightdir + f'CFSv2_{resol}km_conservative.nc'
    print(weight_file)
    ds_out = xr.Dataset(
            {
                "lat": (["lat"], lats),
                "lon": (["lon"], lons),
            }
        )
    run_xe= xe.Regridder(cfsv2, ds_out, "conservative", periodic=True, 
                         filename=weight_file)
    

    
