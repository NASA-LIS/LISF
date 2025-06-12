import os
import sys
import numpy as np
import  xarray as xr
import xesmf as xe
import argparse
import sys
import yaml
from ghis2s.shared.utils import get_domain_info
from ghis2s.bcsd.bcsd_library import convert_forecast_data_to_netcdf as cfdn

CFSV2_file = '/discover/nobackup/projects/lis/MET_FORCING/CFSv2//Oper_TS/2024/20241217/tmp2m.01.2024121706.daily.grb2'
def get_land_mask(config, cfsv2_data):
    """
    Create the CFSv2 land mask using the LDT LANDMASK
    We use xESMF's conservative_normed interpolation method 
    to ensure any LDT resolution will have an active land grid cell on the VFSv2 grid.
    """
    # Get the CFSv2 coordinates
    cfsv2_lat = cfsv2_data.latitude.values
    cfsv2_lon = cfsv2_data.longitude.values
    
    # Open the LDT file and read the land mask
    ldt_xr = xr.open_dataset(config['SETUP']['supplementarydir'] + '/lis_darun/' + config['SETUP']['ldtinputfile'])
    ldt_xr = ldt_xr.rename({'north_south': 'latitude', 'east_west': 'longitude'})
    ldt_xr = ldt_xr.assign_coords(
        longitude=np.round(ldt_xr.lon.values[0, :], 6),
        latitude=np.round(ldt_xr.lat.values[:, 0], 6))

    # Extract the land mask (0=ocean, 1=land)
    ldt_mask_da = ldt_xr['LANDMASK']
    
    if ldt_mask_da.longitude.min() < 0 and cfsv2_lon.min() >= 0:
        print("Converting LDT longitudes from [-180, 180] to [0, 360]")
        # Convert longitudes using xarray's functionality
        ldt_mask_da = ldt_mask_da.assign_coords(longitude=((ldt_mask_da.longitude + 360) % 360)).sortby('longitude')
        print(f"Converted longitudes: {ldt_mask_da.longitude.min().values:.2f} to {ldt_mask_da.longitude.max().values:.2f}")
    
    # da to ds
    ldt_mask_ds = xr.Dataset({"LANDMASK": ldt_mask_da})

    # target CFSv2 grid
    cfsv2_grid = xr.Dataset({
        "latitude": (["latitude"], cfsv2_lat),
        "longitude": (["longitude"], cfsv2_lon),
    })
    
    mask_regridder = xe.Regridder(ldt_mask_ds, cfsv2_grid, "conservative_normed", periodic=True)
    regridded_landmask = mask_regridder(ldt_mask_ds.LANDMASK)
    
    # Create output dataset
    regridded_ds = xr.Dataset({"LANDMASK": (regridded_landmask > 0).astype(np.uint8)})
    print(f"Regridded mask range: {regridded_ds.LANDMASK.min().values:.4f} to {regridded_ds.LANDMASK.max().values:.4f}")
    print(f"Land fraction (>0): {(regridded_ds.LANDMASK > 0).sum().values / regridded_ds.LANDMASK.size:.4f}")
    
    return regridded_ds

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
    land_mask = get_land_mask(config, cfsv2)
    land_mask.to_netcdf(weightdir + f'CFSv2_{resol}km_landmask.nc4', format="NETCDF4",
                        encoding = {"LANDMASK": {"zlib":True, "complevel":6, "shuffle":True, "missing_value":-9999.}})

    #any_land = land_mask.LANDMASK > 0
    #
    #sample_field = cfsv2.T2M.isel(step=0)
    #masked_sample = sample_field.where(any_land)

    #masked_cfsv2 = xr.Dataset({
    #    "T2M": masked_sample,  
    #})

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
    

    
