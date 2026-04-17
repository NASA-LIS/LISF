''' This script is good to create weight files SD '''
import argparse
import numpy as np
import  xarray as xr
import xesmf as xe
import yaml
from ghis2s.shared.utils import get_domain_info
from ghis2s.bcsd.bcsd_library import convert_forecast_data_to_netcdf as cfdn

CFSV2_file = '/discover/nobackup/projects/lis/MET_FORCING/CFSv2//Oper_TS/2024/' \
    '20241217/tmp2m.01.2024121706.daily.grb2'
NMME_file = '/discover/nobackup/projects/usaf_lis/GHI_S2S/NMME//CanSIPS-IC4/' \
    'prec.CanESM5.mon_Jan.2025.nc'
GEOS5_file = '/discover/nobackup/projects/usaf_lis/smahanam/MET_FORCING/GEOSv3/' \
    'sfc_tavg_3hr_glo_L720x361_sfc/200210/ens01/geos_s2s_v3.200210.nc'

def create_geos_source_mask_from_tsurf(geos_ds):
    """
    Create a source grid mask for xESMF using FRLAND from GEOS surface file.
    Cells with FRLAND > 0 are marked as land (mask=1).
    """
    print(f"Creating GEOS source mask from FRLAND in: {GEOS5_file}")

    # Load GEOS surface file
    #geos_ds = xr.open_dataset(GEOS5_file)

    # Get FRLAND (fraction of land)
    tsurf = geos_ds['TSURF'].isel(time=0)

    print("\nTSURF statistics:")
    print(f"  Min: {tsurf.min().values:.6f}")
    print(f"  Max: {tsurf.max().values:.6f}")
    print(f"  Mean: {tsurf.mean().values:.6f}")

    # Create binary land mask: 1 where TSURF > 200 and tsurf < 500, else 0
    land_mask = ((tsurf > 200) & (tsurf < 500)).astype(np.uint8)

    print("\nMask statistics:")
    print(f"  Total cells: {land_mask.size}")
    print(f"  Land cells (TSURF>200 and TSURF< 500): {land_mask.sum().values}")
    print(f"  Ocean cells (TSURF=0): {(land_mask == 0).sum().values}")
    print(f"  Land fraction: {land_mask.sum().values / land_mask.size:.4f}")

    # Check the problematic meridian
    if 0.0 in geos_ds['lon'].values:
        lon_0_mask = land_mask.sel(lon=0.0).values
        print("\nMask at lon=0.0°:")
        print(f"  Total: {len(lon_0_mask)}")
        print(f"  Land: {lon_0_mask.sum()}")
        print(f"  Ocean: {(lon_0_mask == 0).sum()}")

        # Show latitude range where there's land at lon=0
        land_lats = geos_ds['lat'].values[lon_0_mask == 1]
        if len(land_lats) > 0:
            print(f"  Land latitude range: {land_lats.min():.1f}° to {land_lats.max():.1f}°")

    # Create dataset with CLEAN coordinates
    mask_ds = xr.Dataset(
        {
            "LANDMASK": (["lat", "lon"], land_mask.values,
                        {'long_name': 'GEOS land mask from TSURF',
                         'units': '1',
                         'description': '1 where TSURF>0 (land), 0 where TSURF=0 (ocean)'})
        },
        coords={
            "lat": (["lat"], geos_ds['lat'].values, 
                   {'units': 'degrees_north', 'long_name': 'latitude'}),
            "lon": (["lon"], geos_ds['lon'].values, 
                   {'units': 'degrees_east', 'long_name': 'longitude'})
        }
    )

    return mask_ds

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
    ldt_xr = xr.open_dataset(config['SETUP']['supplementarydir'] + '/lis_darun/' +
                             config['SETUP']['ldtinputfile'])
    ldt_xr = ldt_xr.rename({'north_south': 'latitude', 'east_west': 'longitude'})
    ldt_xr = ldt_xr.assign_coords(
        longitude=np.round(ldt_xr.lon.values[0, :], 6),
        latitude=np.round(ldt_xr.lat.values[:, 0], 6))

    # Extract the land mask (0=ocean, 1=land)
    ldt_mask_da = ldt_xr['LANDMASK']

    if ldt_mask_da.longitude.min() < 0 and cfsv2_lon.min() >= 0:
        print("Converting LDT longitudes from [-180, 180] to [0, 360]")
        # Convert longitudes using xarray's functionality
        ldt_mask_da = \
            ldt_mask_da.assign_coords(longitude=
                                      ((ldt_mask_da.longitude + 360) % 360)).sortby('longitude')

    # da to ds
    ldt_mask_ds = xr.Dataset({"LANDMASK": ldt_mask_da})

    # target CFSv2 grid
    cfsv2_grid = xr.Dataset({
        "latitude": (["latitude"], cfsv2_lat, {'units': 'degrees_north'}),
        "longitude": (["longitude"], cfsv2_lon, {'units': 'degrees_east'}),
    })

    mask_regridder = xe.Regridder(ldt_mask_ds, cfsv2_grid, "conservative_normed", periodic=True)
    regridded_landmask = mask_regridder(ldt_mask_ds.LANDMASK)

    # Create output dataset
    regridded_ds = xr.Dataset(
        {"LANDMASK": (regridded_landmask > 0).astype(np.uint8)},
        coords={
            "latitude": (["latitude"], cfsv2_lat, {'units': 'degrees_north',
                                                   'long_name': 'latitude'}),
            "longitude": (["longitude"], cfsv2_lon, {'units': 'degrees_east',
                                                     'long_name': 'longitude'})
        }
    )

    return regridded_ds

if __name__ == "__main__":
    # Parse command arguements
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, help='config file name')
    parser.add_argument('-f', '--forcing', required=True, help='config file name')
    args = parser.parse_args()

    with open(args.config_file, 'r', encoding="utf-8") as file:
        config = yaml.safe_load(file)

    weightdir = config['SETUP']['supplementarydir'] + '/bcsd_fcst/'
    lats, lons = get_domain_info(args.config_file, coord=True)
    resol = round((lats[1] - lats[0])*100)

    if args.forcing == 'NMME':
        #NMME
        force = xr.open_dataset(NMME_file, decode_times=False)
        force = force.rename({'Y': 'latitude', 'X': 'longitude'})
        land_mask = get_land_mask(config, force)

    if args.forcing == 'CFSv2':
        # CFSv2
        force = cfdn.wgrib2_to_netcdf(CFSV2_file)
        land_mask = get_land_mask(config, force)

    if args.forcing == 'GEOSv3':
        # GEOSv3
        force = xr.open_dataset(GEOS5_file)
        force = force.rename({'lat': 'latitude', 'lon': 'longitude'})
        #land_mask = create_geos_source_mask_from_tsurf(force)
        land_mask = get_land_mask(config, force)
        land_mask = land_mask.rename({'latitude': 'lat', 'longitude': 'lon'})

    land_mask.to_netcdf(weightdir + f'{args.forcing}_{resol}km_landmask.nc4', format="NETCDF4",
                        encoding = {"LANDMASK": {'dtype': 'uint8', "zlib":True, "complevel":6, "shuffle":True, "missing_value":255}})

    # bilinear_land
    force_masked = force.copy()
    force_masked['mask'] = land_mask['LANDMASK']
    weight_file = weightdir + f'{args.forcing}_{resol}km_bilinear_land.nc'
    print(weight_file)
    ds_out = xr.Dataset(
            {
                "lat": (["lat"], lats),
                "lon": (["lon"], lons),
            }
        )
    run_xe= xe.Regridder(force_masked, ds_out, "bilinear", periodic=True,
                         filename=weight_file, extrap_method='nearest_s2d')
    # conservative without masking
    weight_file = weightdir + f'{args.forcing}_{resol}km_conservative.nc'
    print(weight_file)
    ds_out = xr.Dataset(
            {
                "lat": (["lat"], lats),
                "lon": (["lon"], lons),
            }
        )
    run_xe= xe.Regridder(force, ds_out, "conservative", periodic=True,
                         filename=weight_file)
