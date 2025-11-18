import os
import sys
import numpy as np
import xarray as xr
import xesmf as xe
import yaml
from ghis2s.shared.utils import get_domain_info
from ghis2s.bcsd.bcsd_library.bcsd_function import apply_regridding_with_mask

'''
setenv PYTHONPATH /discover/nobackup/projects/usaf_lis/smahanam/S2S/LISF-1/lis/utils/usaf/S2S/
obs_in
netcdf LWGAB_obs_clim {
dimensions:
        DIST = 13 ;
        time = 30 ;
        latitude = 720 ;
        longitude = 1440 ;
variables:
        double clim(DIST, time, latitude, longitude) ;
                clim:_FillValue = NaN ;
                clim:missing_value = NaN ;
        int64 DIST(DIST) ;
        int64 time(time) ;
        float latitude(latitude) ;
                latitude:_FillValue = NaNf ;
        float longitude(longitude) ;
                longitude:_FillValue = NaNf ;
}

fcst_in
netcdf LWS_fcst_clim {
dimensions:
        DIST = 10 ;
        time = 360 ;
        latitude = 720 ;
        longitude = 1440 ;
variables:
        double clim(DIST, time, latitude, longitude) ;
                clim:_FillValue = NaN ;
                clim:missing_value = NaN ;
        int64 DIST(DIST) ;
        int64 time(time) ;
        float latitude(latitude) ;
                latitude:_FillValue = NaNf ;
        float longitude(longitude) ;
                longitude:_FillValue = NaNf ;
}

'''

CMDARGS = str(sys.argv)

resol = str(sys.argv[1])
ldtin_file = '/discover/nobackup/projects/usaf_lis/GHI_S2S/supplementary_files/lis_darun/lis_input.s2s_global.noahmp401_hymap2.25km.nc'
ldtout_file = f'/discover/nobackup/projects/usaf_lis/GHI_S2S/supplementary_files/lis_darun/lis_input.s2s_global.noahmp401.merit.{resol}_test.nc'
config_file = f'/discover/nobackup/projects/ghilis/smahanam/E2ES_{resol}/s2s_config_global_fcast'
USAF_clim_in = '/discover/nobackup/projects/ghilis/S2S/GLOBAL/E2ES_rerun_08012023/hindcast/bcsd_fcst/USAF-LIS7.3rc8_25km/raw/Climatology/'
USAF_clim_out = f'/discover/nobackup/projects/ghilis/smahanam/E2ES_{resol}/hindcast/bcsd_fcst/USAF-LIS7.3rc8_{resol}/raw/Climatology/'
#CFSv2_clim_in = '/discover/nobackup/projects/ghilis/S2S/GLOBAL/E2ES_rerun_08012023/hindcast/bcsd_fcst/CFSv2_25km/raw/Climatology/mar01/'
#CFSv2_clim_out = f'/discover/nobackup/projects/ghilis/S2S/GLOBAL/E2ES_{resol}/hindcast/bcsd_fcst/CFSv2_{resol}/raw/Climatology/mar01/'

fcst_vars = ['LWS', 'PRECTOT', 'PS', 'Q2M', 'SLRSF', 'T2M', 'WIND10M']
obs_vars = ['LWGAB', 'PRECTOT', 'PS', 'QV2M', 'SWGDN', 'T2M', 'U10M']
#obs_vars = ['LWGAB', 'PRECTOT', 'PS', 'QV2M', 'SWGDN', 'T2M', 'WIND10M']
fcst_template = '{}/{}_fcst_clim.nc'
obs_template = '{}/{}_obs_clim.nc'

ldt_xr = xr.open_dataset(ldtin_file)
ldt_xr = ldt_xr.rename({'north_south': 'latitude', 'east_west': 'longitude'})
ldtout_xr = xr.open_dataset(ldtout_file)
ldtout_xr = ldtout_xr.rename({'north_south': 'latitude', 'east_west': 'longitude'})

# get lat/lon of the output grid
LATS, LONS = get_domain_info(config_file, coord=True)

ds_out = xr.Dataset({
    "latitude": (["latitude"], LATS),
    "longitude": (["longitude"], LONS),
})

# fcst_var regridding - First file (LWS)
#infile = fcst_template.format(CFSv2_clim_in, 'LWS')
#outfile = fcst_template.format(CFSv2_clim_out, 'LWS')

#ds_in = xr.open_dataset(infile)
#ds_in['mask'] = ldt_xr['LANDMASK']

#regridder = xe.Regridder(ds_in, ds_out, "conservative", periodic=True)
#result = regridder(ds_in)  # Pass ds_in, not ds_out

#result.to_netcdf(outfile, format="NETCDF4", encoding={'clim': encoding})

# Process remaining forecast variables
#for v in range(1, len(fcst_vars)):
#    infile = fcst_template.format(CFSv2_clim_in, fcst_vars[v])
#    outfile = fcst_template.format(CFSv2_clim_out, fcst_vars[v])
    
#    ds_in = xr.open_dataset(infile)
#    ds_in['mask'] = ldt_xr['LANDMASK']  # Add mask for each file
    
    # Create new regridder for each input (if grids differ) or reuse if same
#    regridder = xe.Regridder(ds_in, ds_out, "conservative", periodic=True)
#    result = regridder(ds_in)  # Pass ds_in, not ds_out
#    result.to_netcdf(outfile, format="NETCDF4", encoding={'clim': encoding})

# Similar pattern for obs_vars if needed
encoding = {"zlib":True, "complevel":6, "shuffle":True, "missing_value": np.nan, "_FillValue": np.nan}
for v in range(len(obs_vars)):
    infile = obs_template.format(USAF_clim_in, obs_vars[v])  # Assuming USAF for obs
    outfile = obs_template.format(USAF_clim_out, obs_vars[v])
    
    ds_in = xr.open_dataset(infile)
    ds_in['mask'] = ldt_xr['LANDMASK']
    ds_out['mask'] = ldtout_xr['LANDMASK']
    regridder = xe.Regridder(ds_in, ds_out, "nearest_s2d", periodic=True)
    result = regridder(ds_in)
    result.to_netcdf(outfile, format="NETCDF4", encoding={'clim': encoding})
    
    
    
    

