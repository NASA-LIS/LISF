dset ^OUTPUT/SURFACEMODEL/%y4/%y4%m2%d2/LIS_HIST_%y4%m2%d2%h200.d01.nc
options template
dtype netcdf
undef -9999.0
xdef 1929 linear -124.935 0.03
ydef  929 linear   25.065 0.03
*tdef 672 linear 01Z03may2013 1hr
*tdef 224 linear 03Z03may2013 3hr
tdef 28 linear 00Z04may2013 1dy
zdef 1 linear 1 1
vars 5
lat=>lat 1 y,x #(north_south, east_west)
lon=>lon 1 y,x #(north_south, east_west)
Landmask_inst=>Landmask_inst 1 y,x #(north_south, east_west)
Landcover_inst=>Landcover_inst 1 y,x #(north_south, east_west)
Greenness_inst=>gvf 1 y,x #(north_south, east_west)
endvars
