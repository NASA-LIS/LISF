dset ^OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 36 linear -5.875 0.25
ydef 40 linear 5.125 0.25
zdef 1 linear 1 1
tdef 2 linear 00z03may2001 1dy
vars 3
lat=>lat 1 y,x description
lon=>lon 1 y,x description
LAI_tavg=>LAI_tavg 1 y,x description
endvars
