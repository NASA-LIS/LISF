dset ^OUTPUT/ROUTING/%y4%m2/LIS_HIST_%y4%m2%d2%h200.d01.nc
dtype netcdf
options template
undef -9999
xdef 120 linear -89.45 0.1
ydef 85 linear 34.05 0.1
zdef 1 linear 1 1
* dummy tdef
tdef 2137 linear 00z01feb2025 1hr
vars 6
lat=>lat 1 y,x description
lon=>lon 1 y,x description
Streamflow_inst=>Streamflow_inst 1 y,x description
RiverStor_inst=>RiverStor_inst 1 y,x description
RiverDepth_inst=>RiverDepth_inst 1 y,x description
SurfElev_inst=>SurfElev_inst 1 y,x description
endvars
