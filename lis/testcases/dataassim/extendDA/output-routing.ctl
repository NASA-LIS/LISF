dset ^OUTPUT/ROUTING/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 36 linear -5.875 0.25
ydef 40 linear 5.125 0.25
zdef 1 linear 1 1
tdef 2 linear 00z03may2001 1dy
vars 8
lat=>lat 1 y,x description
lon=>lon 1 y,x description
Streamflow_tavg=>Streamflow_tavg 1 y,x description
RiverStor_tavg=>RiverStor_tavg 1 y,x description
RiverDepth_tavg=>RiverDepth_tavg 1 y,x description
RiverFlowVelocity_tavg=>RiverFlowVelocity_tavg 1 y,x description
FloodStor_tavg=>FloodStor_tavg 1 y,x description
SurfElev_tavg=>SurfElev_tavg 1 y,x description
endvars
