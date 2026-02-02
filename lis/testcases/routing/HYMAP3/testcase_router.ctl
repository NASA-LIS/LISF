dset ^TARGET_OUTPUT/ROUTING/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 75 linear -112.95 0.1
ydef 85 linear 35.05 0.1
zdef 1 linear 1 1
* dummy tdef
tdef 366 linear 00z01jan1990 1dy
vars 10
Streamflow_tavg=>Streamflow_tavg 1 y,x description
RiverStor_tavg=>RiverStor_tavg 1 y,x description
RiverDepth_tavg=>RiverDepth_tavg 1 y,x description
RiverFlowVelocity_tavg=>RiverFlowVelocity_tavg 1 y,x description
FloodQ_tavg=>FloodQ_tavg 1 y,x description
FloodStor_tavg=>FloodStor_tavg 1 y,x description
FloodVelocity_tavg=>FloodVelocity_tavg 1 y,x description
FloodedFrac_tavg=>FloodedFrac_tavg 1 y,x description
SurfElev_tavg=>SurfElev_tavg 1 y,x description
SWS_tavg=>SWS_tavg 1 y,x description
endvars
