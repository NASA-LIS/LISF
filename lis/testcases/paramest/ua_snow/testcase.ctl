dset ^TARGET_OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 355 linear -108.565 0.01
ydef 320 linear 37.525 0.01
zdef 1 linear 1 1
* dummy tdef
tdef 4 linear 0z02jan2007 1dy
vars 13
lat=>lat 1 y,x description
lon=>lon 1 y,x description
Snowf_tavg=>Snowf_tavg 1 y,x description
Rainf_tavg=>Rainf_tavg 1 y,x description
Qs_tavg=>Qs_tavg 1 y,x description
Qsb_tavg=>Qsb_tavg 1 y,x description
SWE_tavg=>SWE_tavg 1 y,x description
SnowDepth_tavg=>SnowDepth_tavg 1 y,x description
SoilMoist_tavg=>SoilMoist_tavg1 0 0,y,x description
SoilMoist_tavg=>SoilMoist_tavg2 0 1,y,x description
SoilMoist_tavg=>SoilMoist_tavg3 0 2,y,x description
SoilMoist_tavg=>SoilMoist_tavg4 0 3,y,x description
TotalPrecip_tavg=>TotalPrecip_tavg 1 y,x description
endvars
