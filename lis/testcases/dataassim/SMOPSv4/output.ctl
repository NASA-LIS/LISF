dset ^OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d2%h200.d01.nc
dtype netcdf
options template
undef -9999
xdef 2560 linear -179.9297 0.140625
ydef 1920 linear -89.95312 0.09375
zdef 1 linear 1 1
* dummy tdef
tdef 3 linear 00z01may2024 12hr
vars 10
lat=>lat 1 y,x description
lon=>lon 1 y,x description
SoilMoist_tavg=>SoilMoist_tavg1 0 0,y,x description
SoilMoist_tavg=>SoilMoist_tavg2 0 1,y,x description
SoilMoist_tavg=>SoilMoist_tavg3 0 2,y,x description
SoilMoist_tavg=>SoilMoist_tavg4 0 3,y,x description
SoilMoist_inst=>SoilMoist_inst1 0 0,y,x description
SoilMoist_inst=>SoilMoist_inst2 0 1,y,x description
SoilMoist_inst=>SoilMoist_inst3 0 2,y,x description
SoilMoist_inst=>SoilMoist_inst4 0 3,y,x description
endvars
