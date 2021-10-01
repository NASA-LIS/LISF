dset ^TARGET_OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 414 linear -125.0859 0.140625
ydef 300 linear 24.98438 0.09375
zdef 1 linear 1 1
* dummy tdef
tdef 182 linear 00z01feb2016 1dy
vars 14
lat=>lat 1 y,x description
lon=>lon 1 y,x description
Evap_tavg=>Evap_tavg 1 y,x description
Qs_tavg=>Qs_tavg 1 y,x description
Qsb_tavg=>Qsb_tavg 1 y,x description
SoilMoist_tavg=>SoilMoist_tavg1 0 0,y,x description
SoilMoist_tavg=>SoilMoist_tavg2 0 1,y,x description
SoilMoist_tavg=>SoilMoist_tavg3 0 2,y,x description
SoilMoist_tavg=>SoilMoist_tavg4 0 3,y,x description
SoilMoist_inst=>SoilMoist_inst1 0 0,y,x description
SoilMoist_inst=>SoilMoist_inst2 0 1,y,x description
SoilMoist_inst=>SoilMoist_inst3 0 2,y,x description
SoilMoist_inst=>SoilMoist_inst4 0 3,y,x description
TotalPrecip_tavg=>TotalPrecip_tavg 1 y,x description
endvars
