dset ^TARGET_OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 464 linear -124.9375 0.125
ydef 224 linear 25.0625 0.125
zdef 1 linear 1 1
* dummy tdef
tdef 59 linear 00z02jan1979 1dy
vars 30
lat=>lat 1 y,x description
lon=>lon 1 y,x description
Qle_tavg=>Qle_tavg 1 y,x description
Qh_tavg=>Qh_tavg 1 y,x description
Evap_acc=>Evap_acc 1 y,x description
Qs_acc=>Qs_acc 1 y,x description
Qsb_acc=>Qsb_acc 1 y,x description
Qsm_acc=>Qsm_acc 1 y,x description
SWE_inst=>SWE_inst 1 y,x description
SnowDepth_inst=>SnowDepth_inst 1 y,x description
SoilMoist_tavg=>SoilMoist_tavg1 0 0,y,x description
SoilMoist_tavg=>SoilMoist_tavg2 0 1,y,x description
SoilMoist_tavg=>SoilMoist_tavg3 0 2,y,x description
SoilMoist_tavg=>SoilMoist_tavg4 0 3,y,x description
SoilMoist_inst=>SoilMoist_inst1 0 0,y,x description
SoilMoist_inst=>SoilMoist_inst2 0 1,y,x description
SoilMoist_inst=>SoilMoist_inst3 0 2,y,x description
SoilMoist_inst=>SoilMoist_inst4 0 3,y,x description
SoilTemp_tavg=>SoilTemp_tavg1 0 0,y,x description
SoilTemp_tavg=>SoilTemp_tavg2 0 1,y,x description
SoilTemp_tavg=>SoilTemp_tavg3 0 2,y,x description
SoilTemp_tavg=>SoilTemp_tavg4 0 3,y,x description
ECanop_acc=>ECanop_acc 1 y,x description
WaterTableD_tavg=>WaterTableD_tavg 1 y,x description
Snowcover_inst=>Snowcover_inst 1 y,x description
GPP_tavg=>GPP_tavg 1 y,x description
NPP_tavg=>NPP_tavg 1 y,x description
SWdown_f_inst=>SWdown_f_inst 1 y,x description
LAI_tavg=>LAI_tavg 1 y,x description
TotalPrecip_acc=>TotalPrecip_acc 1 y,x description
endvars
