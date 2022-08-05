dset ^OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
* ./LVT lvt.config
dtype netcdf
options template
undef -9999
xdef 1560 linear -124.995 0.01
ydef 2110 linear 31.905 0.01
zdef 1 linear 1 1
* dummy tdef
tdef 31 linear 00z02oct2017 1dy
vars 24
Qle_tavg=>Qle_tavg 1 y,x description
Qh_tavg=>Qh_tavg 1 y,x description
Evap_tavg=>Evap_tavg 1 y,x description
Qs_tavg=>Qs_tavg 1 y,x description
Qsb_tavg=>Qsb_tavg 1 y,x description
RadT_tavg=>RadT_tavg 1 y,x description
SWE_tavg=>SWE_tavg 1 y,x description
SnowDepth_tavg=>SnowDepth_tavg 1 y,x description
SoilMoist_tavg=>SoilMoist_tavg1 0 0,y,x description
SoilMoist_tavg=>SoilMoist_tavg2 0 1,y,x description
SoilMoist_tavg=>SoilMoist_tavg3 0 2,y,x description
SoilMoist_tavg=>SoilMoist_tavg4 0 3,y,x description
TVeg_tavg=>TVeg_tavg 1 y,x description
ESoil_tavg=>ESoil_tavg 1 y,x description
TWS_tavg=>TWS_tavg 1 y,x description
GWS_tavg=>GWS_tavg 1 y,x description
GPP_tavg=>GPP_tavg 1 y,x description
Wind_f_tavg=>Wind_f_tavg 1 y,x description
Tair_f_tavg=>Tair_f_tavg 1 y,x description
Qair_f_tavg=>Qair_f_tavg 1 y,x description
Psurf_f_tavg=>Psurf_f_tavg 1 y,x description
LAI_tavg=>LAI_tavg 1 y,x description
RHMin_tavg=>RHMin_tavg 1 y,x description
TotalPrecip_tavg=>TotalPrecip_tavg 1 y,x description
endvars
