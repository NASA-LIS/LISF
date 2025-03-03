dset ^OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d2%h200.d01.nc
dtype netcdf
options template
undef -9999
xdef 250 linear -114.95 0.1
ydef 135 linear 36.55 0.1
zdef 1 linear 1 1
tdef 8 linear 03z02jan2022 3hr
vars 22
lat=>lat 1 y,x description
lon=>lon 1 y,x description
Qle_tavg=>Qle_tavg 1 y,x description
Qh_tavg=>Qh_tavg 1 y,x description
Evap_tavg=>Evap_tavg 1 y,x description
Qs_tavg=>Qs_tavg 1 y,x description
Qsb_tavg=>Qsb_tavg 1 y,x description
SWE_tavg=>SWE_tavg 1 y,x description
SnowDepth_tavg=>SnowDepth_tavg 1 y,x description
SoilMoist_tavg=>SoilMoist_tavg1 0 0,y,x description
SoilMoist_tavg=>SoilMoist_tavg2 0 1,y,x description
SoilMoist_tavg=>SoilMoist_tavg3 0 2,y,x description
SoilMoist_tavg=>SoilMoist_tavg4 0 3,y,x description
Wind_f_tavg=>Wind_f_tavg 1 y,x description
Rainf_f_tavg=>Rainf_f_tavg 1 y,x description
Tair_f_tavg=>Tair_f_tavg 1 y,x description
Qair_f_tavg=>Qair_f_tavg 1 y,x description
Psurf_f_tavg=>Psurf_f_tavg 1 y,x description
SWdown_f_tavg=>SWdown_f_tavg 1 y,x description
LWdown_f_tavg=>LWdown_f_tavg 1 y,x description
LAI_inst=>LAI_inst 1 y,x description
TotalPrecip_tavg=>TotalPrecip_tavg 1 y,x description
endvars
