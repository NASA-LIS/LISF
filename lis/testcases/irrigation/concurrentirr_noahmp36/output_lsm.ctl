dset ^OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 1440 linear -179.875 0.25
ydef 600 linear -59.875 0.25
zdef 1 linear 1 1
tdef 4 linear 00z02jan2000 1dy
vars 42
lat=>lat 1 y,x description
lon=>lon 1 y,x description
Swnet_tavg=>Swnet_tavg 1 y,x description
Lwnet_tavg=>Lwnet_tavg 1 y,x description
Qle_tavg=>Qle_tavg 1 y,x description
Qh_tavg=>Qh_tavg 1 y,x description
Qg_tavg=>Qg_tavg 1 y,x description
Snowf_tavg=>Snowf_tavg 1 y,x description
Rainf_tavg=>Rainf_tavg 1 y,x description
Evap_tavg=>Evap_tavg 1 y,x description
Qs_tavg=>Qs_tavg 1 y,x description
Qsb_tavg=>Qsb_tavg 1 y,x description
Qsm_tavg=>Qsm_tavg 1 y,x description
AvgSurfT_inst=>AvgSurfT_inst 1 y,x description
Albedo_inst=>Albedo_inst 1 y,x description
SWE_inst=>SWE_inst 1 y,x description
SnowDepth_inst=>SnowDepth_inst 1 y,x description
SoilMoist_tavg=>SoilMoist_tavg1 0 0,y,x description
SoilMoist_tavg=>SoilMoist_tavg2 0 1,y,x description
SoilMoist_tavg=>SoilMoist_tavg3 0 2,y,x description
SoilMoist_tavg=>SoilMoist_tavg4 0 3,y,x description
SoilTemp_inst=>SoilTemp_inst1 0 0,y,x description
SoilTemp_inst=>SoilTemp_inst2 0 1,y,x description
SoilTemp_inst=>SoilTemp_inst3 0 2,y,x description
SoilTemp_inst=>SoilTemp_inst4 0 3,y,x description
SmLiqFrac_tavg=>SmLiqFrac_tavg1 0 0,y,x description
SmLiqFrac_tavg=>SmLiqFrac_tavg2 0 1,y,x description
SmLiqFrac_tavg=>SmLiqFrac_tavg3 0 2,y,x description
SmLiqFrac_tavg=>SmLiqFrac_tavg4 0 3,y,x description
ECanop_tavg=>ECanop_tavg 1 y,x description
TVeg_tavg=>TVeg_tavg 1 y,x description
ESoil_tavg=>ESoil_tavg 1 y,x description
CanopInt_inst=>CanopInt_inst 1 y,x description
SubSnow_tavg=>SubSnow_tavg 1 y,x description
Snowcover_inst=>Snowcover_inst 1 y,x description
Wind_f_tavg=>Wind_f_tavg 1 y,x description
Rainf_f_tavg=>Rainf_f_tavg 1 y,x description
Tair_f_tavg=>Tair_f_tavg 1 y,x description
Qair_f_tavg=>Qair_f_tavg 1 y,x description
Psurf_f_tavg=>Psurf_f_tavg 1 y,x description
SWdown_f_tavg=>SWdown_f_tavg 1 y,x description
LWdown_f_tavg=>LWdown_f_tavg 1 y,x description
endvars
