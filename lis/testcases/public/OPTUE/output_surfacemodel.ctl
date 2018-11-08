DSET ^OUTPUT/SURFACEMODEL/LIS_HIST_%y4%m2%d2%h200.d01.nc
DTYPE   netcdf
OPTIONS template
UNDEF -9999.0
*XDEF east_west    20 LINEAR     -99.875   0.25
*YDEF north_south  20 LINEAR      34.125   0.25
*ZDEF SoilMoist_profiles  4 LEVELS  1 2 3 4
*TDEF time       2184 LINEAR  01Z01APR2008  1hr
*
* To view the soil temperature profile, comment out
*  the 'ZDEF SoilMoist_profiles ...' line above and
*  uncomment the line below:
*ZDEF SoilTemp_profiles  4 LEVELS  1 2 3 4
*
XDEF 20 LINEAR     -99.875   0.25
YDEF 20 LINEAR      34.125   0.25
ZDEF 4 LEVELS  1 2 3 4
TDEF 2184 LINEAR  01Z01APR2008  1hr
VARS 17
lat=>lat                         1    y,x lat
lon=>lon                         1    y,x lon
Rainf_tavg=>Rainf_tavg           1    y,x Rainf_tavg
AvgSurfT_tavg=>AvgSurfT_tavg     1    y,x AvgSurfT_tavg
SnowDepth_tavg=>SnowDepth_tavg   1    y,x SnowDepth_tavg
SoilMoist_inst=>SoilMoist_inst1  0  0,y,x SoilMoist_inst
SoilMoist_inst=>SoilMoist_inst2  0  1,y,x SoilMoist_inst
SoilMoist_inst=>SoilMoist_inst3  0  2,y,x SoilMoist_inst
SoilMoist_inst=>SoilMoist_inst4  0  3,y,x SoilMoist_inst
SoilTemp_inst=>SoilTemp_inst1    0  0,y,x SoilTemp_inst
SoilTemp_inst=>SoilTemp_inst2    0  1,y,x SoilTemp_inst
SoilTemp_inst=>SoilTemp_inst3    0  2,y,x SoilTemp_inst
SoilTemp_inst=>SoilTemp_inst4    0  3,y,x SoilTemp_inst
CanopInt_inst=>CanopInt_inst     1    y,x CanopInt_inst
Rainf_f_inst=>Rainf_f_inst       1    y,x Rainf_f_inst
LAI_inst=>LAI_inst               1    y,x LAI_inst
Greenness_inst=>Greenness_inst   1    y,x Greenness_inst
ENDVARS
