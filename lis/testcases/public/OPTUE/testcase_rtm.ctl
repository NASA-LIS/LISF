DSET ^TARGET_OUTPUT/RTM/LIS_HIST_%y4%m2%d2%h200.d01.nc
DTYPE   netcdf
OPTIONS template
UNDEF -9999.0
*XDEF east_west    20 LINEAR     -99.875   0.25
*YDEF north_south  20 LINEAR      34.125   0.25
*ZDEF RTM_emissivity_profiles  12 LEVELS  1 2 3 4 5 6 7 8 9 10 11 12
*TDEF time       2184 LINEAR  01Z01APR2008  1hr
XDEF 20 LINEAR     -99.875   0.25
YDEF 20 LINEAR      34.125   0.25
ZDEF 1 LINEAR 1 1
TDEF 2184 LINEAR  01Z01APR2008  1hr
VARS 14
lat=>lat                             1    y,x Latitude  [-]
lon=>lon                             1    y,x Longitude [-]
RTM_emissivity_inst=>RTM_emis_inst1  0  0,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst2  0  1,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst3  0  2,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst4  0  3,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst5  0  4,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst6  0  5,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst7  0  6,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst8  0  7,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst9  0  8,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst10 0  9,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst11 0 10,y,x RTM emissivity
RTM_emissivity_inst=>RTM_emis_inst12 0 11,y,x RTM emissivity
ENDVARS
