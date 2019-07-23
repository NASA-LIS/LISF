dset ^TARGET_OUTPUT/RTM/%y4/%y4%m2%d2/LIS_HIST_%y4%m2%d2%h2%n2.d01.nc 
dtype netcdf
options template
UNDEF -9999.0
xdef 58 linear -124.5 1.0
ydef 28 linear 25.5 1.0
zdef 1 linear 1 1
TDEF 120 LINEAR 01z29mar2010 01hr
vars 6
RTM_emissivity_inst=>RTM_emiss_inst1 0 0,y,x rtm_emissivity_(H-pol)
RTM_emissivity_inst=>RTM_emiss_inst2 0 1,y,x rtm_emissivity_(V-pol)
RTM_Tb_inst=>RTM_Tb_inst1 0 0,y,x rtm_brightness_temperature_(H-pol)_K
RTM_Tb_inst=>RTM_Tb_inst2 0 1,y,x rtm_brightness_temperature_(V-pol)_K
lat=>lat 1 y,x latitude
lon=>lon 1 y,x longitude
endvars
