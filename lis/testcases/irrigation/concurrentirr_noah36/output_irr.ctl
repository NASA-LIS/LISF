dset ^OUTPUT/IRRIGATION/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 1440 linear -179.875 0.25
ydef 600 linear -59.875 0.25
zdef 1 linear 1 1
tdef 9 linear 2jan2000 1dy
*** SDF_name(15 characters)=>grads_name
vars 13
IrrigatedWater_tavg=>irrigatedwater_  0  y,x  irrigated water amount
IrrigWaterInput_tavg=>irrigwaterinput  0  y,x  input irrigation water amount
IrrigationScale_tavg=>irrigationscale  0  y,x  Scaled irrigation fraction to landcover types
IrrigtypeSp_tavg=>irrigtypesp_tav  0  y,x  Sprinkler irrigation water amount
IrrigtypeDr_tavg=>irrigtypedr_tav  0  y,x  Drip irrigation water amount
IrrigtypeFl_tavg=>irrigtypefl_tav  0  y,x  Flood irrigation water amount
IrrigPaddy_tavg=>irrigpaddy_tavg  0  y,x  irrigated paddy fraction
totR_tavg=>totr_tavg  0  y,x  total subsurface runoff
recR_tavg=>recr_tavg  0  y,x  recycled subsurface runoff
outR_tavg=>outr_tavg  0  y,x  net output subsurface runoff
QRI_tavg=>qri_tavg  0  y,x  runoff recycled to irrigation
lat=>lat 1 y,x description
lon=>lon 1 y,x description
endvars
