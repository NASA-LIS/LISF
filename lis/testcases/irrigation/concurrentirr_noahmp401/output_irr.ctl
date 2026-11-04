dset ^OUTPUT/IRRIGATION/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
dtype netcdf
options template
undef -9999
xdef 1440 linear -179.875 0.25
ydef 600 linear -59.875 0.25
zdef 1 linear 1 1
tdef 4 linear 00z02jan2000 1dy
vars 13
lat=>lat 1 y,x description
lon=>lon 1 y,x description
IrrigatedWater_tavg=>IrrigatedWater_tavg 1 y,x description
IrrigWaterInput_tavg=>IrrigWaterInput_tavg 1 y,x description
IrrigationScale_tavg=>IrrigationScale_tavg 1 y,x description
IrrigtypeSp_tavg=>IrrigtypeSp_tavg 1 y,x description
IrrigtypeDr_tavg=>IrrigtypeDr_tavg 1 y,x description
IrrigtypeFl_tavg=>IrrigtypeFl_tavg 1 y,x description
IrrigPaddy_tavg=>IrrigPaddy_tavg 1 y,x description
totR_tavg=>totR_tavg 1 y,x description
recR_tavg=>recR_tavg 1 y,x description
outR_tavg=>outR_tavg 1 y,x description
QRI_tavg=>QRI_tavg 1 y,x description
endvars
