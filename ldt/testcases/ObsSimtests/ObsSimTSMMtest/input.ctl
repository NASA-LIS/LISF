dset ^OUTPUT/%y4%m2/SimObs_%y4%m2%d20000.nc
dtype netcdf
options template
undef -9999
xdef 355 linear -108.565 0.01
ydef 320 linear 37.525 0.01
zdef 1 linear 1 1
* dummy tdef
tdef 29 linear 00z02feb2020 1dy
vars 4
lat=>lat 1 y,x description
lon=>lon 1 y,x description
SWE_tavg=>SWE_TAVG 1 y,x description
SWE_inst=>SWE_INST 1 y,x description
endvars
