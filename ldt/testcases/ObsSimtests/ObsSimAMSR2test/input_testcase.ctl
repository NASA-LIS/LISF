dset ^TARGET_OUTPUT/%y4%m2/SimObs_%y4%m2%d20000.nc
dtype netcdf
options template
undef -9999
xdef 13 linear -108.445 0.25
ydef 12 linear 37.645 0.25
zdef 1 linear 1 1
* dummy tdef
tdef 29 linear 00z02feb2020 1dy
vars 4
lat=>lat 1 y,x description
lon=>lon 1 y,x description
SWE_tavg=>SWE_TAVG 1 y,x description
SWE_inst=>SWE_INST 1 y,x description
endvars
