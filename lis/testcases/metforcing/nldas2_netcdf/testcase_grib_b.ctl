dset ^TARGET_OUTPUT/grib_forb/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d2%h2%n2.d01.nc
dtype netcdf
options template
undef -9999
xdef 464 linear -124.9375 0.125
ydef 224 linear 25.0625 0.125
zdef 1 linear 1 1
* dummy tdef
tdef 24 linear 14:15Z02jan1979 15mn
vars 10
Wind_f_tavg=>Wind_f_tavg 1 y,x description
Rainf_f_tavg=>Rainf_f_tavg 1 y,x description
CRainf_f_tavg=>CRainf_f_tavg 1 y,x description
Tair_f_tavg=>Tair_f_tavg 1 y,x description
Qair_f_tavg=>Qair_f_tavg 1 y,x description
Psurf_f_tavg=>Psurf_f_tavg 1 y,x description
SWdown_f_tavg=>SWdown_f_tavg 1 y,x description
LWdown_f_tavg=>LWdown_f_tavg 1 y,x description
FHeight_f_tavg=>FHeight_f_tavg 1 y,x description
Ch_f_tavg=>Ch_f_tavg 1 y,x description
endvars
