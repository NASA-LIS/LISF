DSET ^OUTPUT/SURFACEMODEL/202312/LIS_HIST_%y4%m2%d2%h2%n2.d01.nc
dtype netcdf
options template
title LIS output
undef -9999
xdef 101 linear -122.5 0.01
ydef 101 linear 36. 0.01
zdef 1 linear 1 1
tdef 6 LINEAR 01:00Z31dec2023 1hr
vars 5
Tair_f_tavg=>Tair_f_tavg 1 y,x description
Qair_f_tavg=>Qair_f_tavg 1 y,x description
Psurf_f_tavg=>Psurf_f_tavg 1 y,x description
LWdown_f_tavg=>LWdown_f_tavg 1 y,x description
Lapserate_tavg=>Lapserate_tavg 1 y,x description
endvars
