DSET ^TARGET_OUTPUT/SURFACEMODEL/202312/LIS_HIST_%y4%m2%d2%h2%n2.d01.nc
DTYPE netcdf
OPTIONS template
TITLE LIS output
UNDEF -9999.0
XDEF 101 LINEAR -122.5 0.01
YDEF 101 LINEAR 36.0 0.01
ZDEF 1 LINEAR 1 1
TDEF 6 LINEAR 01:00Z31dec2023 1hr
VARS 5
Tair_f_tavg       1 99 ** air temperature K
Qair_f_tavg       1 99 ** specific humidity kg/kg
Psurf_f_tavg      1 99 ** surface pressure Pa
LWdown_f_tavg     1 99 ** downward longwave radiation W/m2
Lapserate_tavg    1 99 ** lapse rate K/km
ENDVARS
