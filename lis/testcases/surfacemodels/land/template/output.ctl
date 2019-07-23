DSET       ^OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d2%h2%n2.d01.gs4r
TITLE       Noah3.3 output for Bondville testcase
OPTIONS     template
OPTIONS     sequential
OPTIONS     big_endian
UNDEF       -9999.
XDEF        229 LINEAR        -124.875    0.25
YDEF        109 LINEAR          25.875    0.25
ZDEF          1 LINEAR             1.0     1.0
TDEF         16 LINEAR 03:00Z29oct2002     3hr
VARS          7
Wind_f        1  99  ** Near-surface wind magnitude [m sec-1]
Rainf_f       1  99  ** Rainfall rate [kg m-2 sec-1]
Tair_f        1  99  ** Near-surface air temperature [K]
Qair_f        1  99  ** Near-surface specific humidity [kg kg-1]
PSurf_f       1  99  ** Surface pressure [Pa]
SWdown_f      1  99  ** Surface incident shortwave radiation [W m-2]
LWdown_f      1  99  ** Surface incident longwave radiation [W m-2]
ENDVARS
