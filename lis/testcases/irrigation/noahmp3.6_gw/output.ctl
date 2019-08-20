DSET ^Sprinkler/GWon/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
options template
dtype netcdf
TITLE LIS output
UNDEF -9999.0
xdef 464 linear -124.9375 0.125
ydef 224 linear 25.0625 0.125
zdef 1 linear 1 1
TDEF 35 LINEAR 00z03jun2001 01dy
VARS 34
lat=>lat 1 y,x latitude
lon=>lon 1 y,x longitude
TWS_tavg=>TWS_tavg 1 y,x TWS mm
SWnet         1  99  ** Net shortwave radiation [W m-2]
LWnet         1  99  ** Net longwave radiation [W m-2]
Qle           1  99  ** Latent heat flux [W m-2]
Qh            1  99  ** Sensible heat flux [W m-2]
Qg            1  99  ** Ground heat flux [W m-2]
Snowf         1  99  ** Snowfall rate [kg m-2 sec-1]
Rainf         1  99  ** Rainfall rate [kg m-2 sec-1]
Evap          1  99  ** Total evapotranspiration [kg m-2 sec-1]
Qs            1  99  ** Surface runoff [kg m-2 sec-1]
Qsb           1  99  ** Subsurface runoff [kg m-2 sec-1]
AvgSurfT      1  99  ** Average surface temperature [K]
Albedo        1  99  ** Surface albedo [0-1]
SWE           1  99  ** Snow water equivalent [kg m-2]
SnowDepth     1  99  ** Snow depth [m]
SoilMoist1    1  99  ** Soil water content layer 1 [kg m-2]
SoilMoist2    1  99  ** Soil water content layer 2 [kg m-2]
SoilMoist3    1  99  ** Soil water content layer 3 [kg m-2]
SoilMoist4    1  99  ** Soil water content layer 4 [kg m-2]
SoilTemp1     1  99  ** Soil temperature layer 1 [K]
SoilTemp2     1  99  ** Soil temperature layer 2 [K]
SoilTemp3     1  99  ** Soil temperature layer 3 [K]
SoilTemp4     1  99  ** Soil temperature layer 4 [K]
CanopInt      1  99  ** Canopy interception [kg m-2]
Snowcover     1  99  ** Snow cover fraction [0-1]
Wind_f        1  99  ** Near-surface wind magnitude [m sec-1]
Rainf_f       1  99  ** Rainfall rate [kg m-2 sec-1]
Tair_f        1  99  ** Near-surface air temperature [K]
Qair_f        1  99  ** Near-surface specific humidity [kg kg-1]
PSurf_f       1  99  ** Surface pressure [Pa]
SWdown_f      1  99  ** Surface incident shortwave radiation [W m-2]
LWdown_f      1  99  ** Surface incident longwave radiation [W m-2]
Greenness     1  99  ** Green vegetation (shade) fraction [-]
ENDVARS
