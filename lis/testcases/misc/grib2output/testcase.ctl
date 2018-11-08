dset ^TARGET_OUTPUT/SURFACEMODEL/%y4/%y4%m2%d2/LIS_HIST_%y4%m2%d2%h2%n2.d01.gr2
index ^testcase.idx
options template
undef 9.999E+20
title OUTPUT/testcase/SURFACEMODEL/2002/20021029/LIS_HIST_200210290300.d01.gr2
* produced by g2ctl v0.1.2
* command line options: OUTPUT/testcase/SURFACEMODEL/2002/20021029/LIS_HIST_200210290300.d01.gr2
* griddef=1:0:(229 x 109):grid_template=0:winds(N/S): lat-lon grid:(229 x 109) units 1e-06 input WE:SN output WE:SN res 48 lat 25.875000 to 52.875000 by 0.250000 lon 235.125000 to 292.125000 by 0.250000 #points=24961:winds(N/S)
*  gribmap -v -0 -i testcase.ctl
*
* units specification is disc,cat,num,sp,sp2
*    for instantaneous fields specify disc,cat,num
*    for time-averaged fields specify disc,cat,num,0

dtype grib2
ydef 109 linear 25.875000 0.25
xdef 229 linear -124.875000 0.250000
tdef 16 linear 03:00Z29oct2002 3hr
zdef 1 linear 1 1
vars 33
Tairfsfc   0,1,0   0,0,0 ** surface desc [unit]
Qlesfc   0,1,0   0,0,10,0 ** surface desc [unit]
Qhsfc   0,1,0   0,0,11,0 ** surface desc [unit]
Emissfsfc   0,1,0   0,0,124 ** surface desc [unit]
Albedosfc   0,1,0   0,19,1 ** surface desc [unit]
Qairfsfc   0,1,0   0,1,0 ** surface desc [unit]
SnowDepthsfc   0,1,0   0,1,11 ** surface desc [unit]
Snowcoversfc   0,1,0   0,1,42 ** surface desc [unit]
Rainffsfc   0,1,0   0,1,52 ** surface desc [unit]
SWEsfc   0,1,0   0,1,60 ** surface desc [unit]
Rainfsfc   0,1,0   0,1,65,0 ** surface desc [unit]
Snowfsfc   0,1,0   0,1,66,0 ** surface desc [unit]
Evapsfc   0,1,0   0,1,79,0 ** surface desc [unit]
Windfsfc   0,1,0   0,2,1 ** surface desc [unit]
Psurffsfc   0,1,0   0,3,0 ** surface desc [unit]
SWdownfsfc   0,1,0   0,4,7 ** surface desc [unit]
Swnetsfc   0,1,0   0,4,9,0 ** surface desc [unit]
LWdownfsfc   0,1,0   0,5,3 ** surface desc [unit]
Lwnetsfc   0,1,0   0,5,5,0 ** surface desc [unit]
Qsbsfc   0,1,0   1,0,5,0 ** surface desc [unit]
Qssfc   0,1,0   1,0,6,0 ** surface desc [unit]
CanopIntsfc   0,1,0   2,0,13 ** surface desc [unit]
Qgsfc   0,1,0   2,0,24,0 ** surface desc [unit]
Greennesssfc   0,1,0   2,0,87 ** surface desc [unit]
SoilTemp0_10cm  0,106,0,0.1   2,3,18 ** 0-0.1 m below ground desc [unit]
SoilTemp10_40cm  0,106,0.1,0.4   2,3,18 ** 0.1-0.4 m below ground desc [unit]
SoilTemp40_100cm  0,106,0.4,1   2,3,18 ** 0.4-1 m below ground desc [unit]
SoilTemp100_200cm  0,106,1,2   2,3,18 ** 1-2 m below ground desc [unit]
SoilMoist0_10cm  0,106,0,0.1   2,3,19 ** 0-0.1 m below ground desc [unit]
SoilMoist10_40cm  0,106,0.1,0.4   2,3,19 ** 0.1-0.4 m below ground desc [unit]
SoilMoist40_100cm  0,106,0.4,1   2,3,19 ** 0.4-1 m below ground desc [unit]
SoilMoist100_200cm  0,106,1,2   2,3,19 ** 1-2 m below ground desc [unit]
AvgSurfTsfc   0,1,0   2,3,201 ** surface desc [unit]
ENDVARS
