DSET ^OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.gs4r 
options template
options sequential
options big_endian
TITLE LIS output
UNDEF -9999.0
xdef 81 linear -95 0.25
ydef 81 linear 20 0.25
zdef 1 linear 1 1
TDEF 4 LINEAR 00z02jan2012 01dy
VARS 7
Wind        1  99  ** 22 insert description here m/s                
Tair        1  99  ** 25 insert description here K                  
Qair        1  99  ** 26 insert description here kg/kg              
PSurf       1  99  ** 27 insert description here Pa                 
SWdown      1  99  ** 28 insert description here W/m2               
LWdown      1  99  ** 29 insert description here W/m2               
Elevation   1  99  ** 29 insert description here W/m2               
ENDVARS
