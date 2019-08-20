DSET ^OUTPUT/SURFACEMODEL/%y4/%y4%m2%d2/LIS_HIST_%y4%m2%d2%h2%n2.d01.gs4r 
options template
options sequential
options big_endian
TITLE LIS output
UNDEF -9999.0
XDEF   232 LINEAR  -124.875  0.25
YDEF   112 LINEAR    25.125  0.25
ZDEF 1 LINEAR 1 1
TDEF 48 LINEAR 02Z29oct2002 1hr
VARS 8
Wind        1  99  ** 22 insert description here m/s                
Rainf       1  99  ** 23 insert description here kg/m2s             
Snowf       1  99  ** 24 insert description here kg/m2s             
Tair        1  99  ** 25 insert description here K                  
Qair        1  99  ** 26 insert description here kg/kg              
PSurf       1  99  ** 27 insert description here Pa                 
SWdown      1  99  ** 28 insert description here W/m2               
LWdown      1  99  ** 29 insert description here W/m2               
ENDVARS
