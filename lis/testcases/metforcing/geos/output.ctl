DSET ^OUTPUT/SURFACEMODEL/%y4/%y4%m2%d2/LIS_HIST_%y4%m2%d2%h2%n2.d01.gs4r
options template
options sequential
options big_endian
TITLE LIS output
UNDEF -9999.0
XDEF    1440 LINEAR  -179.875  0.25
YDEF     600 LINEAR   -59.875  0.25
ZDEF 1 LINEAR 1 1
TDEF 24 LINEAR 02Z30sep2002 1hr
VARS 7
Wind      1  99  **  1 insert description here m/s
Rainf     1  99  **  2 insert description here kg/m2s
Tair      1  99  **  3 insert description here K
Qair      1  99  **  4 insert description here kg/kg
PSurf     1  99  **  5 insert description here Pa
SWdown    1  99  **  6 insert description here W/m2
LWdown    1  99  **  7 insert description here W/m2
ENDVARS
