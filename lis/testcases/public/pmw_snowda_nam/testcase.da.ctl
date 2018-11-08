DSET ^OUTPUT/da_testcase/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d2%h2%n2.d01.nc
DTYPE netcdf 
options template
TITLE LIS output
UNDEF -9999.0
XDEF   860 LINEAR -172.975  0.05
YDEF   515 LINEAR   50.025  0.05
ZDEF 1 LINEAR 1 1
TDEF 91 LINEAR 00z2mar2011 1dy
VARS 5         
Qs_tavg=>qs     1  y,x surface runoff (mm/s)                  
Qsb_tavg=>qsb   1  y,x subsurface runoff (mm/s)              
SWE_inst=>swe   1  y,x snow water equivalent (mm)                 
Snowcover_inst=>snowcover       1  y,x snow cover fration (%)               
SnowDepth_inst=>snowdepth       1  y,x snow depth (m)               
ENDVARS
