DSET ^TARGET_OUTPUT/SURFACEMODEL/%y4%m2/LIS_HIST_%y4%m2%d20000.d01.nc
options template
dtype netcdf
TITLE LIS output
UNDEF -9999.0
xdef 891 linear -115.925 0.01
ydef 987 linear 29.065 0.01
zdef 1 linear 1 1
TDEF 35 LINEAR 00z02oct2011 01dy
VARS 3
lat=>lat 1 y,x latitude
lon=>lon 1 y,x longitude
Rainf_f_tavg=>Rainf_f_tavg 1 y,x rain mm/s
ENDVARS
