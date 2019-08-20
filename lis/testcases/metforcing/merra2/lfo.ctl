*DSET ^input/MET_FORCING/MERRA2/%y4%m2%d2/MERRA2_300.tavg1_2d_lfo_Nx.%y4%m2%d2.nc4
DSET ^input/MET_FORCING/MERRA2/MERRA2_300/stage/Y%y4/M%m2/MERRA2_300.tavg1_2d_lfo_Nx.%y4%m2%d2.nc4
DTYPE netcdf
options template
TITLE LIS output
UNDEF 1e+15
XDEF   576 LINEAR -180.0  0.625
YDEF   361 LINEAR  -90.0  0.5
ZDEF 1 LINEAR 1 1
TDEF 24 LINEAR 00z01nov2005 1hr
vars 8
LWGAB=>LWGAB 1 t,y,x #(time, lat, lon)
PARDF=>PARDF 1 t,y,x #(time, lat, lon)
PARDR=>PARDR 1 t,y,x #(time, lat, lon)
PRECCUCORR=>PRECCUCORR 1 t,y,x #(time, lat, lon)
PRECLSCORR=>PRECLSCORR 1 t,y,x #(time, lat, lon)
PRECSNOCORR=>PRECSNOCORR 1 t,y,x #(time, lat, lon)
SWGDN=>SWGDN 1 t,y,x #(time, lat, lon)
SWLAND=>SWLAND 1 t,y,x #(time, lat, lon)
endvars
