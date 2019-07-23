DSET ^forcing.nc
DTYPE netcdf
TITLE LIS output
UNDEF -9999.0
XDEF   575 LINEAR -179.6875  0.625
YDEF   360 LINEAR  -89.75    0.5
ZDEF 1 LINEAR 1 1
TDEF 24 LINEAR 00z01nov2005 1hr
*T2M=>T2M 1 t,y,x #(time, lat, lon)
*QV2M=>QV2M 1 t,y,x #(time, lat, lon)
*U10M=>U10M 1 t,y,x #(time, lat, lon)
*V10M=>V10M 1 t,y,x #(time, lat, lon)
*PRECTOT=>PRECTOT 1 t,y,x #(time, lat, lon)
vars 14
PS=>PS 1 t,y,x #(time, lat, lon)
TLML=>TLML 1 t,y,x #(time, lat, lon)
QLML=>QLML 1 t,y,x #(time, lat, lon)
ULML=>ULML 1 t,y,x #(time, lat, lon)
VLML=>VLML 1 t,y,x #(time, lat, lon)
PRECTOTCORR=>PRECTOTCORR 1 t,y,x #(time, lat, lon)
PRECCON=>PRECCON 1 t,y,x #(time, lat, lon)
PRECSNO=>PRECSNO 1 t,y,x #(time, lat, lon)
HLML=>HLML 1 t,y,x #(time, lat, lon)
SWGDN=>SWGDN 1 t,y,x #(time, lat, lon)
LWGAB=>LWGAB 1 t,y,x #(time, lat, lon)
SWLAND=>SWLAND 1 t,y,x #(time, lat, lon)
PARDR=>PARDR 1 t,y,x #(time, lat, lon)
PARDF=>PARDF 1 t,y,x #(time, lat, lon)
endvars
