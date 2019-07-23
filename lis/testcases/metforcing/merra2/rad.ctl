*DSET ^input/MET_FORCING/MERRA2/%y4%m2%d2/MERRA2_300.tavg1_2d_rad_Nx.%y4%m2%d2.nc4
DSET ^input/MET_FORCING/MERRA2/MERRA2_300/stage/Y%y4/M%m2/MERRA2_300.tavg1_2d_rad_Nx.%y4%m2%d2.nc4
DTYPE netcdf
options template
TITLE LIS output
UNDEF 1e+15
XDEF   576 LINEAR -180.0  0.625
YDEF   361 LINEAR  -90.0  0.5
ZDEF 1 LINEAR 1 1
TDEF 24 LINEAR 00z01nov2005 1hr
vars 36
ALBEDO=>ALBEDO            1  t,y,x #(time, lat, lon)
ALBNIRDF=>ALBNIRDF        1  t,y,x #(time, lat, lon)
ALBNIRDR=>ALBNIRDR        1  t,y,x #(time, lat, lon)
ALBVISDF=>ALBVISDF        1  t,y,x #(time, lat, lon)
ALBVISDR=>ALBVISDR        1  t,y,x #(time, lat, lon)
CLDHGH=>CLDHGH            1  t,y,x #(time, lat, lon)
CLDLOW=>CLDLOW            1  t,y,x #(time, lat, lon)
CLDMID=>CLDMID            1  t,y,x #(time, lat, lon)
CLDTOT=>CLDTOT            1  t,y,x #(time, lat, lon)
EMIS=>EMIS                1  t,y,x #(time, lat, lon)
LWGAB=>LWGAB              1  t,y,x #(time, lat, lon)
LWGABCLR=>LWGABCLR        1  t,y,x #(time, lat, lon)
LWGABCLRCLN=>LWGABCLRCLN  1  t,y,x #(time, lat, lon)
LWGEM=>LWGEM              1  t,y,x #(time, lat, lon)
LWGNT=>LWGNT              1  t,y,x #(time, lat, lon)
LWGNTCLR=>LWGNTCLR        1  t,y,x #(time, lat, lon)
LWGNTCLRCLN=>LWGNTCLRCLN  1  t,y,x #(time, lat, lon)
LWTUP=>LWTUP              1  t,y,x #(time, lat, lon)
LWTUPCLR=>LWTUPCLR        1  t,y,x #(time, lat, lon)
LWTUPCLRCLN=>LWTUPCLRCLN  1  t,y,x #(time, lat, lon)
SWGDN=>SWGDN              1  t,y,x #(time, lat, lon)
SWGDNCLR=>SWGDNCLR        1  t,y,x #(time, lat, lon)
SWGNT=>SWGNT              1  t,y,x #(time, lat, lon)
SWGNTCLN=>SWGNTCLN        1  t,y,x #(time, lat, lon)
SWGNTCLR=>SWGNTCLR        1  t,y,x #(time, lat, lon)
SWGNTCLRCLN=>SWGNTCLRCLN  1  t,y,x #(time, lat, lon)
SWTDN=>SWTDN              1  t,y,x #(time, lat, lon)
SWTNT=>SWTNT              1  t,y,x #(time, lat, lon)
SWTNTCLN=>SWTNTCLN        1  t,y,x #(time, lat, lon)
SWTNTCLR=>SWTNTCLR        1  t,y,x #(time, lat, lon)
SWTNTCLRCLN=>SWTNTCLRCLN  1  t,y,x #(time, lat, lon)
TAUHGH=>TAUHGH            1  t,y,x #(time, lat, lon)
TAULOW=>TAULOW            1  t,y,x #(time, lat, lon)
TAUMID=>TAUMID            1  t,y,x #(time, lat, lon)
TAUTOT=>TAUTOT            1  t,y,x #(time, lat, lon)
TS=>TS                    1  t,y,x #(time, lat, lon)
endvars
