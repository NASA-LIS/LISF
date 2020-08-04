dset ^lis_input.d01.nc
dtype netcdf
options template
undef -9999
* Lambert conformal:
*pdef 215 110 lcc 37.50 -120.20 1 1 37 38 -120 1000 1000 
*xdef 215 linear -121.20 0.02
*ydef 110 linear   37.00 0.02
* Latlon grid-equivalent:
xdef 215 linear -121.20 0.01
ydef 110 linear   37.50 0.01
zdef 1 linear 1 1
* dummy tdef
tdef 1 linear 00z01jan2012 1hr
vars 110
LANDMASK=>LANDMASK 1 y,x description
DOMAINMASK=>DOMAINMASK 1 y,x description
SURFACETYPE=>SURFACETYPE1 0 0,y,x description
SURFACETYPE=>SURFACETYPE2 0 1,y,x description
SURFACETYPE=>SURFACETYPE3 0 2,y,x description
SURFACETYPE=>SURFACETYPE4 0 3,y,x description
SURFACETYPE=>SURFACETYPE5 0 4,y,x description
SURFACETYPE=>SURFACETYPE6 0 5,y,x description
SURFACETYPE=>SURFACETYPE7 0 6,y,x description
SURFACETYPE=>SURFACETYPE8 0 7,y,x description
SURFACETYPE=>SURFACETYPE9 0 8,y,x description
SURFACETYPE=>SURFACETYPE10 0 9,y,x description
SURFACETYPE=>SURFACETYPE11 0 10,y,x description
SURFACETYPE=>SURFACETYPE12 0 11,y,x description
SURFACETYPE=>SURFACETYPE13 0 12,y,x description
SURFACETYPE=>SURFACETYPE14 0 13,y,x description
SURFACETYPE=>SURFACETYPE15 0 14,y,x description
SURFACETYPE=>SURFACETYPE16 0 15,y,x description
SURFACETYPE=>SURFACETYPE17 0 16,y,x description
SURFACETYPE=>SURFACETYPE18 0 17,y,x description
SURFACETYPE=>SURFACETYPE19 0 18,y,x description
SURFACETYPE=>SURFACETYPE20 0 19,y,x description
LANDCOVER=>LANDCOVER1 0 0,y,x description
LANDCOVER=>LANDCOVER2 0 1,y,x description
LANDCOVER=>LANDCOVER3 0 2,y,x description
LANDCOVER=>LANDCOVER4 0 3,y,x description
LANDCOVER=>LANDCOVER5 0 4,y,x description
LANDCOVER=>LANDCOVER6 0 5,y,x description
LANDCOVER=>LANDCOVER7 0 6,y,x description
LANDCOVER=>LANDCOVER8 0 7,y,x description
LANDCOVER=>LANDCOVER9 0 8,y,x description
LANDCOVER=>LANDCOVER10 0 9,y,x description
LANDCOVER=>LANDCOVER11 0 10,y,x description
LANDCOVER=>LANDCOVER12 0 11,y,x description
LANDCOVER=>LANDCOVER13 0 12,y,x description
LANDCOVER=>LANDCOVER14 0 13,y,x description
LANDCOVER=>LANDCOVER15 0 14,y,x description
LANDCOVER=>LANDCOVER16 0 15,y,x description
LANDCOVER=>LANDCOVER17 0 16,y,x description
LANDCOVER=>LANDCOVER18 0 17,y,x description
LANDCOVER=>LANDCOVER19 0 18,y,x description
LANDCOVER=>LANDCOVER20 0 19,y,x description
TEXTURE=>TEXTURE1 0 0,y,x description
TEXTURE=>TEXTURE2 0 1,y,x description
TEXTURE=>TEXTURE3 0 2,y,x description
TEXTURE=>TEXTURE4 0 3,y,x description
TEXTURE=>TEXTURE5 0 4,y,x description
TEXTURE=>TEXTURE6 0 5,y,x description
TEXTURE=>TEXTURE7 0 6,y,x description
TEXTURE=>TEXTURE8 0 7,y,x description
TEXTURE=>TEXTURE9 0 8,y,x description
TEXTURE=>TEXTURE10 0 9,y,x description
TEXTURE=>TEXTURE11 0 10,y,x description
TEXTURE=>TEXTURE12 0 11,y,x description
TEXTURE=>TEXTURE13 0 12,y,x description
TEXTURE=>TEXTURE14 0 13,y,x description
TEXTURE=>TEXTURE15 0 14,y,x description
TEXTURE=>TEXTURE16 0 15,y,x description
ELEVFGRD=>ELEVFGRD 1 y,x description
ELEVATION=>ELEVATION 1 y,x description
ASPECTFGRD=>ASPECTFGRD 1 y,x description
ASPECT=>ASPECT 1 y,x description
SLOPEFGRD=>SLOPEFGRD 1 y,x description
SLOPE=>SLOPE 1 y,x description
GREENNESS=>GREENNESS1 0 0,y,x description
GREENNESS=>GREENNESS2 0 1,y,x description
GREENNESS=>GREENNESS3 0 2,y,x description
GREENNESS=>GREENNESS4 0 3,y,x description
GREENNESS=>GREENNESS5 0 4,y,x description
GREENNESS=>GREENNESS6 0 5,y,x description
GREENNESS=>GREENNESS7 0 6,y,x description
GREENNESS=>GREENNESS8 0 7,y,x description
GREENNESS=>GREENNESS9 0 8,y,x description
GREENNESS=>GREENNESS10 0 9,y,x description
GREENNESS=>GREENNESS11 0 10,y,x description
GREENNESS=>GREENNESS12 0 11,y,x description
SHDMIN=>SHDMIN 1 y,x description
SHDMAX=>SHDMAX 1 y,x description
ALBEDO=>ALBEDO1 0 0,y,x description
ALBEDO=>ALBEDO2 0 1,y,x description
ALBEDO=>ALBEDO3 0 2,y,x description
ALBEDO=>ALBEDO4 0 3,y,x description
ALBEDO=>ALBEDO5 0 4,y,x description
ALBEDO=>ALBEDO6 0 5,y,x description
ALBEDO=>ALBEDO7 0 6,y,x description
ALBEDO=>ALBEDO8 0 7,y,x description
ALBEDO=>ALBEDO9 0 8,y,x description
ALBEDO=>ALBEDO10 0 9,y,x description
ALBEDO=>ALBEDO11 0 10,y,x description
ALBEDO=>ALBEDO12 0 11,y,x description
MXSNALBEDO=>MXSNALBEDO 1 y,x description
TBOT=>TBOT 1 y,x description
SLOPETYPE=>SLOPETYPE 1 y,x description
NOAHMP36_PBLH=>NOAHMP36PBLH 1 y,x description
PPTCLIM=>PPTCLIM1 0 0,y,x description
PPTCLIM=>PPTCLIM2 0 1,y,x description
PPTCLIM=>PPTCLIM3 0 2,y,x description
PPTCLIM=>PPTCLIM4 0 3,y,x description
PPTCLIM=>PPTCLIM5 0 4,y,x description
PPTCLIM=>PPTCLIM6 0 5,y,x description
PPTCLIM=>PPTCLIM7 0 6,y,x description
PPTCLIM=>PPTCLIM8 0 7,y,x description
PPTCLIM=>PPTCLIM9 0 8,y,x description
PPTCLIM=>PPTCLIM10 0 9,y,x description
PPTCLIM=>PPTCLIM11 0 10,y,x description
PPTCLIM=>PPTCLIM12 0 11,y,x description
ELEV_NLDAS2=>ELEVNLDAS2 0 y,x description
ELEVDIFF_NLDAS2=>ELEVDIFFNLDAS2 0 y,x description
lat=>lat 1 y,x description
lon=>lon 1 y,x description
endvars
