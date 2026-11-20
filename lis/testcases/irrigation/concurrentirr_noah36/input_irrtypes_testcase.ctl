dset  ^TARGET_OUTPUT/lis_input.d01.nc
dtype netcdf
undef -9999
xdef 1440 linear -179.875 0.25
ydef 600 linear -59.875 0.25
zdef 1 linear 1 1
tdef 1 linear 2jan2000 1dy
vars 18
DOMAINMASK=>domainmask  0  y,x  DOMAINMASK
LANDMASK=>landmask  0  y,x  LANDMASK
ELEVFGRD=>elevfgrd  0  y,x  Elevation Area Fraction
ELEVATION=>elevation  0  y,x  GTOPO30 (Native) elevation
SHDMIN=>shdmin  0  y,x  NCEP (Native) min greenness
SHDMAX=>shdmax  0  y,x  NCEP (Native) max greenness
MXSNALBEDO=>mxsnalbedo  0  y,x  NCEP (Native) max snow albedo
TBOT=>tbot  0  y,x  Noah LSM bottom temperature
SLOPETYPE=>slopetype  0  y,x  Noah LSM slope type
IRRIGTYPE=>irrigtype1  0  0,y,x  AQUASTAT and USGS reported Sprinkler Irrigation type fraction
IRRIGTYPE=>irrigtype2  0  1,y,x  AQUASTAT and USGS reported Drip Irrigation type fraction
IRRIGTYPE=>irrigtype3  0  2,y,x  AQUASTAT and USGS reported Flood Irrigation type fraction
IRRIGATEDAREA=>irrigatedarea  0  y,x  GIA (Meier et al, 2018) Irrigated gridcell area
COUNTY=>county  0  y,x  US Census Boundary 2015 GEOID
COUNTRY=>country  0  y,x  GIS political boundary data
ELEV_PRINCETON=>elev_princeton  0  y,x  Forcing elevation for PRINCETON
lat=>lat 1 y,x description
lon=>lon 1 y,x description
endvars

