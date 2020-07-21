!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_PRIV_rcMod 
!BOP
!
! !MODULE: LIS_PRIV_rcMod
!
! !DESCRIPTION:
!  
!  Module for specifying model independent variables in LIS. 
!  The module does not contain any variables that are specific to the 
!  extensible components in LIS. The module specifies variables 
!  for overall run control, options for the choice of parameter datasets, 
!  data assimilation choices, and the variables controlling the time
!  management in LIS. 
!
!  The variables specified in this module include: 
!  \begin{description}
!  \item[runmode]
!   choice of running mode in LIS 
!  \item[nnest]
!   number of nests or running instances (1 or higher)
!  \item[ntiles]
!   This array stores the size of tilespace of the running domain,
!   for each processor, for each nest
!  \item[glbntiles]
!   This array stores the size of the tilespace of the overall running domain 
!   for each nest (including the halo regions)
!  \item[glbntiles\_red]
!   This array stores the size of the tilespace of the overall running domain 
!   for each nest (excluding the halo regions)
!  \item[ngrid]
!   This array stores the size of gridspace of the running domain, 
!   for each processor, for each nest
!  \item[glbngrid]
!   This array stores the size of the gridspace of the overall running domain
!   for  each nest 
!  \item[glbngrid\_red]
!   This array stores the size of the gridspace of the overall running domain 
!   for each nest (excluding the halo regions)
!  \item[gnc]
!   Array containing the East-West grid dimension of the overall running 
!   grid for each nest
!  \item[gnr]
!   Array containing the North-South grid dimension of the overall running 
!   grid for each nest
!  \item[lnc]
!   Array containing the East-West grid dimension of the running 
!   grid for each processor, for each nest (including the halo regions)
!  \item[lnr]
!   Array containing the North-South grid dimension of the running 
!   grid for each processor, for each nest (including the halo regions)
!  \item[lnc\_red]
!   Array containing the East-West grid dimension of the running 
!   grid for each processor, for each nest (excluding the halo regions)
!  \item[lnr\_red]
!   Array containing the North-South grid dimension of the running 
!   grid for each processor, for each nest (excluding the halo regions)
!  \item[ncatg]
!  Total number of land tiles in the catchment-based parameter files per nest.
!  \item[lis\_map\_proj]
!   Choice of map projection used in LIS. 
!  \item[lsm]
!   Choice of the land surface model in LIS. 
!  \item[param\_proj]
!   Map projection type of parameter datasets. 
!  \item[texturemap]
!   boolean value to denote if soil texture map should be used (0- do not
!   use, 1-use). If this value is set to 1, the sand, clay and silt fraction
!   maps will not be read. Otherwise the texture data is read the fraction
!   data will not be used. 
!  \item[vegsrc]
!   Choice of the source of landcover data source
!  \item[soilsrc]
!   Choice of the source of soil parameter data source
!  \item[colorsrc]
!   Choice of the soil color data source(0-do not use)
!  \item[elevsrc]
!   Choice of the elevation data source (0-do not use)
!  \item[slopesrc]
!   Choice of the slope data source (0-do not use)
!  \item[aspectsrc]
!   Choice of the aspect data source (0-do not use)
!  \item[curvsrc]
!   Choice of the curvature data source (0-do not use)
!  \item[albedosrc]
!   Choice of the albedo data source (0-do not use)
!  \item[porositysrc]
!   Choice of the porosity data source (0-do not use)
!  \item[psisatsrc]
!   Choice of the saturated matric potential data source (0-do not use)
!  \item[ksatsrc]
!   Choice of the saturated hydraulic conductivity data source (0-do not use)
!  \item[bexpsrc]
!   Choice of the b parmeter data source (0-do not use)
!  \item[quartzsrc]
!   Choice of the quartz data source (0-do not use)
!  \item[laisrc]
!   Choice of the LAI data source (0-do not use)
!  \item[saisrc]
!   Choice of the SAI data source (0-do not use)
!  \item[tbotsrc]
!   Choice of the TBOT data source (0-do not use)
!  \item[shdmaxsrc]
!   Choice of the SHDMAX data source (0-do not use)
!  \item[shdminsrc]
!   Choice of the SHDMIN data source (0-do not use)
!  \item[slopetypesrc]
!   Choice of the slope type data source (0-do not use)
!  \item[snowsrc]
!   whether snow data sources are being stored (0 - do not use)
!  \item[surface\_maxt]
!   Maximum number of tiles per grid to be considered for subgrid tiling
!  \item[surface\_minp]
!   Minumum cutoff percentage of vegetation distribution for subgrid
!   tiling 
!  \item[decompose\_by\_processes]
!   Logical flag indicating whether to decompose the running domain based
!   on the number of processes (.true.) or to decompose based on a specified
!   layout (.false.)
!  \item[npesx]
!   Number of processors in the East-West dimensions for processor layout
!  \item[npesy]
!   Number of processors in the North-South dimensions for processor layout
!  \item[halox]
!   Halo size (in grid points) in the East-West dimensions for processor layout
!  \item[haloy]
!   Halo size (in grid points) in the North-South dimensions for processor layout
!  \item[udef]
!   Value to be used as undefined variable
!  \item[gridDesc]
!  Array describing the running grid and the parameter grid specification, for each nest.
!  The size of the array is 50 for a particular nest. The first 30 locations of this 
!  array are used to specify the running grid and the rest are used to specify parameter
!  grid. Note that though run domains can be specified in several map projections, 
!  specification of the parameter domains are currently only allowed for lat/lon, 
!  gaussian, and UTM projections. For e.g., if the user specifies a running domain in lambert, 
!  the parameter domain is expected to be in one of the supported grids for parameter
!  domain specification. The array indices for different map projections are defined as follows: 
!
!   \begin{description}
!   \item[Lat/lon or geographic projection] 
!    gridDesc(1) = 0 ; indicates lat/lon projection \newline
!    gridDesc(2) = number of columns in the domain ; \newline
!    gridDesc(3) = number of rows in the domain ; \newline
!    gridDesc(4) = latitude of the lower left corner of the domain \newline
!    gridDesc(5) = longitude of the lower left corner of the domain \newline
!    gridDesc(6) = 128 ; not used \newline
!    gridDesc(7) =  latitude of the upper right corner of the domain \newline
!    gridDesc(8) =  longitude of the upper right corner of the domain \newline
!    gridDesc(9) =  spatial resolution (in degrees) along the E-W dimension \newline
!    gridDesc(10) = spatial resolution (in degrees) along the N-S dimension  \newline
!    gridDesc(11) = 64 ; not used  \newline
!    gridDesc(20) = 255 ; used to specify the ordering of data (non-divisible by 
!                  32 indicates E-W ordering else N-S ordering \newline
!    gridDesc(30) = 0 ; indicates lat/lon projection
!    gridDesc(32) = number of columns in the domain ; \newline
!    gridDesc(33) = number of rows in the domain ; \newline
!    gridDesc(34) = latitude of the lower left corner of the domain \newline
!    gridDesc(35) = longitude of the lower left corner of the domain \newline
!    gridDesc(36) = 128 ; not used \newline
!    gridDesc(37) =  latitude of the upper right corner of the domain \newline
!    gridDesc(38) =  longitude of the upper right corner of the domain \newline
!    gridDesc(39) =  spatial resolution (in degrees) along the E-W dimension \newline
!    gridDesc(40) = spatial resolution (in degrees) along the N-S dimension  \newline
! 
!   \item[Mercator projection]
!    gridDesc(1) = 1 ; indicates mercator projection \newline
!    gridDesc(2) = number of columns in the domain ; \newline
!    gridDesc(3) = number of rows in the domain ; \newline
!    gridDesc(4) = latitude of the lower left corner of the domain \newline
!    gridDesc(5) = longitude of the lower left corner of the domain \newline
!    gridDesc(6) = 8 ; not used \newline
!    gridDesc(7) =  0.0; not used \newline
!    gridDesc(8) =  spatial resolution (in km) along E-W dimension \newline
!    gridDesc(9) =  spatial reolution (in km) along N-S dimension \newline
!    gridDesc(10) = true latitude 1 \newline
!    gridDesc(11) = standard longitude  \newline
!    gridDesc(20) = 255 ; used to specify the ordering of data (non-divisible by 
!                  32 indicates E-W ordering else N-S ordering \newline
!    gridDesc(30) = TBD \newline
!    gridDesc(32) = TBD \newline
!    gridDesc(33) = TBD \newline
!    gridDesc(34) = TBD \newline
!    gridDesc(35) = TBD \newline
!    gridDesc(36) = TBD \newline
!    gridDesc(37) = TBD \newline
!    gridDesc(38) = TBD \newline
!    gridDesc(39) = TBD \newline
!    gridDesc(40) = TBD \newline
! 
!   \item[Lambert conformal projection] 
!    gridDesc(1) = 3 ; indicates lambert projection \newline
!    gridDesc(2) = number of columns in the domain ; \newline
!    gridDesc(3) = number of rows in the domain ; \newline
!    gridDesc(4) = latitude of the lower left corner of the domain \newline
!    gridDesc(5) = longitude of the lower left corner of the domain \newline
!    gridDesc(6) = 8 ; not used \newline
!    gridDesc(7) =  true latitude 2 \newline
!    gridDesc(8) =  spatial resolution (in km) along E-W dimension \newline
!    gridDesc(9) =  spatial reolution (in km) along N-S dimension \newline
!    gridDesc(10) = true latitude 1 \newline
!    gridDesc(11) = standard longitude  \newline
!    gridDesc(20) = 255 ; used to specify the ordering of data (non-divisible by 
!                  32 indicates E-W ordering else N-S ordering \newline
!    gridDesc(30) = TBD \newline
!    gridDesc(32) = TBD \newline
!    gridDesc(33) = TBD \newline
!    gridDesc(34) = TBD \newline
!    gridDesc(35) = TBD \newline
!    gridDesc(36) = TBD \newline
!    gridDesc(37) = TBD \newline
!    gridDesc(38) = TBD \newline
!    gridDesc(39) = TBD \newline
!    gridDesc(40) = TBD \newline
!
!   \item[Gaussian projection] 
!    gridDesc(1) = 4 ; indicates gaussian projection \newline
!    gridDesc(2) = number of columns in the domain ; \newline
!    gridDesc(3) = number of rows in the domain ; \newline
!    gridDesc(4) = latitude of the lower left corner of the domain \newline
!    gridDesc(5) = longitude of the lower left corner of the domain \newline
!    gridDesc(6) = 8 ; not used \newline
!    gridDesc(7) = latitude of the upper right corner of the domain  \newline
!    gridDesc(8) =  longitude of the upper right corner of the domain \newline
!    gridDesc(9) =  spatial resolution (in degrees) along the E-W dimension \newline
!    gridDesc(10) = number of latitude circles \newline
!    gridDesc(11) = 64 ; not used  \newline
!    gridDesc(20) = 255 ; used to specify the ordering of data (non-divisible by 
!                  32 indicates E-W ordering else N-S ordering \newline
!    gridDesc(41) = 4 \newline
!    gridDesc(42) = number of columns in the domain  \newline
!    gridDesc(43) = number of rows in the domain \newline
!    gridDesc(44) = latitude of the lower left corner point  \newline
!    gridDesc(45) = longitude of the lower left corner point  \newline
!    gridDesc(46) = 128; not used \newline
!    gridDesc(47) = latitude of the upper right corner point  \newline
!    gridDesc(48) = latitude of the upper right corner point  \newline
!    gridDesc(49) = spatial resolution in the E-W dimension \newline
!    gridDesc(50) = number of latitude circles \newline
!   \item[polar stereographic projection] 
!    gridDesc(1) = 5 ; indicates polar stereographic projection \newline
!    gridDesc(2) = number of columns in the domain ; \newline
!    gridDesc(3) = number of rows in the domain ; \newline
!    gridDesc(4) = latitude of the lower left corner of the domain \newline
!    gridDesc(5) = longitude of the lower left corner of the domain \newline
!    gridDesc(6) = 8 ; not used \newline
!    gridDesc(7) =  orientation of the grid \newline
!    gridDesc(8) =  spatial resolution (in km) along E-W dimension \newline
!    gridDesc(9) =  spatial reolution (in km) along N-S dimension \newline
!    gridDesc(10) = true latitude \newline
!    gridDesc(11) = standard longitude  \newline
!    gridDesc(20) = 255 ; used to specify the ordering of data (non-divisible by 
!                  32 indicates E-W ordering else N-S ordering \newline
!    gridDesc(30) = TBD \newline
!    gridDesc(32) = TBD \newline
!    gridDesc(33) = TBD \newline
!    gridDesc(34) = TBD \newline
!    gridDesc(35) = TBD \newline
!    gridDesc(36) = TBD \newline
!    gridDesc(37) = TBD \newline
!    gridDesc(38) = TBD \newline
!    gridDesc(39) = TBD \newline
!    gridDesc(40) = TBD \newline
!   \item[UTM projection] 
!    gridDesc(1) = 7 ; indicates UTM projection \newline
!    gridDesc(2) = number of columns in the domain ; \newline
!    gridDesc(3) = number of rows in the domain ; \newline
!    gridDesc(4) = northing of the lower left corner of the domain \newline
!    gridDesc(5) = easting of the lower left corner of the domain \newline
!    gridDesc(6) = 128 ; not used \newline
!    gridDesc(7) = northing of the upper right corner of the domain \newline
!    gridDesc(8) = easting of the upper right corner of the domain \newline
!    gridDesc(9) = spatial resolution (in meters) \newline
!    gridDesc(10) = UTM zone \newline
!    gridDesc(11) = 64 ; not used \newline
!    gridDesc(20) = 255 ; used to specify the ordering of data (non-divisible by 
!                  32 indicates E-W ordering else N-S ordering \newline
!    gridDesc(34) = UTM zone \newline
!    gridDesc(35) = northing of the lower left corner point \newline
!    gridDesc(37) = easting of the lower left corner point \newline
!    gridDesc(38) = number of columns in the domain \newline
!    gridDesc(39) = number of rows in the domain \newline
!    gridDesc(40) = spatial resolution (in meters) \newline
!   \end{description}
!  \item[soil\_gridDesc, topo\_gridDesc, lc\_gridDesc]
!  Array describing the soils, topography, and landmask/landcover
!  data grid, respectively, for each nest. The specification 
!  of these arrays for different map projections are 
!  as follows: 
!   \begin{description}
!    \item[Lat/lon or geographic projection] 
!     gridDesc(1) = latitude of the lower left corner of the domain \newline
!     gridDesc(2) = longitude of the lower left corner of the domain \newline
!     gridDesc(3) =  latitude of the upper right corner of the domain \newline
!     gridDesc(4) =  longitude of the upper right corner of the domain \newline
!     gridDesc(5) =  spatial resolution (in degrees) along the E-W dimension \newline
!     gridDesc(6) = spatial resolution (in degrees) along the N-S dimension  \newline
!   \item[Gaussian projection] 
!     gridDesc(1) = latitude of the lower left corner of the domain \newline
!     gridDesc(2) = longitude of the lower left corner of the domain \newline
!     gridDesc(3) = latitude of the upper right corner of the domain \newline
!     gridDesc(4) = longitude of the upper right corner of the domain \newline
!     gridDesc(5) = spatial resolution (in degrees) along the E-W dimension \newline
!     gridDesc(6) = number of latitude circles  \newline
!   \item[UTM projection] 
!     gridDesc(1) = UTM zone
!     gridDesc(2) = northing of the lower left corner of the domain \newline
!     gridDesc(3) = easting of the lower left corner of the domain \newline
!     gridDesc(4) = number of columns in the domain \newline
!     gridDesc(5) = number of rows in the domain \newline
!     gridDesc(6) = spatial resolution (in meters) \newline
!   \end{description}
!  \item[nf]
!   Number of forcing variables
!  \item[nmetforc]
!   Number of meteorological forcing datasets used in a LIS simulation.
!  \item[nperforc]
!   Number of ensemble members that correspond to a single forcing source. 
!  \item[metforc\_blend\_alg]
!   Blending algorithm used for merging different forcing datasets. 
!  \item[metforc]
!   Choice of meteorological forcings
!  \item[met\_ecor]
!   Choice of topographical downscaling method for met forcing
!  \item[met\_nf]
!   Number of forcing variables in each met forcing dataset. 
!  \item[met\_interp]
!   Spatial interpolation option for met forcing (1-bilinear, 2-conservative, 3-neighbor)
!  \item[met\_upscale]
!   Spatial upscaling option for met forcing, currently only supports
!   upscaling by averaging
!  \item[met\_tinterp]
!   Temporal interpolation option for met forcing (1-next, 2-uber-next)
!  \item[rstflag]
!   boolean array denoting if the current simulation is a restart or not,
!   for each nest
!  \item[gridchange]
!   boolean array denoting if the native grid for the meteorological forcing
!   need to be changed. 
!  \item[shortflag]
!   Shortwave radiation source flag (1-instantaneous, 2-time averaged)
!  \item[longflag]
!   Longwave radiation source flag (1-instantaneous, 2-time averaged)
!  \item[perturb\_forcing]
!  Choice of forcing perturbation algorithm (0-do not use)
!  \item[pertforcInterval]
!  Forcing perturbation frequency (in seconds)
!  \item[perturb\_obs]
!  Choice of observation perturbation algorithm (0-do not use)
!  \item[pertobsInterval]
!  Observation perturbation frequency (in seconds)
!  \item[perturb\_state]
!  Choice of state perturbation algorithm (0-do not use)
!  \item[pertstateInterval]
!  State perturbation frequency (in seconds)
!  \item[pertrestart]
!  ??
!  \item[pertrestartInterval]
!  ??
!  \item[pertrestartfile]
!  ??
!  \item[nt]
!  Number of vegetation classes in the landcover dataset
!  \item[bareclass]
!  Index of bare class in the landcover dataset
!  \item[urbanclass]
!  Index of urban class in the landcover dataset
!  \item[snowclass]
!  Index of snow class in the landcover dataset
!  \item[waterclass]
!  Index of water class in the landcover dataset
!  \item[laiflag]
!  flag to control the LAI data read
!  \item[saiflag]
!  flag to contro the SAI data read
!  \item[mfile]
!  Name of the landmask file, for each nest
!  \item[vfile]
!  Name of the landcover file, for each nest
!  \item[vfile\_form]
!  
!  \item[safile]
!  Name of the sand fraction data file, for each nest
!  \item[clfile]
!  Name of the clay fraction data file, for each nest
!  \item[sifile]
!  Name of the silt fraction data file, for each nest
!  \item[txtfile]
!  Name of the soil texture data file, for each nest
!  \item[pofile]
!  Name of the porosity data file, for each nest
!  \item[psisatfile]
!  Name of the saturated matric potential data file, for each nest
!  \item[ksatfile]
!  Name of the saturated hydraulic conductivity data file, for each nest
!  \item[bexpfile]
!  Name of the b parameter data file, for each nest
!  \item[qzfile]
!  Name of the quartz data file, for each nest
!  \item[iscfile]
!  Name of the soil color data file, for each nest
!  \item[elevfile]
!  Name of the elevation data file, for each nest
!  \item[slfile]
!  Name of the slope data file, for each nest
!  \item[aspfile]
!  Name of the aspect data file, for each nest
!  \item[curvfile]
!  Name of the curvature data file, for each nest
!  \item[albfile]
!  Name of the climatology albedo data file, for each nest
!  \item[mxsnal]
!  Name of the max snow albedo data file, for each nest
!  \item[tbotfile]
!  Name of the bottom temperature data file, for each nest
!  \item[shdmaxfile]
!  Name of the maximum greenness data file, for each nest
!  \item[shdminfile]
!  Name of the minimum greenness data file, for each nest
!  \item[slopetypefile]
!  Name of the slope type data file, for each nest
!  \item[tile\_coord\_file]
!  Name of the catchment-based tile coordinate file, for each nest
!  \item[tile\_veg\_file]
!  Name of the catchment-based vegetation type file, for each nest
!  \item[outputSpecFile]
!  Name of the model output specification file
!  \item[output\_at\_specifictime]
!  Option to specify if output is to be written only at a particular
!  time (instead of at synoptic intervals)
!  \item[albInterval]
!  Frequency of albedo climatology in months
!  \item[laiInterval]
!  Frequency of LAI climatology in months
!  \item[laitime]
!  Variable to keep track of LAI data interval
!  \item[saitime]
!  Variable to keep track of SAI data interval
!  \item[wopt]
!  Output methodology (none, 1d tilespace,2d gridspace, 1d gridspace)
!  \item[wout]
!  Output data format (binary, GRIB-1, netcdf)
!  \item[wsingle]
!  Option to write each variable to a separate file 
!  \item[wstyle]
!  Output file naming style (3 level hierarchy, 5 level, WMO convention)
!  \item[sout]
!  Option to specify whether to write the statistics file
!  \item[grib\_table]
!  GRIB Table version number
!  \item[grib\_center\_id]
!  GRIB center id
!  \item[grib\_subcenter\_id]
!  GRIB subcenter id
!  \item[grib\_process\_id]
!  GRIB process id
!  \item[grib\_grid\_id]
!  GRIB grid id
!  \item[startcode]
!  Start mode (restart, coldstart)
!  \item[plevel]
!  Logging level
!  \item[odir]
!  Output directory
!  \item[dfile]
!  Diagnostic output file 
!  \item[sdoy]
!  Starting julian day
!  \item[sss]
!  Starting second 
!  \item[smn]
!  Starting minute
!  \item[shr]
!  Starting hour
!  \item[sda]
!  Starting day
!  \item[smo]
!  Starting month
!  \item[syr]
!  Starting year
!  \item[endcode]
!  End mode (1-specific date)
!  \item[edoy]
!  Ending julian day
!  \item[ess]
!  Ending second 
!  \item[emn]
!  Ending minute
!  \item[ehr]
!  Ending hour
!  \item[eda]
!  Ending day
!  \item[emo]
!  Ending month
!  \item[eyr]
!  Ending year
!  \item[endtime]
!  Flag to indicate the end of simulation (1-end of simulation)
!  \item[etime]
!  Ending time
!  \item[egmt]
!  Ending time in GMT
!  \item[nts]
!  Array containing the model timestep for each nest
!  \item[ts]
!  Timestep for the clock (minimum timestep of different nests)
!  \item[tscount]
!  Timestep count for each nest
!  \item[doy]
!  Current julian day
!  \item[ss]
!  Current second 
!  \item[mn]
!  Current minute
!  \item[hr]
!  Current hour
!  \item[da]
!  Current day
!  \item[mo]
!  Current month
!  \item[yr]
!  Current year
!  \item[time]
!  Current time
!  \item[gmt]
!  Current time in GMT
!  \item[ndas]
!  Number of data assimilation instances
!  \item[nperts]
!  Number of perturbation instances
!  \item[daalg]
!  Choice of data assimilation algorithm
!  \item[biasalg]
!  Choice of bias estimation algorithm
!  \item[incroption]
!  Option to specify whether to apply both bias and analysis increments
!  (1- apply both, 0 - apply bias correction only)
!  \item[daset]
!  Choice of Assimilation set: (variable to be assimilated + observation
!  dataset)
!  \item[nstvars]
!  Number of state variables to be assimilated
!  \item[nobtypes]
!  Number of observation types to be assimilated
!  \item[daoutInterval]
!  output interval for DA diagnostics
!  \item[nensem]
!  Number of ensembles per grid
!  \item[forcvarlistFile]
!  name of forcing variable list file
!  \item[forcattribFile]
!  name of forcing attributes file
!  \item[forcpertattribFile]
!  name of forcing perturbation attributes file
!  \item[progpertattribFile]
!  name of prognostic variable attributes file
!  \item[progpertattribFile]
!  name of prognostic variable perturbation attributes file
!  \item[obsattribFile]
!  name of observation attributes file
!  \item[obspertattribFile]
!  name of observation perturbation attributes file
!  \item[biasOptionsFile]
!  name of bias specification attributes file
!  \item[nforcepert]
!  Number of perturbation
!  \item[wensems]
!  flag to specify if ensemble members are to be output
!  \item[wobs]
!  flag to specify if observations are to be output
!  \item[winnov]
!  flag to specify if normalized innovations are to be output
!  \item[optUEAlg]
!  Choice of optimization/uncertainty estimation algorithm
!  \item[optUEset]
!   Dataset to calibrate to
!  \item[optUEtype]
!  Choice of optimization/uncertainty estimation
!  type (parameter estimation, tuning DA params, etc.)
!  \item[julbeg]
!  Beginning julian time of the current LIS cycle
!  \item[julend]
!  Ending julian time of the current LIS cycle
! \item[security\_class]
!   Security classification code for the LIS application
! \item[distribution\_class]
!   Distribution classification code for the LIS application
! \item[data\_category]
!   Data category code for the LIS application
! \item[area\_of\_data]
!  Geographical code for the LIS domain
! \item[lis\_config\_file]
!  Name of the run-time LIS configuration file (default = lis.config)
!  \item[rtm]
!   Choice of radiative transfer model
!  \item[nappmodel]
!   Total number of application models
!  \item[landslidemodel]
!   Choice of landslide model. 
!  \end{description}
! !REVISION HISTORY:
!  12 Apr 2001: Urszula Jambor; Added domain,lsm,& force namefile paramters 
!  30 Jul 2001: Matt Rodell; Add new soil parameter variables
!  14 Nov 2002; Sujay Kumar; Optimized version for LIS
!  14 Oct 2003; Sujay Kumar; Removed LSM specific variables. 
!  19 Jan 2007; Chuck Alonge; Added Flag to output parameters
!  17 Jan 2011: David Mocko, added max/min greenness & slope type
!
!EOP
  implicit none
  type lisrcdec
     character*50           :: runmode 
     integer                :: nnest 
    
     integer                :: max_model_types 
     integer                :: nsf_model_types
     integer, allocatable       :: sf_model_type(:)
     character*50, allocatable  :: sf_model_type_name(:)
     integer, allocatable       :: sf_model_type_select(:)
     character*50, allocatable  :: sf_model_type_name_select(:)

     integer                :: lsm_index
     integer                :: lake_index 
     integer                :: glacier_index 
     integer                :: wetland_index 
     integer                :: openwater_index 

     integer, allocatable       :: ntiles(:)
     integer, allocatable       :: glbntiles(:)
     integer, allocatable       :: glbntiles_red(:)
!tilespace in each surface model, which is essentially a "patch"
     integer, allocatable       :: npatch(:,:)
     integer, allocatable       :: glbnpatch(:,:)
     integer, allocatable       :: glbnpatch_red(:,:)
     
     integer, allocatable       :: ngrid(:) 
     integer, allocatable       :: obs_ngrid(:) 
     integer, allocatable       :: glbngrid(:)
     integer, allocatable       :: obs_glbngrid(:)
     integer, allocatable       :: obs_glbngrid_red(:)
     integer, allocatable       :: glbngrid_red(:)
     integer, allocatable       :: gnc(:)
     integer, allocatable       :: gnr(:)
     integer, allocatable       :: gnc_b(:)
     integer, allocatable       :: gnr_b(:)
     integer, allocatable       :: lnc(:)
     integer, allocatable       :: lnr(:)
     integer, allocatable       :: lnc_b(:)
     integer, allocatable       :: lnr_b(:)
     integer, allocatable       :: lnc_red(:)
     integer, allocatable       :: lnr_red(:)
     integer, allocatable       :: ncatg(:)

     integer, allocatable       :: obs_lnc(:)
     integer, allocatable       :: obs_lnr(:)
     integer, allocatable       :: obs_lnc_red(:)
     integer, allocatable       :: obs_lnr_red(:)
     integer, allocatable       :: obs_gnc(:)
     integer, allocatable       :: obs_gnr(:)
     integer, allocatable       :: obs_halox(:)
     integer, allocatable       :: obs_haloy(:)

     character*50               :: lis_map_proj
     character*50, allocatable  :: lis_obs_map_proj(:)
     character*50               :: lsm
     character*50               :: lakemodel
     character*50               :: glaciermodel
     character*50               :: openwatermodel
     character*50               :: param_proj
    
     character*100, allocatable :: paramfile(:)
     character*100, allocatable :: obsdomainfile(:)
     character*50, allocatable  :: usemaskmap(:)
     character*50, allocatable  :: uselcmap(:)
     character*50           :: lcscheme
     character*50           :: cropscheme
     integer                :: numbercrops
     character*50, allocatable  :: usetexturemap(:)
     character*50, allocatable  :: usesoilfractionmap(:)
     character*50, allocatable  :: usesoilcolormap(:)
     character*50, allocatable  :: useelevationmap(:)
     character*50, allocatable  :: useslopemap(:)
     character*50, allocatable  :: useaspectmap(:)
     character*50, allocatable  :: usecurvaturemap(:)
     character*50, allocatable  :: uselaimap(:)
     character*50, allocatable  :: usesaimap(:)
     character*50, allocatable  :: usealbedomap(:)
     character*50, allocatable  :: usemxsnalbmap(:)
     character*50, allocatable  :: usegreennessmap(:)
     character*50, allocatable  :: useemissmap(:)
     character*50, allocatable  :: useroughnessmap(:)
     character*50, allocatable  :: useporositymap(:)
     character*50, allocatable  :: usepsisatmap(:)
     character*50, allocatable  :: useksatmap(:)
     character*50, allocatable  :: usebexpmap(:)
     character*50, allocatable  :: usequartzmap(:)
     integer, allocatable       :: usesnowmap(:)
     
     integer, allocatable       :: soilsrc(:)
     integer, allocatable       :: snowsrc(:)
     integer                :: surface_maxt
     real                   :: surface_minp    
     integer                :: soilt_maxt
     real                   :: soilt_minp    
     integer                :: soilf_maxt
     real                   :: soilf_minp    
     integer                :: elev_maxt
     real                   :: elev_minp    
     integer                :: slope_maxt
     real                   :: slope_minp    
     integer                :: aspect_maxt
     real                   :: aspect_minp    

     logical                :: decompose_by_processes
     integer                :: npesx
     integer                :: npesy
     integer                :: halox
     integer                :: haloy
     real                   :: udef              
     real, allocatable      :: gridDesc(:,:)     
     real, allocatable      :: obs_gridDesc(:,:)     
     real, allocatable      :: minLat(:),maxLat(:)
     real, allocatable      :: minLon(:),maxLon(:)
!     real, allocatable          :: soil_gridDesc(:,:)
!     real, allocatable          :: topo_gridDesc(:,:)
!     real, allocatable          :: lc_gridDesc(:,:)  

     integer                :: nf   
     integer                :: nmetforc
     character*50           :: metforc_blend_alg
     character*50, allocatable  :: metforc(:)
     character*50, allocatable  :: met_ecor(:)       
     integer,      allocatable  :: metforc_ensmem(:)
     integer, allocatable       :: pcp_downscale(:)       
     integer, allocatable       :: met_nf(:)
     integer, allocatable       :: met_nensem(:)
     integer, allocatable       :: met_nperforc(:)
     
     character*50, allocatable  :: met_interp(:)
     character*50, allocatable  :: met_upscale(:)
     character*50, allocatable  :: met_tinterp(:)
     character*50, allocatable  :: met_proj(:)

     integer, allocatable       :: rstflag(:) 
     integer, allocatable       :: gridchange(:)
     integer                :: shortflag
     integer                :: longflag 

     character*50           :: perturb_forcing    
     real                   :: pertforcInterval
     character*50, allocatable  :: perturb_obs(:)
     real, allocatable          :: pertobsInterval(:)
     character*50, allocatable  :: perturb_state(:)
     real, allocatable          :: pertstateInterval(:)
     character*50           :: pertrestart
     real                   :: pertrestartInterval
     character*100, allocatable :: pertrestartfile(:)
     integer                :: pert_bias_corr

     integer                :: nsurfacetypes
     integer                :: nvegtypes
     integer                :: nsoiltypes
     integer                :: nsoilfbands
     integer                :: nelevbands
     integer                :: nslopebands
     integer                :: naspectbands

     integer                :: bareclass 
     integer                :: urbanclass
     integer                :: snowclass 
     integer                :: waterclass
     integer                :: wetlandclass
     integer                :: glacierclass
     integer                :: cropclass
     integer                :: laiflag  
     integer                :: saiflag       
     character*100, allocatable :: mfile(:)  
     character*100, allocatable :: vfile(:)
     integer,       allocatable :: vfile_form(:)
     character*100, allocatable :: safile(:) 
     character*100, allocatable :: clfile(:) 
     character*100, allocatable :: sifile(:) 
     character*100, allocatable :: txtfile(:)
     character*100, allocatable :: pofile(:)
     character*100, allocatable :: psisatfile(:) 
     character*100, allocatable :: ksatfile(:) 
     character*100, allocatable :: bexpfile(:) 
     character*100, allocatable :: qzfile(:)   
     character*100, allocatable :: iscfile(:) 
     character*100, allocatable :: elevfile(:)
     character*100, allocatable :: slfile(:)
     character*100, allocatable :: aspfile(:)
     character*100, allocatable :: curvfile(:)
     character*100, allocatable :: albfile(:)  
     character*100, allocatable :: mxsnal(:)   
     character*100, allocatable :: tbotfile(:) 
     integer                :: tbot_terrain_adj
     integer                :: tbot_update_lag
     integer                :: tbot_lagday
     character*100, allocatable :: shdmaxfile(:) 
     character*100, allocatable :: shdminfile(:) 
     character*100, allocatable :: slopetypefile(:) 
     character*100, allocatable :: tile_coord_file(:) 
     character*100, allocatable :: tile_veg_file(:)
     character*100, allocatable :: outputSpecFile(:)
     integer                :: output_at_specifictime
     integer                :: albInterval
     integer                :: laiInterval
     real*8                 :: laitime  
     real*8                 :: saitime  
     character*50           :: wopt     
     integer                :: wopt_rst
     character*50           :: wout           
     integer                :: wsingle
     character*50           :: wstyle
     logical                :: sout           
     integer                :: grib_table
     integer                :: grib_center_id
     integer                :: grib_subcenter_id
     integer                :: grib_grid_id
     integer                :: grib_process_id
     character*50           :: grib_packing_type
     character*50           :: startcode
     integer                :: plevel
     character*100          :: odir
     character*100          :: dfile      
     integer                :: sdoy        
     integer                :: sss         
     integer                :: smn         
     integer                :: shr         
     integer                :: sda         
     integer                :: smo         
     integer                :: syr         
     integer                :: endcode        
     integer                :: ess            
     integer                :: emn            
     integer                :: edoy           
     integer                :: ehr            
     integer                :: eda            
     integer                :: emo            
     integer                :: eyr            
     integer                :: endtime        
     real*8                 :: etime               
     real                   :: egmt
     real                   :: twInterval

     integer                :: monthCount
     integer                :: alrm_prev_mo

     integer                :: ess1            
     integer                :: emn1            
     integer                :: edoy1           
     integer                :: ehr1            
     integer                :: eda1            
     integer                :: emo1            
     integer                :: eyr1
     real*8                 :: etime1
     real                   :: egmt1
     real,    allocatable   :: nts(:)     
     real                   :: ts         
     integer, allocatable   :: tscount(:)
     integer                :: doy
     integer                :: yr
     integer                :: mo
     integer                :: da
     integer                :: hr
     integer                :: mn
     integer                :: ss 
     integer                :: ms
     real*8                 :: time      
     real                   :: gmt
     integer                :: ndas
     integer                :: nperts
     character*50, allocatable  :: daalg(:)
     
     integer, allocatable       :: useANNinDA(:)
     character*100, allocatable :: ANNdaFile(:)

     character*50, allocatable  :: biasalg(:)
     character*50, allocatable  :: biasrst(:)
     real, allocatable          :: biasrstInterval(:)
     integer, allocatable       :: incroption(:)
     character*50, allocatable  :: daset(:)
     integer, allocatable       :: nstvars(:)
     integer, allocatable       :: nobtypes(:)
     real, allocatable          :: daoutInterval(:)
     integer, allocatable       :: nensem(:)
     character*100, allocatable :: biasrstfile(:)
     character*100          :: forcvarlistFile
     character*100          :: forcattribFile
     character*100          :: forcpertattribFile
     character*100, allocatable :: progpertattribFile(:)
     character*100, allocatable :: progattribFile(:)
     character*100, allocatable :: obspertattribFile(:)
     character*100, allocatable :: obsattribFile(:)
     character*100, allocatable :: biasOptionsFile(:)
     character*50, allocatable  :: dascaloption(:)

     integer                :: nforcepert
     integer, allocatable   :: wensems(:)
     integer, allocatable   :: wobs(:)
     integer, allocatable   :: winnov(:)
     integer, allocatable   :: DAincrMode(:)
     integer, allocatable   :: iterationId(:)

     character*50           :: optUEAlg
     character*50           :: optUEset
     character*50           :: optUEtype
     character*50           :: objfuncmethod
     integer                :: decSpaceInitMode
     integer                :: wpeobs

     integer                :: julbeg
     integer                :: julend
     logical                :: use_twelve
     logical                :: reset_flag
     logical                :: run_model

     character*20           :: security_class
     character*20           :: distribution_class
     character*20           :: data_category
     character*20           :: area_of_data
     character*100          :: lis_config_file='lis.config'
     character*100          :: institution = 'NASA GSFC'
!RTM related variables
     character*50           :: rtm
!landslide model related variables
     integer                :: nappmodel
     character*50           :: landslidemodel ! SY   

     character*50           :: routingmodel
     character*50           :: runoffdatasource

     character*50           :: irrigation_type
     real                   :: irrigation_thresh !BZ
     integer                :: irrigation_mxsoildpth

     integer                :: forecastMode
     logical                :: zterp_correction

     real                   :: irrigation_GVFparam1   !WN
     real                   :: irrigation_GVFparam2   !WN
     integer                :: irrigation_GWabstraction !JE 
 
  end type lisrcdec
  
end module LIS_PRIV_rcMod
