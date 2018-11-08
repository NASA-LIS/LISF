subroutine ruc_driver_37(                                                  &
  year, month, day, hour, minute, ktime,                                   &
  lwdown, swdown, psurf, rainf, snowf, tair, qair, wind_e, wind_n,         &
  ivgtyp, isltyp, dt, sldpth2, zlvl, zlvl_wind,                            &
  local_param, rdlai2d, usemonalb, iz0tlnd, sfcdif_option,                 &
  landuse_tbl_name, soil_tbl_name, gen_tbl_name,                           &
  landuse_scheme_name, soil_scheme_name,                                   &
  nsoil, water_class_num, ice_class_num, urban_class_num,                  &
  albedo_monthly, shdfac_monthly, z0brd_monthly, lai_monthly,              &
  albbck1, tbot1, snoalb1, emiss1, ch_1, cm_1, sneqv1, snowh1,             &
  snowc1, canwat1, alb1, smc2, sh2o2, stc2, smfr3, keepfr3,                &
  tskin1, qvg1, qsfc1, qcg1, qsg1, snt75cm1, tsnav1, soilm1, smroot1,      &
  shdfac1, shdmin1, shdmax1, qh1, qle1, eta1, qg1, rhosnf1, precipfr1,     &
  snowfallac1, qmax2,qmin2,psis1,ksat1,bclh1,qwrtz1,wilt1,ref1,            &
  snthresh1, acsnow1, sfcevp1, snomlt1, dew1, drip1, qs1,                  &
  qsb1, snflx1, edir1, ec1, ett1, esnow1)                                   
  !
  ! module module_ruclsm_utility contains a few utility routines needed for setting
  ! up some data for ruc lsm.  

  !this driver program makes use of module subroutines
  ! caltmp and calhum.
  !

  use module_Noahlsm_utility_ruc37
  !
  ! module_sf_ruclsm contains the ruc lsm physics code.  this driver program
  ! makes use of module subroutines lsmruc, soilvegin, and sfcdif_off.
  !

  use module_sf_ruclsm_ruc37

  !
  ! module_sfcdif_wrf contains the myjsfc version of sfcdif.  this driver program
  ! makes use of subroutines myjsfcinit and sfcdif_myj.
  !

   use module_sfcdif_wrf_ruc37

  !
  ! kwm_date_utilities contains handy subroutines for manipulating and
  ! computing date/time information.
  !

  !use kwm_date_utilities

  implicit none
  integer, intent(in) :: year 
  integer, intent(in) :: month 
  integer, intent(in) :: day 
  integer, intent(in) :: hour 
  integer, intent(in) :: minute
  integer, intent(in) :: ktime   ! a counter for the timesteps in the main loop timeloop.

  real, intent(in)    :: lwdown           ! downward longwave radiation flux at surface (w m-2) [forcing]
  real, intent(in)    :: swdown           ! downward shortwave radiation flux at surface (w m-2) [forcing]
  real, intent(in)    :: psurf            ! surface atmospheric pressure (pa) [forcing]
  real, intent(in)    :: rainf            ! rainfall rate (kg m-2 s-1) [forcing]
  real, intent(in)    :: snowf            ! snowfall rate (kg m-2 s-1) [forcing]
  real, intent(in)    :: tair             ! air temperature (k) [forcing]
  real, intent(in)    :: qair             ! surface specific humidity (kg kg-1) [forcing]
  real, intent(in)    :: wind_e           ! eastward wind speed (m s-1) [forcing]
  real, intent(in)    :: wind_n           ! northward wind speed (m s-1) [forcing]
  
  integer, intent(in) :: ivgtyp           ! vegetation category
  integer, intent(in) :: isltyp           ! soil category 
  real, intent(in)    :: dt               ! time step (seconds).
  real, intent(in)    :: sldpth2(nsoil)   ! thicknesses of each soil level (m) 
  logical, intent(in) :: local_param      ! .true. to use table values for albbck, shdfac, and z0brd; .false. to use values for albbck, shdfac, and z0brd as set in this driver routine
  logical, intent(in) :: rdlai2d          ! if rdlai2d == .true., then the xlai value that we pass to lsmruc will be used. if rdlai2d == .false., then xlai will be computed within lsmruc, from table minimum and maximum values in vegparm.tbl, and the current green vegetation fraction.
  logical, intent(in) :: usemonalb        ! if usemonalb == .true., then the alb value passed to lsmruc will be used as the background snow-free albedo term.  if usemonalb == .false., then alb will be computed within lsmruc from minimum and maximum values in vegparm.tbl, and the current green vegetation fraction.
  integer,intent(in)  :: iz0tlnd          ! option to turn on (iz0tlnd=1) or off (iz0tlnd=0) the vegetation-category-dependent calculation of the zilitinkivich coefficient czil in the sfcdif subroutines.
  integer, intent(in) :: sfcdif_option    ! option to use previous (sfcdif_option=0) or updated (sfcdif_option=1) version of sfcdif subroutine.
  
  character(len=256), intent(in) :: landuse_tbl_name      ! noah model landuse parameter table
  character(len=256), intent(in) :: soil_tbl_name         ! noah model soil parameter table
  character(len=256), intent(in) :: gen_tbl_name          ! noah model general parameter table
  character(len=256), intent(in) :: landuse_scheme_name   ! landuse classficiation scheme
  character(len=256), intent(in) :: soil_scheme_name      ! soil classification scheme  
  integer, intent(in) :: nsoil            ! number of soil levels.
  
  integer, intent(in) :: water_class_num  ! number of water category in llanduse classification
  integer, intent(in) :: ice_class_num    ! number of ice category in llanduse classification
  integer, intent(in) :: urban_class_num    ! number of urban category in llanduse classification

  
  real, dimension(12) :: albedo_monthly   ! monthly values of background (i.e., snow-free) albedo ( fraction [0.0-1.0] )
  real, dimension(12) :: shdfac_monthly   ! monthly values for green vegetation fraction ( fraction [0.0-1.0] )
  real, dimension(12) :: z0brd_monthly    ! monthly values for background (i.e., snow-free) roughness length ( m )
  real, dimension(12) :: lai_monthly      ! monthly values for leaf area index ( dimensionless )
  
  
  real, intent(inout)    :: albbck1          ! background snow-free albedo (0.0-1.0).
  real, intent(in)    :: tbot1            ! deep-soil time-invariant temperature (k).  representing sort of a mean annual air temperature.
  real, intent(in)    :: snoalb1          ! maximum snow albedo over deep snow (0.0-1.0)

  real, intent(inout) :: emiss1           ! surface emissivity (0.0 - 1.0).
  real, intent(inout) :: ch_1              ! exchange coefficient for head and moisture (m s-1). 
  real, intent(inout) :: cm_1              ! exchange coefficient for momentum (m s-1).  
  real, intent(inout) :: sneqv1           ! water equivalent of accumulated snow depth (m).
  real, intent(inout) :: snowh1           ! physical snow depth (m). 
  real, intent(inout) :: snowc1           ! fractional snow cover ( fraction [0.0-1.0] )
  real, intent(inout) :: canwat1          ! canopy moisture content (kg m-2)
  real, intent(inout) :: alb1             ! surface albedo including possible snow-cover effect.  this is set in lsmruc,
  
  real, intent(inout) :: smc2(nsoil)      ! total soil moisture content (m3 m-3)
  real, intent(inout) :: sh2o2(nsoil)     ! liquid soil moisture content (m3 m-3) 
  real, intent(inout) :: stc2(nsoil)      ! soil temperature (k)    
  real, intent(inout) :: smfr3(nsoil)     ! soil ice content (m3 m-3)    
  real, intent(inout) :: keepfr3(nsoil)   ! frozen soil flag
  real, intent(inout) :: tskin1           ! skin temperature (k)
  real, intent(inout) :: qvg1             ! mixing ratio at the surface ( kg kg{-1} )
  real, intent(inout) :: qsfc1            ! specific humidity at the surface ( kg kg{-1} )
  real, intent(inout) :: qcg1             ! cloud water mixing ratio at the surface ( kg kg{-1} ) 
  real, intent(inout) :: qsg1             ! surface water vapor mixing ratio at satration (kg kg-1) 
  real, intent(inout) :: snt75cm1         ! snow temperature at 7.5 cm depth (k)
  real, intent(inout) :: tsnav1           ! average snow temperathre in k
  real, intent(inout) :: soilm1           ! total soil column moisture content, frozen and unfrozen ( m ) 
  real, intent(inout) :: smroot1          ! available soil moisture in the root zone ( fraction [smcwlt-smcmax] 
  
  real, intent(out)   :: qh1              ! sensible heat flux ( w m{-2} )
  real, intent(out)   :: qle1             ! latent heat flux (evapotranspiration) ( w m{-2} )
  real, intent(out)   :: eta1             ! latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
  real, intent(out)   :: qg1              ! soil heat flux ( w m{-2} )
  real, intent(out)   :: rhosnf1          ! density of frozen precipitation (kg m{-3}) 
  real, intent(out)   :: precipfr1        ! time-step frozen precipitation (kg m{-2}) 
  real, intent(inout)   :: snowfallac1      ! run total snowfall accumulation (kg m{-2})
  real, intent(inout)   :: acsnow1          ! run total frozen precipitation (kg m{-2})
  real, intent(inout)   :: sfcevp1          ! run total evaporation flux  (kg m{-2}) 
  real, intent(out)   :: snomlt1          ! snow melt water ( m )
  real, intent(out)   :: snthresh1        ! snow depth threshold ( m )
  real, intent(out)   :: dew1             ! dewfall (or frostfall for t<273.15) ( m ) 
        
  real, intent(out)   :: drip1            ! throughfall of precipitation from canopy (kg m{-2} s{-1}) 
  real, intent(out)   :: qs1              ! surface runoff, not infiltrating the soil ( m s{-1} )
  real, intent(out)   :: qsb1             ! subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
  real, intent(out)   :: snflx1           ! snow heat flux (w/m^2: negative, if downward from surface)
  real, intent(out)   :: edir1            ! latent heat flux component: direct soil evaporation ( w m{-2} )
  real, intent(out)   :: ec1              ! latent heat flux component: canopy water evaporation ( w m{-2} )
  real, intent(out)   :: ett1             ! latent heat flux component: total plant transpiration ( w m{-2} )
  real, intent(out)   :: esnow1           ! sublimation from snowpack (w m{-2})   
  real, intent(out)   :: shdfac1          ! vegetation fraction (0-1)
  real, intent(out)   :: shdmin1          ! min vegetation fraction (0-1)
  real, intent(out)   :: shdmax1          ! max vegetation fraction (0-1)
  real, intent(out)   :: qmax2 
  real, intent(out)   :: qmin2  
  real, intent(out)   :: psis1
  real, intent(out)   :: ksat1 
  real, intent(out)   :: bclh1  
  real, intent(out)   :: qwrtz1 
  real, intent(out)   :: wilt1
  real, intent(out)   :: ref1 

  
  character(len=32)   :: llanduse  ! land-use dataset.  valid values are :"usgs" (usgs 24/27 category dataset) and "modified_igbp_modis_noah" (modis 20-category dataset)
  character(len=32)   :: lsoil    ! soil-category dateset.  only "stas" (statsgo dataset) supported.
  integer, parameter  :: ime = 1
  integer, parameter  :: jme = 1
  integer, parameter  :: kme = 1


  character(len=12)  :: nowdate ! the date of each time step, ( yyyymmddhhmm )
  !
  ! various arguments to subroutine lsmruc:
  !

  integer :: iswater          ! number of water category in llanduse classification
  integer :: isice            ! number of ice category in llanduse classification
  logical :: frpcpn     ! flag: .true. - mixed phase precipitation
  real, dimension (1:ime,1:jme) :: ffrozp     ! fraction of precip which is frozen (0.0 - 1.0).
  integer, dimension (1:ime,1:jme) :: ice        ! flag for sea-ice (1) or land (0).
  real, dimension (1:ime,1:jme)    :: xice       ! flag for sea-ice (1) or land (0).
  integer :: fractional_seaice  ! 0 no fractional, 1 - yes
  real    :: xice_threshold  ! threshold below which no ice
  integer :: isurban    ! vegetation category for urban land class.
  real    :: zlvl       ! height at which atmospheric forcing variables are taken to be valid (m)
  real    :: zlvl_wind  ! height at which the wind forcing variable is taken to be valid (m)
  logical :: local      ! 
  logical :: myj        ! use myj scheme to compute exchange coefficients
  real    :: sfcu       ! west-to-east component of the surface wind (m s-1)
  real    :: sfcv       ! south-to-north component of the surface wind (m s-1)
  real, dimension (1:ime,1:jme)    :: lwdn       ! downward longwave radiation flux at surface (w m-2) [forcing]
  real, dimension (1:ime,1:jme)    :: soldn      ! downward shortwave radiation flux at surface (w m-2) [forcing]
  real, dimension (1:ime,1:jme)    :: solnet     ! net downward shortwave radiation flux at the surface (w m-2)
  real, dimension (1:ime,1:jme)    :: sfcprs     ! surface atmospheric pressure (pa) [forcing]
  real, dimension (1:ime,1:jme)    :: prcp       ! precipitation rate (kg m-2 s-1) [forcing]
  real, dimension (1:ime,1:jme)    :: tabs       ! air temperature (k) [forcing]
  real, dimension (1:ime,1:jme)    :: q2         ! surface specific humidity (kg kg-1) [forcing]
  real, dimension (1:ime,1:jme)    :: sfcspd     ! surface wind speed (m s-1) [forcing]
  real, dimension (1:ime,1:jme)    :: th2        ! potential temperature at level zlvl (k)
  real, dimension (1:ime,1:jme)    :: t1v       ! virtual skin temperature (k).  used in sfcdif_off for computing cm and ch, but not passed to lsmruc
  real, dimension (1:ime,1:jme)    :: th2v       ! virtual potential temperature at level zlvl (k).  used in sfcdif_off for computing cm and ch, but not passed to lsmruc
  real, dimension (1:ime,1:jme)    :: rho        ! air density (dummy value output from caltmp, needs to be passed to lsmruc).
  real, dimension (1:ime,1:jme)    :: q2sat      ! saturated specific humidity (kg kg-1)
  real, dimension (1:ime,1:jme)    :: dqsdt2     ! slope of the saturated specific humidity curve w.r.t. temperature.
  integer, dimension (1:ime,1:jme) :: vegtyp     ! vegetation category.
  integer, dimension (1:ime,1:jme) :: soiltyp    ! soil category.
  real, dimension (1:ime,1:jme)    :: shdfac     ! shade factor (0.0-1.0).
  real, dimension (1:ime,1:jme)    :: qmax
  real, dimension (1:ime,1:jme)    :: qmin
  real, dimension (1:ime,1:jme)    :: psis
  real, dimension (1:ime,1:jme)    :: ksat
  real, dimension (1:ime,1:jme)    :: bclh
  real, dimension (1:ime,1:jme)    :: qwrtz
  real, dimension (1:ime,1:jme)    :: wilt
  real, dimension (1:ime,1:jme)    :: ref
  real, dimension (1:ime,1:jme)    :: shdmin     ! minimum shade factor (0.0-1.0).
  real, dimension (1:ime,1:jme)    :: shdmax     ! maximum shade factor (0.0-1.0).
  real, dimension (1:ime,1:jme)    :: albbck     ! background snow-free albedo (0.0-1.0).
  real, dimension (1:ime,1:jme)    :: snoalb     ! maximum snow albedo over deep snow (0.0-1.0)
  real, dimension (1:ime,1:jme)    :: tbot       ! deep-soil time-invariant temperature (k).  representing sort of a mean annual air temperature.
  real, dimension (1:ime,1:jme)    :: z0brd      ! background z0 value (m).
  real, dimension (1:ime,1:jme)    :: z0         ! roughness length (m)
  real, dimension (1:ime,1:jme)    :: emissi     ! surface emissivity (0.0 - 1.0).  this includes the snow-cover effect.
  real, dimension (1:ime,1:jme)    :: cmc        ! canopy moisture content (kg m-2)
  real, dimension (1:ime,1:jme)    :: t1         ! skin temperature (k)
  real, allocatable, dimension(:) :: et   ! plant transpiration from each soil level.
  real, allocatable, dimension(:) :: smav ! soil moisture availability at each level, fraction between
                                          ! smcwlt (smav=0.0) and smcmax (smav=1.0)

  real, allocatable, dimension(:,:,:) :: tslb ! 3-d soil temperature (k)
  real, allocatable, dimension(:,:,:) :: soilmois ! 3-d soil moisture
  real, allocatable, dimension(:,:,:) :: sliqw ! 3-d liquid soil moisture
  real, allocatable, dimension(:,:,:) :: smfr3d  ! 3-d soil ice content (m3 m-3)
  real, allocatable, dimension(:,:,:) :: keepfr3dflag  ! 3-d frozen soil flag

  real, dimension (1:ime,1:jme)    :: rhosnf
  real, dimension (1:ime,1:jme)    :: snowfallac
  real, dimension (1:ime,1:jme)    :: graupelncv
  real, dimension (1:ime,1:jme)    :: snowncv
  real, dimension (1:ime,1:jme)    :: rainncv
  real, dimension (1:ime,1:jme)    :: precipfr
  real, dimension (1:ime,1:jme)    :: qcatm
  real, dimension (1:ime,1:jme)    :: qsg
  real, dimension (1:ime,1:jme)    :: soilt1
  real, dimension (1:ime,1:jme)    :: tsnav
  real, dimension (1:ime,1:jme)    :: xland
  real, dimension (1:ime,1:jme)    :: sfcevp
  real, dimension (1:ime,1:jme)    :: acsnow
  real, dimension (1:ime,1:jme)    :: mavail
  real, dimension (1:ime,1:jme)    :: fltot
  real, dimension (1:ime,1:jme)    :: snflx
  real, dimension (1:ime,1:jme)    :: snowh      ! physical snow depth.
  real, dimension (1:ime,1:jme)    :: snow       ! snow water equivalent in kg/m^2/s
  real, dimension (1:ime,1:jme)    :: sneqv      ! water equivalent of accumulated snow depth (m).
  real, dimension (1:ime,1:jme)    :: albedo     ! surface albedo including possible snow-cover effect.  this is set in lsmruc, overriding any value given; it should perhaps be intent(out) from lsmruc.
  real, dimension (1:ime,1:jme)    :: ch         ! exchange coefficient for head and moisture (m s-1).  an initial value is needed for sfcdif_off.
  real, dimension (1:ime,1:jme)    :: cm         ! exchange coefficient for momentum (m s-1).  an initial value is needed for sfcdif_off.
  real, dimension (1:ime,1:jme)    :: eta        ! latent heat flux (evapotranspiration) ( w m{-2} )
  real, dimension (1:ime,1:jme)    :: sheat      ! sensible heat flux ( w m{-2} )
  real, dimension (1:ime,1:jme)    :: etakin     ! latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
  real, dimension (1:ime,1:jme)    :: ec         ! latent heat flux component: canopy water evaporation ( w m{-2} )
  real, dimension (1:ime,1:jme)    :: edir       ! latent heat flux component: direct soil evaporation ( w m{-2} )
  real, dimension (1:ime,1:jme)    :: ett        ! latent heat flux component: total plant transpiration ( w m{-2} )
  real, dimension (1:ime,1:jme)    :: esnow      ! latent heat flux component: sublimation from (or deposition to) snowpack ( w m{-2} )
  real, dimension (1:ime,1:jme)    :: drip       ! precipitation or dew falling through canopy, in excess of canopy holding capacity ( m )
  real, dimension (1:ime,1:jme)    :: dew        ! dewfall (or frostfall for t<273.15) ( m )
  real, dimension (1:ime,1:jme)    :: ssoil      ! soil heat flux ( w m{-2} )
  real, dimension (1:ime,1:jme)    :: snomlt     ! snow melt water ( m )
  real, dimension (1:ime,1:jme)    :: snthresh   ! snow depth threshold ( m )
  real, dimension (1:ime,1:jme)    :: sncovr     ! fractional snow cover ( fraction [0.0-1.0] )
  real, dimension (1:ime,1:jme)    :: runoff1    ! surface runoff, not infiltrating the soil ( m s{-1} )
  real, dimension (1:ime,1:jme)    :: runoff2    ! subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
  real, dimension (1:ime,1:jme)    :: xlai       ! leaf area index ( dimensionless )
  real, dimension (1:ime,1:jme)    :: soilw      ! available soil moisture in the root zone ( fraction [smcwlt-smcmax] )
  real, dimension (1:ime,1:jme)    :: soilm      ! total soil column moisture content, frozen and unfrozen ( m )
  real, dimension (1:ime,1:jme)    :: q1         ! mixing ratio at the surface ( kg kg{-1} )
  real, dimension (1:ime,1:jme)    :: qsfc       ! specific humidity at the surface ( kg kg{-1} )
  real, dimension (1:ime,1:jme)    :: qcg        ! cloud water mixing ratio at the surface ( kg kg{-1} )
  real    :: ribb       ! bulk richardson number used to limit the dew/frost.
  integer :: nlcat            ! number of landuse categories
  integer :: nscat            ! number of soil categories
  real    :: ustar 

  !
  ! some diagnostics computed from the output of subroutine lsmruc
  !

  real, dimension (1:ime,1:jme) :: res       ! residual of the surface energy balance equation ( w m{-2} )

  !
  ! miscellaneous declarations
  !
  integer, parameter :: iunit = 10       ! fortran unit number for reading initial/forcing conditions file.
  real, external     :: ruc37_month_d          ! external function (follows this main program):  given an array (dimension 12)
  !                                      ! representing monthly values for some parameter, return a value for 
  !                                      ! a specified date.
  real                :: czil            ! zilitinkevich constant, read from genparm.tbl and used to compute surface
  !                                      ! exchange coefficients

  real, dimension (1:ime,1:jme) :: longwave  ! longwave radiation as read from the forcing data, which is immediately 
  !                                          ! adjusted (by the emissivity factor) to set variable lwdn.

  real, dimension (1:ime,1:kme,1:jme)    :: zlvl3d

  real, allocatable :: landusef(:)
  real, allocatable :: soilctop(:)

  integer    :: mosaic_lu
  integer    :: mosaic_soil
  integer :: k

  !
  ! allocate additonal arrays (dimensioned by the number of soil levels) which we will need for lsmruc.
  !
  allocate( soilmois     (1:ime, 1:nsoil, 1:jme) )
  allocate( sliqw        (1:ime, 1:nsoil, 1:jme) )
  allocate( smfr3d       (1:ime, 1:nsoil, 1:jme) )
  allocate( keepfr3dflag (1:ime, 1:nsoil, 1:jme) )
  allocate( tslb         (1:ime, 1:nsoil, 1:jme) )
  allocate( et ( nsoil ) )
  allocate( smav ( nsoil ) )
!  keepfr3dflag = 0.0 

  llanduse   = trim(landuse_scheme_name)  
  lsoil      = trim(soil_scheme_name) 
  iswater    = water_class_num
  isice      = ice_class_num 
  isurban    = urban_class_num 

  local = local_param 
  
  !!!! param 
  vegtyp(1,1)  = ivgtyp
  soiltyp(1,1) = isltyp 
  tbot(1,1)    = tbot1 
  snoalb(1,1)  = snoalb1
  
  !
  ! forcing session 
  ! 

  longwave(1,1)= lwdown                                ! downward longwave radiation flux at surface (w m-2) [forcing]
  soldn(1,1)   = swdown                                ! downward shortwave radiation flux at surface (w m-2) [forcing]
  sfcprs(1,1)  = psurf                                 ! surface atmospheric pressure (pa) [forcing]
  prcp(1,1)    = rainf + snowf                         ! precipitation rate (kg m-2 s-1) [forcing]
  tabs(1,1)    = tair                                  ! air temperature (k) [forcing]
  q2(1,1)      = qair                                  ! surface specific humidity (kg kg-1) [forcing]
  sfcspd(1,1)  = sqrt(wind_e*wind_e + wind_n*wind_n)   ! surface wind speed (m s-1) [forcing]
  sfcu    = wind_e                                ! west-to-east component of the surface wind (m s-1)
  sfcv    = wind_n                                ! south-to-north component of the surface wind (m s-1)
 
  !
  ! read our lookup tables and parameter tables:  vegparm.tbl, soilparm.tbl, genparm.tbl
  !

  ! call ruc37_soilvegprm( llanduse, lsoil, landuse_tbl_name, soil_tbl_name, gen_tbl_name)


  ribb       = 0.0         ! bulk richardson number used to limit the dew/frost.
  
  !
  ! czil is needed for the sfcdif_off step.  this comes from czil_data, as read
  ! from the genparm.tbl file, which is how redprm ultimately gets it as well:
  !
  czil = czil_data

  ! sanity check for snow  and snowh
  if(sneqv1 < 1.e-9 .or. snowh1 < 1.e-12) then
    sneqv1 = 0.0
    snowh1 = 0.0
  endif 

  !!!! state in 
  ch(1,1)     = ch_1
  cm(1,1)     = cm_1 
  sneqv(1,1)  = sneqv1 
  snowh(1,1)  = snowh1   ! [m]
  sncovr(1,1) = snowc1  
  snow(1,1)   = sneqv(1,1)*1.e3  ! [mm]
  albedo(1,1) = alb1 
  albbck(1,1) = albbck1
  soilt1(1,1) = snt75cm1        ! the physical meaning of soilt1 is confusing. it means snow temperature at 7.5 cm depth (k)
  t1(1,1)     = tskin1          ! skin temperature (k)
  q1(1,1)     = qvg1            ! mixing ratio at the surface ( kg kg{-1} )
  qsfc(1,1)   = qsfc1           ! specific humidity at the surface ( kg kg{-1} )
  qcg(1,1)    = qcg1            ! effective cloud water mixing ratio at the surface ( kg kg{-1} )
  qsg(1,1)    = qsg1            ! surface water vapor mixing ratio at satration (kg kg-1) 
  tsnav(1,1)  = tsnav1 - 273.15 ! from k to c 
  soilm(1,1)  = soilm1          ! total soil column moisture content, frozen and unfrozen ( m ) 
  soilw(1,1)  = smroot1         ! available soil moisture in the root zone ( fraction [smcwlt-smcmax]  
  cmc(1,1)    = canwat1 
  emissi(1,1) = emiss1 
  
  ! treat run total variables as state variables 
  snowfallac(1,1) =snowfallac1  
  acsnow(1,1)     =acsnow1      
  sfcevp(1,1)     =sfcevp1      
  do k=1,nsoil
      soilmois(1,k,1)     = smc2(k)
      sliqw   (1,k,1)     = sh2o2(k)
      smfr3d  (1,k,1)     = smfr3(k)
      keepfr3dflag(1,k,1) = keepfr3(k)
      tslb    (1,k,1)     = stc2(k)
  enddo

  if ( sfcdif_option .eq. 1 ) then
      call myjsfcinit()
  endif

  !!!  
  ! ktime = 1  
  write(nowdate, "(i0.4,i0.2,i0.2,i0.2,i0.2)"), year, month, day, hour, minute ! yyyymmddhhmm

  !
  ! update ffrozp for each time step, depending on the air temperature in the forcing data.
  ! ffrozp indicates the fraction of the total precipitation which is considered to be
  ! frozen.
  !

  if ( (prcp(1,1) > 0) .and. (tabs(1,1) < 273.15) ) then
     ffrozp(1,1) = 1.0
  else
     ffrozp(1,1) = 0.0
  endif

  !
  ! at each time step, using the forcing fields (and t1, the skin temperature, which
  ! gets updated by lsmruc), we need to compute a few additional thermodynamic variables.
  ! ultimately, th2, q2sat and dqsdt2 get passed to sflx;
  !             t1v and th2v get used in sfcdif_off but are not used by lsmruc.
  !             rho is used in lsmruc.
  !

  call caltmp(t1(1,1), tabs(1,1), sfcprs(1,1), zlvl, q2(1,1), th2(1,1), t1v(1,1), th2v(1,1), rho(1,1)) ! returns th2, t1v, th2v, rho
  call calhum(tabs(1,1), sfcprs(1,1), q2sat(1,1), dqsdt2(1,1)) ! returns q2sat, dqsdt2

  !
  ! if the usemonalb flag is .true., we want to provide alb from the user-specified
  ! trend through the year, rather than let lsmruc calculate it for us.
  !
  if (usemonalb .eqv. .true. ) then
     albbck(1,1) = ruc37_month_d(albedo_monthly, nowdate)
  else
     albbck(1,1) = 0.18
  endif
 
  !
  ! if the rdlai2d flag is .true., we want to provide xlai from the user-specified
  ! trend through the year, rather than let lsmruc calculate it for us.
  !
 
  if (rdlai2d .eqv. .true.) then
     xlai(1,1) = ruc37_month_d(lai_monthly, nowdate)
  else
     xlai(1,1) = -9.9999996E+35  
  endif
 
  !
  ! shdfac comes from the user-specified trend through the year.  no other option
  ! at the moment
  !
 
  shdfac(1,1) = ruc37_month_d(shdfac_monthly, nowdate)
  z0(1,1)     = 0.05
!  z0(1,1)     = ruc37_month_d(z0brd_monthly, nowdate) 
  shdmin(1,1) = minval(shdfac_monthly)
  shdmax(1,1) = maxval(shdfac_monthly) 
  !
  ! z0brd is computed within lsmruc.  but we need an initial value, so the call to 
  ! sfcdif_myj can do its thing.  subsequent timesteps will recycle the z0brd 
  ! value as returned from lsmruc in the previous timestep.
  !

  if ( sfcdif_option .eq. 1 ) then
     z0brd(1,1)  = z0(1,1)
  endif
  !
  ! sfcdif_off computes mixing lengths for momentum and heat, cm and ch.
  ! z0 is needed for sfcdif_off.  we use the z0 as computed in the previous
  ! timestep of lsmruc, but on the first time step, we need a value of z0.  this
  ! is set above from our background value.  the initial value may not be quite
  ! what we want, but for that one timestep, it should be ok.  additionally,
  ! ch and cm need some values for the initial timestep.  these values are
  ! set above.
  !
  if ( sfcdif_option .eq. 0 ) then
       call sfcdif_off_ruc37( zlvl , z0(1,1) , t1v(1,1) , th2v(1,1) , sfcspd(1,1) , czil , cm(1,1) , ch(1,1))
  else if ( sfcdif_option .eq. 1 ) then
     call sfcdif_myj ( zlvl, zlvl_wind , z0(1,1) , z0brd(1,1) , sfcprs(1,1) , t1(1,1) , tabs(1,1) , q1(1,1)/(1.+q1(1,1)) , &
          q2(1,1) , sfcspd(1,1) , czil , ribb , cm(1,1) , ch(1,1) , vegtyp(1,1) , isurban , iz0tlnd )
  endif
  
  ! additional forcing variables: 
 
  ! solnet is an additional forcing field, created by applying the albedo to soldn.
  ! albedo is returned each time step from lsmruc.  the initial value is perhaps
  ! not quite what we want, but each subsequent timestep should be ok.
  !
  solnet(1,1) = soldn(1,1) * (1.0-albedo(1,1))
 
  !
  ! apply the emissivity factor to the given longwave radiation.
  ! this takes the emissi value from the previous time step, except
  ! for the first time through the loop, when emissi is set above.
  !
  lwdn(1,1) = longwave(1,1) * emissi(1,1)
 
 
  !tgs as default, let's use frpcpn .eq. .false.. reset frpcpn to .true. if there is data
  ! about amounts of different precipitation types (rain, snow, graupel, ice, etc.)
  frpcpn = .false.

  ! temporary set fractional_seaice=0
  fractional_seaice = 0
  if ( fractional_seaice .eq. 0 ) then
     xice_threshold = 0.5
  else if ( fractional_seaice .eq. 1 ) then
     xice_threshold = 0.02
  endif
  
  ! lis only run land models at non-ice grid boxes 
  ice(1,1) = 0  
  xice(1,1) = ice(1,1)

  if ( sfcdif_option .eq. 1 ) then
     myj = .true.
  endif

  xland(1,1) = 1.
  nlcat=lucats
  nscat=slcats
  ! 
  ! tanya's explanation about mosiac_lu, mosaic_soil, landusef and soilctop 
  ! 
  ! i have an option in ruc lsm to compute grid cell averaged surface properties,
  ! using fractional land use and soil data: landusef and soilctop. this option is
  ! controlled by 2 parameters, and if they are set to 0, then this option is
  ! disabled:
  ! 
  !     mosaic_lu=0 mosaic_soil=0
  ! 
  ! for a point integration this option is not used.  but for 2-d applications
  ! these two arrays should be defined as in wrf.
  mosaic_lu=0
  mosaic_soil=0
  allocate(landusef(nlcat))
  allocate(soilctop(nlcat))
  landusef = 0.0 
  soilctop = 0.0 
 
  
  graupelncv(1,1) = 0.
   
! if( snowf > 0.) then
!    snowncv(1,1) = snowf*dt
!    rainncv(1,1) = rainf*dt
! else
  if(tabs(1,1) > 273.15) then
    snowncv(1,1) = 0.
    rainncv(1,1) = prcp(1,1)*dt
  else
    snowncv(1,1) = prcp(1,1)*dt
    rainncv(1,1) = 0.
  endif
! endif
  mavail(1,1) = soilmois(1,1,1)/maxsmc(nscat)
  q2(1,1)=q2(1,1)/(1.-q2(1,1))

  zlvl3d(1,1,1) = zlvl

  call lsmruc(local,dt,ktime,nsoil,                            &
              sldpth2,prcp*dt,snow,snowh,sncovr,ffrozp,frpcpn, &
              rhosnf,graupelncv,snowncv,rainncv,precipfr,      &
              zlvl3d*2.,sfcprs,tabs,q2,qcatm,rho,                &
              longwave,solnet,emissi,ch,                       &
              mavail,cmc,shdfac*100.,albedo,z0,                &
              z0brd,snoalb, albbck, xlai,                      &
              llanduse, landusef, nlcat, mosaic_lu,            &
              mosaic_soil, soilctop, nscat,                    &
              qmax,qmin,psis,ksat,bclh,qwrtz,wilt,ref,         &
              qsfc,qsg,q1,qcg,dew,soilt1,tsnav,                &
              tbot,vegtyp,soiltyp,xland,                       &
              iswater,isice,xice,xice_threshold,               &
              soilmois,sliqw,soilm,soilw,                      &
              tslb,t1,sheat,etakin,eta,                        &
              edir,ec,ett,et,snflx,fltot,                      &
              runoff1,runoff2,drip,esnow,                      &
              sfcevp,ssoil,snowfallac,acsnow,snomlt,           &
              smfr3d,keepfr3dflag,snthresh,                    &
              myj,shdmin*100.,shdmax*100.,rdlai2d,             &
              1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1 )

  !
  ! residual of surface energy balance equation terms
  !
  res(1,1)=fltot(1,1)
  q2(1,1)=q2(1,1)/(1.+q2(1,1))
  
  ! convert drip from m/timestep to kg m{-2} s{-1} (mm s-1)
  drip(1,1) = 1.e3 * drip(1,1) / dt

  !!! state out
  emiss1   = emissi(1,1) 
  ch_1      = ch(1,1)    
  cm_1      = cm(1,1)    
  sneqv1   = snow(1,1)*1e-3  ! kg/m^-2 --> [m]
  snowh1   = snowh(1,1)      ! m
  snowc1   = sncovr(1,1)
  alb1     = albedo(1,1)
  albbck1  = albbck(1,1)
  snt75cm1 = soilt1(1,1)
  tskin1   = t1(1,1)    
  qvg1     = q1(1,1)
  qsfc1    = qsfc(1,1)
  qcg1     = qcg(1,1)   
  qsg1     = qsg(1,1)   
  tsnav1   = tsnav(1,1) + 273.15 
  soilm1   = soilm(1,1) 
  smroot1  = soilw(1,1) 
  canwat1  = cmc(1,1)   
  shdfac1  = shdfac(1,1)
  shdmin1  = shdmin(1,1)
  shdmax1  = shdmax(1,1)
  qmax2    = qmax(1,1)
  qmin2    = qmin(1,1)
  psis1    = psis(1,1)
  ksat1    = ksat(1,1)
  bclh1    = bclh(1,1)
  qwrtz1   = qwrtz(1,1)
  wilt1    = wilt(1,1)
  ref1     = ref(1,1)
  
  do k=1,nsoil
    smc2(k)  = soilmois (1,k,1)
    sh2o2(k) = sliqw    (1,k,1)
    stc2(k)  = tslb     (1,k,1)
    smfr3(k) = smfr3d   (1,k,1)
    keepfr3(k)= keepfr3dflag(1,k,1)
  enddo
 
  !!! output 
  qh1         = sheat(1,1)  ! W m-2
  qle1        = eta(1,1)    ! W m-2 
  eta1        = etakin(1,1) ! Kg/m2s
  qg1         = ssoil(1,1) 
  rhosnf1     = rhosnf(1,1) 
  precipfr1   = precipfr(1,1)
  snowfallac1 = snowfallac(1,1) 
  acsnow1     = acsnow(1,1)
  sfcevp1     = sfcevp(1,1) 
  snomlt1     = snomlt(1,1) 
  snthresh1   = snthresh(1,1)  ! [m]
  dew1        = dew(1,1) 
  
  drip1       = drip(1,1)
  qs1         = runoff1(1,1) * 1000.0 ! from m s-1 to kg m-2 s-1 
  qsb1        = runoff2(1,1) * 1000.0 ! from m s-1 to kg m-2 s-1 
  snflx1      = snflx(1,1)

  edir1       = edir(1,1) 
  ec1         = ec(1,1) 
  ett1        = ett(1,1) 
  esnow1      = esnow(1,1) 


  deallocate(landusef)
  deallocate(soilctop)
  deallocate(soilmois)
  deallocate(sliqw )
  deallocate(smfr3d)
  deallocate(keepfr3dflag)
  deallocate(tslb)
  deallocate(et)
  deallocate(smav)
end subroutine ruc_driver_37

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

real function ruc37_month_d(a12, nowdate) result (nowval)
  !
  ! given a set of 12 values, taken to be valid on the fifteenth of each month (jan through dec)
  ! and a date in the form <yyyymmdd[hhmmss]> ....
  ! 
  ! return a value valid for the day given in <nowdate>, as an interpolation from the 12
  ! monthly values.
  !
  use kwm_date_utilities_ruc37
  implicit none
  real, dimension(12), intent(in) :: a12 ! 12 monthly values, taken to be valid on the 15th of
  !                                      ! the month
  character(len=12), intent(in) :: nowdate ! date, in the form <yyyymmdd[hhmmss]>
  integer :: nowy, nowm, nowd
  integer :: prevm, postm
  real    :: factor
  integer, dimension(12) :: ndays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  !
  ! handle leap year by setting the number of days in february for the year in question.
  !
  read(nowdate(1:8),'(i4,i2,i2)') nowy, nowm, nowd
  ndays(2) = nfeb(nowy)

  !
  ! do interpolation between the fifteenth of two successive months.
  !
  if (nowd .eq. 15) then
     nowval = a12(nowm)
     return
  else if (nowd < 15) then
     postm = nowm
     prevm = nowm - 1
     if (prevm .eq. 0) prevm = 12
     factor = real(ndays(prevm)-15+nowd)/real(ndays(prevm))
  else if (nowd > 15) then
     prevm = nowm
     postm = nowm + 1
     if (postm .eq. 13) postm = 1
     factor = real(nowd-15)/real(ndays(prevm))
  endif

  nowval = a12(prevm)*(1.0-factor) + a12(postm)*factor

end function ruc37_month_d

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! adapted from the wrf subroutine in module_sf_ruclsm.f:
!-----------------------------------------------------------------subroutine ruclsm_soilvegprm( landuse_tbl, soil_tbl, mminlu, mminsl)
subroutine ruc37_soilvegprm(mminlu, mminsl, landuse_tbl, soil_tbl, gen_tbl )

  use module_sf_ruclsm_ruc37
!  use RUC37_lsmMod
  implicit none
  
  character(len=*), intent(in) :: landuse_tbl, soil_tbl, gen_tbl 
  character(len=*), intent(in) :: mminlu, mminsl
  integer :: lumatch, iindex, lc, num_slope
  integer :: ierr
  integer , parameter :: open_ok = 0

  character*128 :: mess , message
  logical :: wrf_dm_on_monitor


!-----specify vegetation related characteristics :
!             albbck: sfc albedo (in percentage)
!       note: recommend to use monthly climatological albedo 
!                 z0: roughness length (m)
!            lemitbl: emmisivity
!             pctbl : plant coefficient
!            shdtbl : shading factor
!       note: recommend to use monthly climatological greeness fractioin 
!            ifortbl: parameter to define rootiung depth
!              rsmin: mimimum stomatal resistance (s m-1) - not used
!                rgl: parameters used in radiation stress function - not used
!                 hs: parameter used in vapor pressure deficit functio
!               snup: threshold snow depth (in water equivalent m) that - not used
!             laitbl: leaf area index (dimensionless)
!       note: recommend to use monthly climatological lai
!             maxalb: upper bound on maximum albedo over deep snow
!       note: recommend to use monthly climatological data of max snow albedo
!               topt: optimum transpiration air temperature. (k) - not used
!             cmcmax: maximum canopy water capacity
!             cmxtbl: max cnpy capacity (m)
!             cfactr: parameter used in the canopy inteception calculati
!                     implies 100% snow cover - not used
!              rsmax: max. stomatal resistance (s m-1) - not used
!               bare: category number for bare soil
!            natural: category number for grassland
!               crop: category number for crops
!
!-----read in vegetaion properties from vegparm.tbl
!

  if ( wrf_dm_on_monitor() ) then

     open(19, file=trim(landuse_tbl),form='formatted',status='old',iostat=ierr)
     if(ierr .ne. open_ok ) then
        write(message,fmt='(a)') &
             'module_sf_ruclsm.f: soil_veg_gen_parm: failure opening ', trim(landuse_tbl)
        call wrf_error_fatal ( message )
     end if


     lumatch=0

     find_lutype : do while (lumatch .eq. 0)
        read (19,*,end=2002)
        read (19,*,end=2002)lutype
        read (19,*)lucats,iindex

        if(lutype.eq.mminlu)then
           write( mess , * ) 'landuse type = ' // trim ( lutype ) // ' found', lucats,' categories'
           call wrf_message( mess )
           lumatch=1
        else
           call wrf_message ( "skipping over lutype = " // trim ( lutype ) )
           do lc = 1, lucats+16
              read(19,*)
           enddo
        endif
     enddo find_lutype
! prevent possible array overwrite, bill bovermann, ibm, may 6, 2008
     if ( size(shdtbl)       < lucats .or. &
          size(rstbl)        < lucats .or. &
          size(rgltbl)       < lucats .or. &
          size(hstbl)        < lucats .or. &
          size(snuptbl)      < lucats .or. &
          size(maxalb)       < lucats .or. &
          size(laitbl)       < lucats .or. & 
          size(z0tbl)        < lucats .or. &
          size(albtbl)       < lucats .or. &
          size(lemitbl )     < lucats ) then
        call wrf_error_fatal('table sizes too small for value of lucats in module_sf_ruclsm.f')
     endif

     if(lutype.eq.mminlu)then
        do lc=1,lucats
           read (19,*)iindex,albtbl(lc),z0tbl(lc),         &
                 lemitbl(lc),pctbl(lc), shdtbl(lc),        &
                 ifortbl(lc),rstbl(lc),rgltbl(lc),         &
                 hstbl(lc),snuptbl(lc),laitbl(lc),maxalb(lc)

        enddo
!
        read (19,*)
        read (19,*)topt_data
        read (19,*)
        read (19,*)cmcmax_data
        read (19,*)
        read (19,*)cfactr_data
        read (19,*)
        read (19,*)rsmax_data
        read (19,*)
        read (19,*)bare
        read (19,*)
        read (19,*)natural
        read (19,*) 
        read (19,*)crop
        read (19,*)
        read (19,*)urban
     endif
!
2002 continue

     close (19)
     if (lumatch .eq. 0) then
        call wrf_error_fatal ("land use dataset '"//mminlu//"' not found in "//trim(landuse_tbl)//".")
     endif
  endif


!
!-----read in soil properties from soilparm.tbl
!
  if ( wrf_dm_on_monitor() ) then
     open(19, file=trim(soil_tbl),form='formatted',status='old',iostat=ierr)
     if(ierr .ne. open_ok ) then
        write(message,fmt='(a)') &
             'module_sf_ruclsm.f: soil_veg_gen_parm: failure opening ', trim(soil_tbl)
        call wrf_error_fatal ( message )
     end if

     write(mess,*) 'input soil texture classificaion = ', trim ( mminsl )
     call wrf_message( mess )

     lumatch=0

     read (19,*)
     read (19,2000,end=2003)sltype
2000 format (a8)
     read (19,*)slcats,iindex
        if(sltype.ne.mminsl)then
          do lc=1,slcats
              read (19,*) iindex,bb(lc),drysmc(lc),hc(lc),maxsmc(lc),&
                        refsmc(lc),satpsi(lc),satdk(lc), satdw(lc),   &
                        wltsmc(lc), qtz(lc)
          enddo
        endif

!     read (19,*)
!     read (19,2000,end=2003)sltype
!     read (19,*)slcats,iindex

!print *,'sltype=',sltype,'slcats=',slcats
 
     if(sltype.eq.mminsl)then
        write( mess , * ) 'soil texture classification = ', trim ( sltype ) , ' found', &
             slcats,' categories'
        call wrf_message ( mess )
        lumatch=1
     endif
!   enddo find_lutype

!      print *,'soil lumatch=',lumatch
! prevent possible array overwrite, bill bovermann, ibm, may 6, 2008
     if ( size(bb    ) < slcats .or. &
          size(drysmc) < slcats .or. &
          size(hc   )  < slcats .or. &
          size(maxsmc) < slcats .or. &
          size(refsmc) < slcats .or. &
          size(satpsi) < slcats .or. &
          size(satdk ) < slcats .or. &
          size(satdw ) < slcats .or. &
          size(wltsmc) < slcats .or. &
          size(qtz   ) < slcats  ) then
        call wrf_error_fatal('table sizes too small for value of slcats in module_sf_rucdrv.f')
     endif
     if(sltype.eq.mminsl)then
        do lc=1,slcats
           read (19,*) iindex,bb(lc),drysmc(lc),hc(lc),maxsmc(lc),   &
                       refsmc(lc),satpsi(lc),satdk(lc), satdw(lc),   &
                       wltsmc(lc), qtz(lc)
        enddo
     endif

2003 continue

     close (19)
  endif

  if(lumatch.eq.0)then
     call wrf_message( 'soil texture in input file does not ' )
     call wrf_message( 'match soilparm table'                 )
     call wrf_error_fatal ( 'inconsistent or missing soilparm file' )
  endif

!
!-----read in general parameters from genparm.tbl
!
  if ( wrf_dm_on_monitor() ) then
     open(19, file=trim(gen_tbl),form='formatted',status='old',iostat=ierr)
     if(ierr .ne. open_ok ) then
        write(message,fmt='(a)') &
             'module_sf_ruclsm.f: soil_veg_gen_parm: failure opening ', trim(gen_tbl)
        call wrf_error_fatal ( message )
     end if

     read (19,*)
     read (19,*)
     read (19,*) num_slope

     slpcats=num_slope
! prevent possible array overwrite, bill bovermann, ibm, may 6, 2008
     if ( size(slope_data) < num_slope ) then
        call wrf_error_fatal('num_slope too large for slope_data array in module_sf_ruclsm')
     endif

     do lc=1,slpcats
        read (19,*)slope_data(lc)
     enddo

     read (19,*)
     read (19,*)sbeta_data
     read (19,*)
     read (19,*)fxexp_data
     read (19,*)
     read (19,*)csoil_data
     read (19,*)
     read (19,*)salp_data
     read (19,*)
     read (19,*)refdk_data
     read (19,*)
     read (19,*)refkdt_data
     read (19,*)
     read (19,*)frzk_data
     read (19,*)
     read (19,*)zbot_data
     read (19,*)
     read (19,*)czil_data
     read (19,*)
     read (19,*)smlow_data
     read (19,*)
     read (19,*)smhigh_data
     read (19,*)
     read (19,*)lvcoef_data
     close (19)
  endif
!-----------------------------------------------------------------
end subroutine ruc37_soilvegprm
!-----------------------------------------------------------------
  
subroutine sfcdif_off_ruc37 (zlm,z0,thz0,thlm,sfcspd,czil,akms,akhs)

! ----------------------------------------------------------------------
! subroutine sfcdif (renamed sfcdif_off to avoid clash with eta pbl)
! ----------------------------------------------------------------------
! calculate surface layer exchange coefficients via iterative process.
! see chen et al (1997, blm)
! ----------------------------------------------------------------------

      implicit none
      real     wwst, wwst2, g, vkrm, excm, beta, btg, elfc, wold, wnew
      real     pihf, epsu2, epsust, epsit, epsa, ztmin, ztmax, hpbl,     &
     & sqvisc
      real     ric, rric, fhneu, rfc, rfac, zz, pslmu, pslms, pslhu,     &
     & pslhs
      real     xx, pspmu, yy, pspms, psphu, psphs, zlm, z0, thz0, thlm
      real     sfcspd, czil, akms, akhs, zilfc, zu, zt, rdz, cxch
      real     dthv, du2, btgh, wstar2, ustar, zslu, zslt, rlogu, rlogt
      real     rlmo, zetalt, zetalu, zetau, zetat, xlu4, xlt4, xu4, xt4
!cc   ......real ztfc

      real     xlu, xlt, xu, xt, psmz, simm, pshz, simh, ustark, rlmn,  &
     &         rlma

      integer  itrmx, ilech, itr
      parameter                                                         &
     &        (wwst = 1.2,wwst2 = wwst * wwst,g = 9.8,vkrm = 0.40,      &
     &         excm = 0.001                                             &
     &        ,beta = 1./270.,btg = beta * g,elfc = vkrm * btg          &
     &                  ,wold =.15,wnew = 1. - wold,itrmx = 05,         &
     &                   pihf = 3.14159265/2.)
      parameter                                                         &
     &         (epsu2 = 1.e-4,epsust = 0.07,epsit = 1.e-4,epsa = 1.e-8  &
     &         ,ztmin = -5.,ztmax = 1.,hpbl = 1000.0                    &
     &          ,sqvisc = 258.2)
      parameter                                                         &
     &       (ric = 0.183,rric = 1.0/ ric,fhneu = 0.8,rfc = 0.191       &
     &        ,rfac = ric / (fhneu * rfc * rfc))

! ----------------------------------------------------------------------
! note: the two code blocks below define functions
! ----------------------------------------------------------------------
! lech's surface functions
! ----------------------------------------------------------------------
      pslmu (zz)= -0.96* log (1.0-4.5* zz)
      pslms (zz)= zz * rric -2.076* (1. -1./ (zz +1.))
      pslhu (zz)= -0.96* log (1.0-4.5* zz)

! ----------------------------------------------------------------------
! paulson's surface functions
! ----------------------------------------------------------------------
      pslhs (zz)= zz * rfac -2.076* (1. -1./ (zz +1.))
      pspmu (xx)= -2.* log ( (xx +1.)*0.5) - log ( (xx * xx +1.)*0.5)   &
     &        +2.* atan (xx)                                            &
     &- pihf
      pspms (yy)= 5.* yy
      psphu (xx)= -2.* log ( (xx * xx +1.)*0.5)

! ----------------------------------------------------------------------
! this routine sfcdif can handle both over open water (sea, ocean) and
! over solid surface (land, sea-ice).
! ----------------------------------------------------------------------
      psphs (yy)= 5.* yy

! ----------------------------------------------------------------------
!     ztfc: ratio of zoh/zom  less or equal than 1
!     c......ztfc=0.1
!     czil: constant c in zilitinkevich, s. s.1995,:note about zt
! ----------------------------------------------------------------------
      ilech = 0

! ----------------------------------------------------------------------
      zilfc = - czil * vkrm * sqvisc
!     c.......zt=z0*ztfc
      zu = z0
      rdz = 1./ zlm
      cxch = excm * rdz
      dthv = thlm - thz0

! ----------------------------------------------------------------------
! beljars correction of ustar
! ----------------------------------------------------------------------
      du2 = max (sfcspd * sfcspd,epsu2)
!cc   if statements to avoid tangent linear problems near zero
      btgh = btg * hpbl
      if (btgh * akhs * dthv .ne. 0.0) then
         wstar2 = wwst2* abs (btgh * akhs * dthv)** (2./3.)
      else
         wstar2 = 0.0
      end if

! ----------------------------------------------------------------------
! zilitinkevitch approach for zt
! ----------------------------------------------------------------------
      ustar = max (sqrt (akms * sqrt (du2+ wstar2)),epsust)

! ----------------------------------------------------------------------
      zt = exp (zilfc * sqrt (ustar * z0))* z0
      zslu = zlm + zu
!     print*,'zslt=',zslt
!     print*,'zlm=',zlm
!     print*,'zt=',zt

      zslt = zlm + zt
      rlogu = log (zslu / zu)

      rlogt = log (zslt / zt)
!     print*,'rlmo=',rlmo
!     print*,'elfc=',elfc
!     print*,'akhs=',akhs
!     print*,'dthv=',dthv
!     print*,'ustar=',ustar

      rlmo = elfc * akhs * dthv / ustar **3
! ----------------------------------------------------------------------
! 1./monin-obukkhov length-scale
! ----------------------------------------------------------------------
      do itr = 1,itrmx
         zetalt = max (zslt * rlmo,ztmin)
         rlmo = zetalt / zslt
         zetalu = zslu * rlmo
         zetau = zu * rlmo

         zetat = zt * rlmo
         if (ilech .eq. 0) then
            if (rlmo .lt. 0.)then
               xlu4 = 1. -16.* zetalu
               xlt4 = 1. -16.* zetalt
               xu4 = 1. -16.* zetau

               xt4 = 1. -16.* zetat
               xlu = sqrt (sqrt (xlu4))
               xlt = sqrt (sqrt (xlt4))
               xu = sqrt (sqrt (xu4))

               xt = sqrt (sqrt (xt4))
!     print*,'-----------1------------'
!     print*,'psmz=',psmz
!     print*,'pspmu(zetau)=',pspmu(zetau)
!     print*,'xu=',xu
!     print*,'------------------------'
               psmz = pspmu (xu)
               simm = pspmu (xlu) - psmz + rlogu
               pshz = psphu (xt)
               simh = psphu (xlt) - pshz + rlogt
            else
               zetalu = min (zetalu,ztmax)
               zetalt = min (zetalt,ztmax)
!     print*,'-----------2------------'
!     print*,'psmz=',psmz
!     print*,'pspms(zetau)=',pspms(zetau)
!     print*,'zetau=',zetau
!     print*,'------------------------'
               psmz = pspms (zetau)
               simm = pspms (zetalu) - psmz + rlogu
               pshz = psphs (zetat)
               simh = psphs (zetalt) - pshz + rlogt
            end if
! ----------------------------------------------------------------------
! lech's functions
! ----------------------------------------------------------------------
         else
            if (rlmo .lt. 0.)then
!     print*,'-----------3------------'
!     print*,'psmz=',psmz
!     print*,'pslmu(zetau)=',pslmu(zetau)
!     print*,'zetau=',zetau
!     print*,'------------------------'
               psmz = pslmu (zetau)
               simm = pslmu (zetalu) - psmz + rlogu
               pshz = pslhu (zetat)
               simh = pslhu (zetalt) - pshz + rlogt
            else
               zetalu = min (zetalu,ztmax)

               zetalt = min (zetalt,ztmax)
!     print*,'-----------4------------'
!     print*,'psmz=',psmz
!     print*,'pslms(zetau)=',pslms(zetau)
!     print*,'zetau=',zetau
!     print*,'------------------------'
               psmz = pslms (zetau)
               simm = pslms (zetalu) - psmz + rlogu
               pshz = pslhs (zetat)
               simh = pslhs (zetalt) - pshz + rlogt
            end if
! ----------------------------------------------------------------------
! beljaars correction for ustar
! ----------------------------------------------------------------------
         end if

! ----------------------------------------------------------------------
! zilitinkevitch fix for zt
! ----------------------------------------------------------------------
         ustar = max (sqrt (akms * sqrt (du2+ wstar2)),epsust)

         zt = exp (zilfc * sqrt (ustar * z0))* z0
         zslt = zlm + zt
!-----------------------------------------------------------------------
         rlogt = log (zslt / zt)
         ustark = ustar * vkrm
         akms = max (ustark / simm,cxch)
!-----------------------------------------------------------------------
! if statements to avoid tangent linear problems near zero
!-----------------------------------------------------------------------
         akhs = max (ustark / simh,cxch)
         if (btgh * akhs * dthv .ne. 0.0) then
            wstar2 = wwst2* abs (btgh * akhs * dthv)** (2./3.)
         else
            wstar2 = 0.0
         end if
!-----------------------------------------------------------------------
         rlmn = elfc * akhs * dthv / ustar **3
!-----------------------------------------------------------------------
!     if(abs((rlmn-rlmo)/rlma).lt.epsit)    go to 110
!-----------------------------------------------------------------------
         rlma = rlmo * wold+ rlmn * wnew
!-----------------------------------------------------------------------
         rlmo = rlma
!     print*,'----------------------------'
!     print*,'sfcdif output !  ! ! ! ! ! ! ! !  !   !    !'

!     print*,'zlm=',zlm
!     print*,'z0=',z0
!     print*,'thz0=',thz0
!     print*,'thlm=',thlm
!     print*,'sfcspd=',sfcspd
!     print*,'czil=',czil
!     print*,'akms=',akms
!     print*,'akhs=',akhs
!     print*,'----------------------------'

      end do
! ----------------------------------------------------------------------
  end subroutine sfcdif_off_ruc37


! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------

!subroutine wrf_error_fatal(string)
!  implicit none
!  character(len=*), intent(in) :: string
!  print*, string
!  stop
!end subroutine wrf_error_fatal

! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------

!subroutine wrf_message(msg)
!  implicit none
!  character(len=*), intent(in) :: msg
!  write(*,'(A)') msg
!end subroutine wrf_message

! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------

logical function wrf_dm_on_monitor() result (return_value)
  implicit none
  return_value = .TRUE.
end function wrf_dm_on_monitor

! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------

subroutine wrf_dm_bcast_real(rval, ival)
  implicit none
  real, intent(in) :: rval
  integer, intent(in) :: ival
end subroutine wrf_dm_bcast_real

! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------

subroutine wrf_dm_bcast_integer(ival1, ival2)
  implicit none
  real, intent(in) :: ival1
  integer, intent(in) :: ival2
end subroutine wrf_dm_bcast_integer

! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------

subroutine wrf_dm_bcast_string(sval, ival)
  implicit none
  character(len=*), intent(in) :: sval
  integer, intent(in) :: ival
end subroutine wrf_dm_bcast_string

! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------
