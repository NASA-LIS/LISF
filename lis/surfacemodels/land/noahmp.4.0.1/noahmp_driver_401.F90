!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Oct 15 2018 Shugong Wang started for implementing Noah-MP 4.0.1 based on the version of 3.6  
! Oct 15 2018 Zhuo Wang modifed for implementing Noah-MP 4.0.1 

#undef LIS_NoahMP_TEST 
#undef WRF_HYDRO  
! !INTERFACE
subroutine noahmp_driver_401(n, ttile, itimestep, &  
                            latitude, longitude,                                        &
                            year    , month   , day     , hour    , minute  ,           &
                            dz8w    ,                                                   & ! new in : model configuration
                            dt      , sldpth  , nsoil   , nsnow   ,                     & ! in : model configuration  
                            vegetype, soiltype, shdfac_monthly    , tbot    ,           & ! in : Vegetation/Soil characteristics 
                            urban_vegetype,                                             & ! in
                            cropcat,  planting, harvest ,season_gdd,                    & ! in : Vegetation/Soil characteristics
                            dveg_opt, crs_opt, btr_opt, run_opt, sfc_opt, frz_opt,      & ! in : User options
                            inf_opt, rad_opt, alb_opt , snf_opt, tbot_opt, stc_opt,     & ! in : User options
                            gla_opt, rsf_opt, soil_opt, pedo_opt, crop_opt, iz0tlnd   , & ! in : new options
                            urban_opt,                                                  & ! in : new options
                            soilcomp, soilcL1, soilcL2, soilcL3, soilcL4,               & ! in : new options
                            tair    , psurf   , wind_e   , wind_n   , qair    ,         & ! in : forcing
                            swdown  , lwdown  , prcp    ,                               & ! in : forcing
                            tsk     , hfx     , qfx     , lh      , grdflx   ,          & ! in/out LSM eqv
                            sfcrunoff, udrunoff, albedo , qsnbot  , subsnow  ,          & ! in/out LSM eqv
                            snowc   , smc     ,     pah ,                               & ! in/out LSM eqv
                            sh2o    , tslb    , sneqv   , snowh   , canwat  , acsnom  , & ! in/out LSM eqv
                            acsnow  , emiss   , rs      ,                               & ! in/out LSM eqv
                            isnow   , tv      , tg      , canice  , canliq  , eah     , & ! in/out Noah MP only
                            tah     , cm      , ch      , fwet    , sneqvo  , albold  , & ! in/out Noah MP only
                            qsnow   , wslake  , zwt     , wa      , wt      , tsno    , & ! in/out Noah MP only
                            zss     , snowice , snowliq , lfmass  , rtmass  , stmass  , & ! in/out Noah MP only
                            wood    , stblcp  , fastcp  , lai     , sai     , tauss   , & ! in/out Noah MP only
                            smoiseq , smcwtd  ,deeprech , rech    ,                     & ! in/out Noah MP only
                            grain   , gdd     , pgs     ,                               & ! in/out Noah MP only for crop model
                            gecros_state,                                               & ! in/out gecros model
                            t2mv    , t2mb    , q2mv    , q2mb    ,                     & ! out Noah MP only
                            trad    , nee     , gpp     , npp     , fveg    , runsf   , & ! out Noah MP only
                            runsb   , ecan    , edir    , etran   ,                     & ! out Noah MP only
                            rainf   , snowf   , fsa     , fira    ,                     & ! out Noah MP only
                            apar    , psn     , sav     , sag     , rssun   , rssha   , & ! out Noah MP only
                            bgap    , wgap    , tgb     , tgv     , chv     , chb     , & ! out Noah MP only
                            shg     , shc     , shb     , evg     , evb     , ghv     , & ! out Noah MP only
                            ghb     , irg     , irc     , irb     , tr      , evc     , & ! out Noah MP only
                            chleaf  , chuc    , chv2    , chb2                          ) ! out Noah MP only

  use module_sf_noahmpdrv_401, only: noahmplsm_401
  use LIS_coreMod, only    : LIS_rc
  use LIS_logMod,  only    : LIS_logunit, LIS_endrun
  use LIS_timeMgrMod, only : LIS_date2time, LIS_tick

  implicit none
  integer, intent(in) :: n                    ! nest id 
  integer, intent(in) :: ttile                ! tile id
  integer, intent(in) :: itimestep            ! timestep number
  real,    intent(in) :: latitude             ! latitude in decimal degree [-]
  real,    intent(in) :: longitude            ! longitude in decimal degree [-]
  integer, intent(in) :: year                 ! year of the current time step [-]
  integer, intent(in) :: month                ! month of the current time step [-] 
  integer, intent(in) :: day                  ! day of the current time step [-]
  integer, intent(in) :: hour                 ! hour of the current time step [-] 
  integer, intent(in) :: minute               ! minute of the current time step [-]
  real,    intent(in) :: dt                   ! timestep [s] 
  
  ! Revised by Zhuo Wang and Shugong Wang
  real,    intent(in) :: dz8w                 ! thickness of atmo layers [m]
  
  integer, intent(in) :: nsoil                ! number of soil layers [-] 
  real,    intent(in) :: sldpth(nsoil)        ! thickness of soil layers [m]  
  integer, intent(in) :: nsnow                ! maximum number of snow layers (e.g. 3) [-] 
  integer, intent(in) :: vegetype             ! vegetation type [-] 
  integer, intent(in) :: soiltype             ! soil type [-] 
  integer, intent(in) :: urban_vegetype       ! urban land cover type index [-] 
  real,    intent(in) :: shdfac_monthly(12)   ! monthly values for green vegetation fraction [-]
  real,    intent(in) :: tbot                 ! deep soil temperature [K]

  ! Crop Model
  integer, intent(in) :: cropcat              ! crop catagory
  real,    intent(in) :: planting             ! planting date
  real,    intent(in) :: harvest              ! harvest date
  real,    intent(in) :: season_gdd           ! growing season GDD
  real,    intent(inout) :: gdd               ! growing degree days XING (based on 10C)
  real,    intent(inout) :: grain             ! mass of grain XING [g/m2]
  integer,    intent(inout) :: pgs

  ! gecros model
  real,    intent(inout) :: gecros_state(60)  !  gecros crop

  integer, intent(in) :: dveg_opt             ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
                                              ! vegetation model ( 1->prescribed [table LAI, shdfac=FVEG]; 2->dynamic; 3->table LAI, calculate FVEG 4->table LAI, shdfac=maximum)
  integer, intent(in) :: crs_opt              ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  integer, intent(in) :: btr_opt              ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  integer, intent(in) :: run_opt              ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  integer, intent(in) :: sfc_opt              ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  integer, intent(in) :: frz_opt              ! supercooled liquid water (1-> NY06; 2->Koren99)
  integer, intent(in) :: inf_opt              ! frozen soil permeability (1-> NY06; 2->Koren99)
  integer, intent(in) :: rad_opt              ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1--Fveg)
  integer, intent(in) :: alb_opt              ! snow surface albedo (1->BATS; 2->CLASS)
  integer, intent(in) :: snf_opt              ! precipitation partitioning between snow and rain
                                              ! (1-Jordan91; 2->BATS: Snow when SFCTMP < TFRZ+2.2;
                                              !  3->Noah: Snow when SFCTMP < TFRZ;
                                              !  4->Use WRF precipitation partitioning )
  integer, intent(in) :: tbot_opt             ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
  integer, intent(in) :: stc_opt              ! snow/soil temperature time scheme (1->semi-implicit; 2->fully implicit)

  ! Added by Zhuo Wang and Shugong Wang
  integer, intent(in) :: gla_opt              ! glacier option (1->phase change; 2->simple)
  integer, intent(in) :: rsf_opt              ! surface resistance (1->Sakaguchi/Zeng; 2->Seller; 3->mod Sellers; 4->1+snow)
  integer, intent(in) :: soil_opt             ! soil configuration option
  integer, intent(in) :: pedo_opt             ! soil pedotransfer function option
  integer, intent(in) :: crop_opt             ! crop model option (0->none; 1->Liu et al.; 2->Gecros)
  integer, intent(in) :: iz0tlnd              ! option of Chen adjustment of Czil (not used)
  integer, intent(in) :: urban_opt            ! urban physics option
  real, intent(in) :: soilcomp(8)             ! soil sand and clay percentage
  real, intent(in) :: soilcL1                 ! soil texture in layer 1
  real, intent(in) :: soilcL2                 ! soil texture in layer 2
  real, intent(in) :: soilcL3                 ! soil texture in layer 3 
  real, intent(in) :: soilcL4                 ! soil texture in layer 4

  real,    intent(in) :: tair                 ! air temperature [K] 
  real,    intent(in) :: psurf                ! air pressure [Pa] 
  real,    intent(in) :: wind_e               ! U wind component [m s-1] 
  real,    intent(in) :: wind_n               ! V wind component [m s-1] 

  real,    intent(in) :: qair                 ! specific humidity [kg/kg]; mv/(mv+md)

  real,    intent(in) :: swdown               ! downward solar radiation [w/m2]
  real,    intent(in) :: lwdown               ! downward longwave radiation [w/m2]
  real,    intent(in) :: prcp                 ! total precipitatin (rainfall+snowfall) [mm]   

  ! Modified by Zhuo Wang and Shugong Wang
  ! in/out (with generic LSM equivalent)
  real, intent(out) :: tsk                  ! surface radiative temperature [K]
  real, intent(out) :: hfx                  ! sensible heat flux [W m-2]
  real, intent(out) :: qfx                  ! latent heat flux [kg s-1 m-2] 
  real, intent(out) :: lh                   ! latent heat flux [W m-2]
  real, intent(out) :: grdflx               ! ground/snow heat flux [W m-2]
  real, intent(inout) :: sfcrunoff            ! accumulated surface runoff [m]
  real, intent(inout) :: udrunoff             ! accumulated sub-surface runoff [m]
  real, intent(out) :: albedo               ! total grid albedo []
  real, intent(out) :: qsnbot               ! snowmelt out the bottom layer [kg s-1 m-2]
  real, intent(out) :: subsnow              ! snow sublimation [kg s-1 m-2]
  real, intent(out) :: pah                  ! precipitation advected heat - total (W/m2)
  real, intent(out) :: snowc                ! snow cover fraction []
  real, intent(inout) :: smc(nsoil)           ! volumetric soil moisture, ice + liquid [m3/m3]
  real, intent(inout) :: sh2o(nsoil)          ! volumetric liquid soil moisture [m3/m3]
  real, intent(inout) :: tslb(nsoil)          ! soil temperature [K] 
  real, intent(inout) :: sneqv                ! snow water equivalent [mm]
  real, intent(inout) :: snowh                ! physical snow depth [m]
  real, intent(inout) :: canwat               ! total canopy water + ice [mm]
  real, intent(inout) :: acsnom               ! accumulated snow melt leaving pack
  real, intent(inout) :: acsnow               ! accumulated snow on grid
  real, intent(out) :: emiss                ! surface bulk emissivity
  real, intent(out) :: rs                   ! Total stomatal resistance (s/m)   
  ! In module_sf_noahmpdrv_401.F, it is defined as inout variable (with generic LSM equivalent), but it is also in the
  ! output list as NoahMP out only. Which one is correct??????/

  ! in/out Noah MP only
  integer, intent(inout) :: isnow             ! actual no. of snow layers 
  real, intent(inout) :: tv                   ! vegetation leaf temperature
  real, intent(inout) :: tg                   ! bulk ground surface temperature
  real, intent(inout) :: canice               ! canopy-intercepted ice (mm)
  real, intent(inout) :: canliq               ! canopy-intercepted liquid water (mm)
  real, intent(inout) :: eah                  ! canopy air vapor pressure (pa)
  real, intent(inout) :: tah                  ! canopy air temperature (k)
  real, intent(inout) :: cm                   ! bulk momentum drag coefficient
  real, intent(inout) :: ch                   ! bulk sensible heat exchange coefficient
  real, intent(inout) :: fwet                 ! wetted or snowed fraction of the canopy (-)
  real, intent(inout) :: sneqvo               ! snow mass at last time step(mm h2o)
  real, intent(inout) :: albold               ! snow albedo at last time step (-)
  real, intent(inout) :: qsnow                ! snowfall on the ground [mm/s]
  real, intent(inout) :: wslake               ! lake water storage [mm]  
  real, intent(inout) :: zwt                  ! water table depth [m]
  real, intent(inout) :: wa                   ! water in the "aquifer" [mm]
  real, intent(inout) :: wt                   ! groundwater storage [mm]
  real, intent(inout) :: tsno(nsnow)          ! snow temperature [K]
  real, intent(inout) :: zss(nsnow+nsoil)     ! snow/soil layer-bottom depth from snow surface [m] 
  real, intent(inout) :: snowice(nsnow)       ! snow layer ice [mm]      
  real, intent(inout) :: snowliq(nsnow)       ! snow layer liquid water [mm]
  real, intent(inout) :: lfmass               ! leaf mass [g/m2]
  real, intent(inout) :: rtmass               ! mass of fine roots [g/m2]
  real, intent(inout) :: stmass               ! stem mass [g/m2]
  real, intent(inout) :: wood                 ! mass of wood (incl. woody roots) [g/m2]
  real, intent(inout) :: stblcp               ! stable carbon in deep soil [g/m2]
  real, intent(inout) :: fastcp               ! short-lived carbon, shallow soil [g/m2]
  real, intent(inout) :: lai                  ! leaf area index
  real, intent(inout) :: sai                  ! stem area index
  real, intent(inout) :: tauss                ! snow age factor
  real, intent(inout) :: smoiseq(nsoil)       ! eq volumetric soil moisture [m3/m3]
  real, intent(inout) :: smcwtd               ! soil moisture content in the layer to the water table when deep
  real, intent(inout) :: deeprech             ! recharge to the water table when deep
  real, intent(inout) :: rech                 ! recharge to the water table (diagnostic)

  ! OUT (with no Noah LSM equivalent)
  real, intent(out) :: t2mv                   ! 2m temperature of vegetation part
  real, intent(out) :: t2mb                   ! 2m temperature of bare ground part
  real, intent(out) :: q2mv                   ! 2m mixing ratio of vegetation part
  real, intent(out) :: q2mb                   ! 2m mixing ratio of bare ground part
  real, intent(out) :: trad                   ! surface radiative temperature (k)
  real, intent(out) :: nee                    ! net ecosys exchange (g/m2/s CO2)
  real, intent(out) :: gpp                    ! gross primary assimilation [g/m2/s C]
  real, intent(out) :: npp                    ! net primary productivity [g/m2/s C]
  real, intent(out) :: fveg                   ! Noah-MP vegetation fraction [-]
  real, intent(out) :: runsf                  ! surface runoff [mm/s]       
  real, intent(out) :: runsb                  ! subsurface runoff [mm/s]
  real, intent(out) :: ecan                   ! evaporation of intercepted water (mm/s)
  real, intent(out) :: edir                   ! soil surface evaporation rate (mm/s]
  real, intent(out) :: etran                  ! transpiration rate (mm/s)
  real, intent(out) :: rainf                  ! rainfall rate (kg/m2s)
  real, intent(out) :: snowf                  ! snowfall rate (kg/m2s) 
  real, intent(out) :: fsa                    ! total absorbed solar radiation (w/m2)
  real, intent(out) :: fira                   ! total net longwave rad (w/m2) [+ to atm]
  real, intent(out) :: apar                   ! photosyn active energy by canopy (w/m2)
  real, intent(out) :: psn                    ! total photosynthesis (umol co2/m2/s) [+]
  real, intent(out) :: sav                    ! solar rad absorbed by veg. (w/m2)
  real, intent(out) :: sag                    ! solar rad absorbed by ground (w/m2)
  real, intent(out) :: rssun                  ! sunlit leaf stomatal resistance (s/m)
  real, intent(out) :: rssha                  ! shaded leaf stomatal resistance (s/m)
  real, intent(out) :: bgap                   ! between gap fraction
  real, intent(out) :: wgap                   ! within gap fraction
  real, intent(out) :: tgb                    ! bare ground temperature [K]
  real, intent(out) :: tgv                    ! under canopy ground temperature[K]
  real, intent(out) :: chv                    ! sensible heat exchange coefficient vegetated
  real, intent(out) :: chb                    ! sensible heat exchange coefficient bare-ground
  real, intent(out) :: shg                    ! veg ground sen. heat [w/m2]   [+ to atm]
  real, intent(out) :: shc                    ! canopy sen. heat [w/m2]   [+ to atm]
  real, intent(out) :: shb                    ! bare sensible heat [w/m2]     [+ to atm]
  real, intent(out) :: evg                    ! veg ground evap. heat [w/m2]  [+ to atm]
  real, intent(out) :: evb                    ! bare soil evaporation [w/m2]  [+ to atm]
  real, intent(out) :: ghv                    ! veg ground heat flux [w/m2]  [+ to soil]
  real, intent(out) :: ghb                    ! bare ground heat flux [w/m2] [+ to soil]
  real, intent(out) :: irg                    ! veg ground net LW rad. [w/m2] [+ to atm]
  real, intent(out) :: irc                    ! canopy net LW rad. [w/m2] [+ to atm]
  real, intent(out) :: irb                    ! bare net longwave rad. [w/m2] [+ to atm]
  real, intent(out) :: tr                     ! transpiration [w/m2]  [+ to atm]
  real, intent(out) :: evc                    ! canopy evaporation heat [w/m2]  [+ to atm]
  real, intent(out) :: chleaf                 ! leaf exchange coefficient
  real, intent(out) :: chuc                   ! under canopy exchange coefficient 
  real, intent(out) :: chv2                   ! veg 2m exchange coefficient
  real, intent(out) :: chb2                   ! bare 2m exchange coefficient

! real, intent(inout) :: sfcheadrt,INFXSRT,soldrain   ! for WRF-Hydro
!--------------------------------------------------------------------------------
  ! external function
  real, external      :: month_d_401

  ! local variables 
  real                :: qsfc                 ! bulk surface specific humidity
  real                :: smstav               ! soil moisture avail. [not used], maintained for Noah consistency
  real                :: smstot               ! total soil water [mm][not used], maintained for Noah consistency 

  ! Added by Zhuo Wang and Shugong Wang
  integer             :: ids,ide,  jds,jde,  kds,kde   ! d -> domain
  integer             :: ims,ime,  jms,jme,  kms,kme   ! m -> memory
  integer             :: its,ite,  jts,jte,  kts,kte   ! t -> tile
  real                :: sr                            ! input  frozen precipitation ratio [-]
  real                :: fpice                         ! output frozen precipitation ratio [-]
! real, dimension(1:60) :: gecros1d     !  gecros crop
! real,                 :: gecros_dd ,gecros_tbem,gecros_emb ,gecros_ema, &
!                          gecros_ds1,gecros_ds2 ,gecros_ds1x,gecros_ds2x
  
  real                :: dx
  character(len=12)   :: nowdate
  integer :: k 
  integer :: ice
  integer :: ist

  ! Revised by Zhuo Wang and Shugong Wang
  real    :: xland                ! = 2 ocean; = 1 land/seaice
  real    :: xice                 ! fraction of grid that is seaice
  real    :: xice_thres           ! fraction of grid determining seaice
! integer :: isc
  integer :: croptype
  real    :: vegfra               ! vegetation fraction []
  real    :: vegmax               ! annual max vegetation fraction []

  integer :: yearlen
  real    :: julian
  real    :: cosz
  real    :: lat, lon
  real    :: q2(2)                   ! water vapor mixing ratio [kg/kg_dry]
  real    :: sfctmp(2)
  real    :: uu(2)
  real    :: vv(2)
  real    :: soldn
  real    :: lwdn
  real    :: sfcprs(2)            ! multiple layer is required 
  real    :: sfcheadrt            ! For WRF-Hydro
  real    :: dz8w3d(2) 
  real, allocatable, dimension(:) :: zsoil
  real, allocatable, dimension(:) :: zsnso
  real, allocatable, dimension(:) :: snice
  real, allocatable, dimension(:) :: snliq

  ! Added by Zhuo Wang and Shugong Wang on 10/30/2018
  real, allocatable, dimension(:) :: tsnow

  integer :: snl_idx

  ! Added by David Mocko on 11/19/2018
  logical :: Bondvillecheck
  integer :: i,local_hour
  integer :: locyr,locmo,locda,lochr,locmn,locss,locdoy
  real*8  :: loctime
  real    :: locgmt,change

  !!!! local variables for dimension match 
  !!!! Added by Zhuo Wang and Shugong Wang on 10/30/2018
  real, dimension(1,1) :: coszin 
  real, dimension(1,1) :: latin
  real, dimension(1,1) :: lonin
  integer, dimension(1,1) :: vegetypein
  integer, dimension(1,1) :: soiltypein 
  real, dimension(1,1) :: vegfrain
  real, dimension(1,1) :: vegmaxin
  real, dimension(1,1) :: tbotin 
  real, dimension(1,1) :: xlandin
  real, dimension(1,1) :: xicein
  integer, dimension(1,1) :: cropcatin
  real, dimension(1,1) :: plantingin
  real, dimension(1,1) :: harvestin
  real, dimension(1,1) :: season_gddin
  real, dimension(1,8,1) :: soilcompin
  real, dimension(1,1) :: soilcL1in
  real, dimension(1,1) :: soilcL2in
  real, dimension(1,1) :: soilcL3in
  real, dimension(1,1) :: soilcL4in
  real, dimension(1,1) :: soldnin
  real, dimension(1,1) :: lwdnin
  real, dimension(1,1) :: prcpin
  real, dimension(1,1) :: srin
  real, dimension(1,1) :: tskinout
  real, dimension(1,1) :: hfxinout 
  real, dimension(1,1) :: qfxinout
  real, dimension(1,1) :: lhinout
  real, dimension(1,1) :: grdflxinout
  real, dimension(1,1) :: smstavinout 
  real, dimension(1,1) :: smstotinout
  real, dimension(1,1) :: sfcrunoffinout 
  real, dimension(1,1) :: udrunoffinout
! real, dimension(1,1) :: albedoinout
  real, dimension(1,1) :: albedoout
  real, dimension(1,1) :: snowcinout
  real, dimension(1,nsoil,1) :: smcinout
  real, dimension(1,nsoil,1) :: sh2oinout
  real, dimension(1,nsoil,1) :: tslbinout
  real, dimension(1,1) :: sneqvinout
  real, dimension(1,1) :: snowhinout
  real, dimension(1,1) :: canwatinout
  real, dimension(1,1) :: acsnominout
  real, dimension(1,1) :: acsnowinout
  real, dimension(1,1) :: emissinout
  real, dimension(1,1) :: qsfcinout
  real, dimension(1,1) :: z0inout
  real, dimension(1,1) :: zntinout
  integer, dimension(1,1) :: isnowinout
  real, dimension(1,1) :: tvinout
  real, dimension(1,1) :: tginout
  real, dimension(1,1) :: caniceinout
  real, dimension(1,1) :: canliqinout
  real, dimension(1,1) :: eahinout 
  real, dimension(1,1) :: tahinout
  real, dimension(1,1) :: cminout
  real, dimension(1,1) :: chinout
  real, dimension(1,1) :: fwetinout
  real, dimension(1,1) :: sneqvoinout
  real, dimension(1,1) :: alboldinout
  real, dimension(1,1) :: qsnowinout
  real, dimension(1,1) :: wslakeinout
  real, dimension(1,1) :: zwtinout
  real, dimension(1,1) :: wainout
  real, dimension(1,1) :: wtinout
  real, dimension(1,-nsnow+1:0,1) :: tsnowinout
  real, dimension(1,-nsnow+1:nsoil,1) :: zsnsoinout
  real, dimension(1,-nsnow+1:0,1) :: sniceinout
  real, dimension(1,-nsnow+1:0,1) :: snliqinout
  real, dimension(1,1) :: lfmassinout 
  real, dimension(1,1) :: rtmassinout
  real, dimension(1,1) :: stmassinout
  real, dimension(1,1) :: woodinout
  real, dimension(1,1) :: stblcpinout
  real, dimension(1,1) :: fastcpinout
  real, dimension(1,1) :: laiinout
  real, dimension(1,1) :: saiinout
  real, dimension(1,1) :: taussinout
  real, dimension(1,nsoil,1) :: smoiseqinout
  real, dimension(1,1) :: smcwtdinout 
  real, dimension(1,1) :: deeprechinout
  real, dimension(1,1) :: rechinout
  real, dimension(1,1) :: graininout
  real, dimension(1,1) :: gddinout 
  integer, dimension(1,1) :: pgsinout
  real, dimension(1,60,1) :: gecros_stateinout
  real, dimension(1,1) :: t2mvout 
  real, dimension(1,1) :: t2mbout
  real, dimension(1,1) :: q2mvout
  real, dimension(1,1) :: q2mbout
  real, dimension(1,1) :: tradout
  real, dimension(1,1) :: neeout
  real, dimension(1,1) :: gppout
  real, dimension(1,1) :: nppout
  real, dimension(1,1) :: fvegout 
  real, dimension(1,1) :: runsfout
  real, dimension(1,1) :: runsbout
  real, dimension(1,1) :: ecanout 
  real, dimension(1,1) :: edirout
  real, dimension(1,1) :: etranout
  real, dimension(1,1) :: fsaout
  real, dimension(1,1) :: firaout
  real, dimension(1,1) :: aparout
  real, dimension(1,1) :: psnout
  real, dimension(1,1) :: savout
  real, dimension(1,1) :: sagout
  real, dimension(1,1) :: rssunout
  real, dimension(1,1) :: rsshaout
  real, dimension(1,1) :: bgapout 
  real, dimension(1,1) :: wgapout
  real, dimension(1,1) :: tgvout
  real, dimension(1,1) :: tgbout
  real, dimension(1,1) :: chvout
  real, dimension(1,1) :: chbout
  real, dimension(1,1) :: shgout
  real, dimension(1,1) :: shcout
  real, dimension(1,1) :: shbout
  real, dimension(1,1) :: evgout
  real, dimension(1,1) :: evbout
  real, dimension(1,1) :: ghvout
  real, dimension(1,1) :: ghbout
  real, dimension(1,1) :: irgout
  real, dimension(1,1) :: ircout
  real, dimension(1,1) :: irbout
  real, dimension(1,1) :: trout
  real, dimension(1,1) :: evcout
  real, dimension(1,1) :: chleafout
  real, dimension(1,1) :: chucout
  real, dimension(1,1) :: chv2out
  real, dimension(1,1) :: chb2out
  real, dimension(1,1) :: rsout

   ids = 1
   ide = 1
   jds = 1
   jde = 1
   kds = 1
   kde = 1
   ims = 1
   ime = 1
   jms = 1
   jme = 1
   kms = 1
   kme = 2
   its = 1
   ite = 1
   jts = 1
   jte = 1
   kts = 1
   kte = 1

   xland = 1
   xice  = 0
   ! Fraction of grid determining seaice (from WRF and HRLDAS)
   xice_thres = 0.5

   sfcheadrt = 0.0
! SR = frozen precip fraction.  For offline, it is set to zero.
! If running coupled to WRF, then SR is set by the WRF model.
   if (snf_opt.ne.4) then
      sr = 0.0
   else
      write(LIS_logunit,*) "[ERR] SR should be set by the WRF model."
      write(LIS_logunit,*) "[ERR] Code needs fixing.  Stopping LIS."
      call LIS_endrun
   endif

   ! dx is horizontal grid spacing. dx is not used past this point,
   ! but is used earlier when run_opt=5 (new groundwater scheme).
   dx = 0.0

  !!!!! print all the options not supported in offline mode 
  if (sfc_opt == 3) then
     stop "(opt_sfc == 3) and (opt_sfc == 4) are not for offline use"
  endif

  ! set ZSOIL 
  allocate(zsoil(nsoil))
  do k = 1,nsoil
     zsoil(k) = sldpth(k)
  enddo

  allocate(zsnso(-nsnow+1:nsoil))
  allocate(snice(-nsnow+1:0))
  allocate(snliq(-nsnow+1:0))
  allocate(tsnow(-nsnow+1:0))   ! Added by Zhuo Wang and Shugong Wang on 10/30/2018

  smstav = 0.0       ! Not used
  smstot = 0.0       ! Not used

  ! state variables 
  zsnso(-nsnow+1:nsoil) = zss(1:nsnow+nsoil) 
  snice(-nsnow+1:0)     = snowice(1:nsnow)
  snliq(-nsnow+1:0)     = snowliq(1:nsnow) 
  tsnow(-nsnow+1:0)     = tsno(1:nsnow) ! Added by Zhuo Wang and Shugong Wang on 10/30/2018

  ! cosz, yearlen, and julian are calculated in subroutine calc_declin_401.
  ! Be careful here!!!, LIS uses GMT; the date for calc_declin_401 should
  ! be local time.  Longitude is need to convert GMT into local time!!!.
  ! This should only be done when using Bondville single_point forcing,
  ! as the forcing file uses local time instead of GMT.
  Bondvillecheck = .false.
  do i=1,LIS_rc%nmetforc
     if (trim(LIS_rc%metforc(i)).eq."Bondville") Bondvillecheck = .true.
  enddo

! For a true benchmark against the HRLDAS Noah-MP-4.0.1 testcase
! from NCAR, set "change = -21600.0".  This line allows the code
! to _incorrectly_ run in the same way as the HRLDAS testcase,
! which runs on local time instead of on UTC time.  The changes
! is 6 hours * 3600.0 seconds per hour, to go from GMT to the
! Bondville location in Illinois, which uses Central time.
  change = 0.0
  if (Bondvillecheck) change = -21600.0
  locyr = LIS_rc%yr
  locmo = LIS_rc%mo
  locda = LIS_rc%da
  lochr = LIS_rc%hr
  locmn = LIS_rc%mn
  locss = LIS_rc%ss

  call LIS_date2time(loctime,locdoy,locgmt,locyr,locmo,locda,lochr,locmn,locss)
  call LIS_tick(loctime,locdoy,locgmt,locyr,locmo,locda,lochr,locmn,locss,change)

  write(nowdate,'(I4.4,4I2.2)') locyr, locmo, locda, lochr, locmn
  call calc_declin_401(nowdate(1:4)//"-"//nowdate(5:6)//"-"//nowdate(7:8)//"_"//nowdate(9:10)//":"//nowdate(11:12)//":00", &
      latitude, longitude, cosz, yearlen, julian)

  if ( dveg_opt == 1 ) then
    ! with dveg_opt==1, shdfac is fed directly to fveg
    vegfra = month_d_401(shdfac_monthly, nowdate)
  else
    ! with dveg_opt==2, fveg is computed from lai and sai, and shdfac is unused
    vegfra = -1.E36   ! To match with HRLDAS initialization
  endif
  vegmax = maxval(shdfac_monthly)

  ! assign forcing variables 
  dz8w3d(:) = dz8w 
  sfctmp(:) = tair 
  sfcprs(:) = psurf 
  q2(:)     = qair/(1.0-qair)   ! Convert specific humidity to water vapor mixing ratio
  soldn  = swdown
  lwdn   = lwdown
  qsfc = qair 
  uu(:) = wind_e
  vv(:) = wind_n 

  !!!! set up local variables for dimension match 
  !!!! Added by Zhuo Wang and Shugong Wang on 10/30/2018 
  coszin(1,1)   = cosz
  latin(1,1)    = latitude
  lonin(1,1)    = longitude
  vegetypein(1,1) = vegetype
  soiltypein(1,1) = soiltype
  vegfrain(1,1)   = vegfra
  vegmaxin(1,1)   = vegmax
  tbotin(1,1)     = tbot
  xlandin(1,1)    = xland
  xicein(1,1)     = xice
  cropcatin(1,1)  = cropcat
  plantingin(1,1) = planting
  harvestin(1,1)  = harvest
  season_gddin(1,1) = season_gdd
  soilcompin(1,:,1) = soilcomp(:)
  soilcL1in(1,1)  = soilcL1
  soilcL2in(1,1)  = soilcL2
  soilcL3in(1,1)  = soilcL3
  soilcL4in(1,1)  = soilcL4
  soldnin(1,1)    = soldn
  lwdnin(1,1)     = lwdn
  prcpin(1,1)     = prcp
  srin(1,1)       = sr
  tskinout(1,1)   = tsk
  hfxinout(1,1)   = hfx
  qfxinout(1,1)   = qfx
  lhinout(1,1)    = lh
  grdflxinout(1,1) = grdflx
  smstavinout(1,1) = smstav
  smstotinout(1,1) = smstot
  sfcrunoffinout(1,1) = sfcrunoff
  udrunoffinout(1,1)  = udrunoff
! albedoinout(1,1) = albedo
  snowcinout(1,1)  = snowc
  smcinout(1,:,1)  = smc(:)
  sh2oinout(1,:,1) = sh2o(:)
  tslbinout(1,:,1) = tslb(:)
  sneqvinout(1,1)  = sneqv  
  snowhinout(1,1)  = snowh
  canwatinout(1,1) = canwat
  acsnominout(1,1) = acsnom
  acsnowinout(1,1) = acsnow
  emissinout(1,1)  = emiss
  qsfcinout(1,1)   = qsfc
! If coupled to WRF, set these variables to realistic values,
! and then pass back to WRF after the call to noahmplsm_401.
  z0inout(1,1)     = 0.0
  zntinout(1,1)    = 0.0
! z0 and znt should be passed to WRF, if coupled. - dmm
  isnowinout(1,1)  = isnow
  tvinout(1,1)     = tv
  tginout(1,1)     = tg
  caniceinout(1,1) = canice
  canliqinout(1,1) = canliq
  eahinout(1,1)    = eah
  tahinout(1,1)    = tah
  cminout(1,1)     = cm
  chinout(1,1)     = ch
  fwetinout(1,1)   = fwet
  sneqvoinout(1,1) = sneqvo
  alboldinout(1,1) = albold
  qsnowinout(1,1)  = qsnow
  wslakeinout(1,1)  = wslake
  zwtinout(1,1)    = zwt
  wainout(1,1)     = wa
  wtinout(1,1)     = wt
  tsnowinout(1,:,1) = tsnow(:)
  zsnsoinout(1,:,1)  = zsnso(:)
  sniceinout(1,:,1) = snice(:)
  snliqinout(1,:,1) = snliq(:)
  lfmassinout(1,1) = lfmass
  rtmassinout(1,1) = rtmass
  stmassinout(1,1) = stmass
  woodinout(1,1)   = wood
  stblcpinout(1,1) = stblcp
  fastcpinout(1,1) = fastcp
  laiinout(1,1)    = lai
  saiinout(1,1)    = sai
  taussinout(1,1)  = tauss
  smoiseqinout(1,:,1) = smoiseq(:)
  smcwtdinout(1,1) = smcwtd
  deeprechinout(1,1)  = deeprech
  rechinout(1,1)      = rech
  graininout(1,1)     = grain
  gddinout(1,1)       = gdd
  pgsinout(1,1)       = pgs
  gecros_stateinout(1,:,1) = gecros_state(:)
  t2mvout(1,1) = t2mv
  t2mbout(1,1) = t2mb
  q2mvout(1,1) = q2mv
  q2mbout(1,1) = q2mb
  tradout(1,1) = trad
  neeout(1,1)  = nee
  gppout(1,1)  = gpp
  nppout(1,1)  = npp
  fvegout(1,1) = fveg
  runsfout(1,1) = runsf
  runsbout(1,1) = runsb
  ecanout(1,1)  = ecan
  edirout(1,1)  = edir
  etranout(1,1) = etran
  fsaout(1,1)   = fsa
  firaout(1,1)  = fira
  aparout(1,1)  = apar
  psnout(1,1)   = psn
  savout(1,1)   = sav
  sagout(1,1)   = sag
  rssunout(1,1) = rssun
  rsshaout(1,1) = rssha
  bgapout(1,1)  = bgap
  wgapout(1,1)  = wgap
  tgvout(1,1)   = tgv
  tgbout(1,1)   = tgb
  chvout(1,1)   = chv
  chbout(1,1)   = chb
  shgout(1,1)   = shg
  shcout(1,1)   = shc
  shbout(1,1)   = shb
  evgout(1,1)   = evg
  evbout(1,1)   = evb
  ghvout(1,1)   = ghv
  ghbout(1,1)   = ghb
  irgout(1,1)   = irg
  ircout(1,1)   = irc
  irbout(1,1)   = irb
  trout(1,1)    = tr
  evcout(1,1)   = evc
  chleafout(1,1) = chleaf
  chucout(1,1)  = chuc
  chv2out(1,1)  = chv2
  chb2out(1,1)  = chb2
  rsout(1,1)    = rs

! Code from module_NoahMP_hrldas_driver.F.  Initial guess only.
  if ((trim(LIS_rc%startcode).eq."coldstart").and.(itimestep.eq.1)) then
      eahinout(1,1) = sfcprs(1) * (q2(1)/(0.622+q2(1)))
      tahinout(1,1) = sfctmp(1)
      cminout(1,1)  = 0.1
      chinout(1,1)  = 0.1
  endif

  call noahmplsm_401  (LIS_rc%udef,  & ! in : LIS undefined value (David Mocko)
                       itimestep,yearlen , julian  , coszin    , latin   , lonin  , & ! in : time/space-related
                       dz8w3d(1), dt      , zsoil   , nsoil   , dx      ,           & ! in : model configuration 
                       vegetypein, soiltypein, vegfrain, vegmaxin, tbotin  ,        & ! in : Vegetation/Soil characteristics
                       xlandin , xicein  , xice_thres,                             & ! in : Vegetation/Soil characteristics
                       cropcatin , plantingin, harvestin ,season_gddin,                    &
                       dveg_opt, crs_opt , btr_opt ,run_opt  , sfc_opt , frz_opt,  & ! in : user options
                       inf_opt , rad_opt , alb_opt ,snf_opt  , tbot_opt, stc_opt,  & ! in : user options
                       gla_opt , rsf_opt , soil_opt,pedo_opt , crop_opt,           & ! in : user options
                       iz0tlnd , urban_opt,                                        & ! in : user options
                       soilcompin, soilcL1in, soilcL2in, soilcL3in, soilcL4in,               & ! in : user options
                       sfctmp(1)  , q2(1)  , uu(1)    , vv(1) , soldnin , lwdnin  , & ! in : forcing 
                       sfcprs(1)  , prcpin  , srin      ,                             & ! in : forcing

                       tskinout, hfxinout, qfxinout, lhinout , grdflxinout  , smstavinout  , & ! in/out LSM eqv 

                       smstotinout  ,sfcrunoffinout, udrunoffinout, albedoout  , qsnbot  , subsnow, & ! in/out LSM eqv 
                       snowcinout   , smcinout     ,           pah,                                 & ! in/out LSM eqv
                       sh2oinout    , tslbinout    , sneqvinout   , snowhinout   , canwatinout  , acsnominout  , & ! in/out LSM eqv
                       acsnowinout  , emissinout   , qsfcinout    , z0inout      , zntinout     ,           & ! in/out LSM eqv

                       isnowinout   , tvinout      , tginout      , caniceinout  , canliqinout  , eahinout     , & ! in/out Noah MP only
                       tahinout     , cminout      , chinout      , fwetinout    , sneqvoinout  , alboldinout  , & ! in/out Noah MP only
                       qsnowinout   , wslakeinout  , zwtinout     , wainout      , wtinout      , tsnowinout    , & ! in/out Noah MP only
                       zsnsoinout     , sniceinout , snliqinout , lfmassinout  , rtmassinout  , stmassinout  , & ! in/out Noah MP only
                       woodinout    , stblcpinout  , fastcpinout  , laiinout     , saiinout     , taussinout   , & ! in/out Noah MP only
                       smoiseqinout , smcwtdinout  ,deeprechinout , rechinout    , graininout   , gddinout     , & ! in/out Noah MP only 
                       pgsinout     ,                                                   & ! in/out Noah MP only
                       gecros_stateinout,                                               & ! in/out gecros model

                       t2mvout , t2mbout , q2mvout , q2mbout ,                     & ! out Noah MP only
                       tradout , neeout  , gppout  , nppout  , fvegout , runsfout, & ! out Noah MP only
                       runsbout, ecanout , edirout , etranout, fsaout  , firaout , & ! out Noah MP only
                       aparout , psnout  , savout  , sagout  , rssunout, rsshaout, & ! out Noah MP only
                       bgapout , wgapout , tgvout  , tgbout  , chvout  , chbout  , & ! out Noah MP only
                       shgout  , shcout  , shbout  , evgout  , evbout  , ghvout  , & ! out Noah MP only
                       ghbout  , irgout  , ircout  , irbout  , trout   , evcout  , & ! out Noah MP only
                       chleafout  , chucout , chv2out , chb2out , rsout , fpice  , & ! out Noah MP only
#ifdef WRF_HYDRO
                       sfcheadrt, INFXSRT, soldrain,                               &
#endif
                       ids,ide,  jds,jde,  kds,kde,                                &
                       ims,ime,  jms,jme,  kms,kme,                                &
                       its,ite,  jts,jte,  kts,kte)

  ! Added by Zhuo Wang and Shugong on 10/30/2018
  tsk = tskinout(1,1)
  hfx = hfxinout(1,1)
  qfx = qfxinout(1,1)
  lh  = lhinout(1,1)
  grdflx = grdflxinout(1,1)
  smstav = smstavinout(1,1)
  smstot = smstotinout(1,1)
  sfcrunoff = sfcrunoffinout(1,1)
  udrunoff  = udrunoffinout(1,1)
! albedo    = albedoinout(1,1)
  albedo    = albedoout(1,1)
  snowc  = snowcinout(1,1)
  smc(:) = smcinout(1,:,1) 
  sh2o(:) = sh2oinout(1,:,1)
  tslb(:) = tslbinout(1,:,1)
  sneqv = sneqvinout(1,1)
  snowh = snowhinout(1,1)
  canwat = canwatinout(1,1)
  acsnom = acsnominout(1,1)
  acsnow = acsnowinout(1,1)
  emiss = emissinout(1,1)
  qsfc = qsfcinout(1,1)
  isnow = isnowinout(1,1)
  tv = tvinout(1,1)
  tg = tginout(1,1)
  canice = caniceinout(1,1)
  canliq = canliqinout(1,1)
  eah = eahinout(1,1)
  tah = tahinout(1,1)
  cm = cminout(1,1)
  ch = chinout(1,1)
  fwet = fwetinout(1,1)
  sneqvo = sneqvoinout(1,1)
  albold = alboldinout(1,1)
  qsnow = qsnowinout(1,1)
  wslake = wslakeinout(1,1)
  zwt = zwtinout(1,1)
  wa = wainout(1,1)
  wt = wtinout(1,1)

  ! Modified by Zhuo Wang on 12/31/2018 
  tsnow(-nsnow+1:0) = tsnowinout(1,-nsnow+1:0,1)
  zsnso(-nsnow+1:nsoil) = zsnsoinout(1,-nsnow+1:nsoil,1)
  snice(-nsnow+1:0) = sniceinout(1,-nsnow+1:0,1)
  snliq(-nsnow+1:0) = snliqinout(1,-nsnow+1:0,1)
  tsno(1:nsnow)       = tsnow(-nsnow+1:0) 
  zss(1:nsnow+nsoil)  = zsnso(-nsnow+1:nsoil) 
  snowice(1:nsnow)    = snice(-nsnow+1:0)     
  snowliq(1:nsnow)    = snliq(-nsnow+1:0)   

  lfmass = lfmassinout(1,1)
  rtmass = rtmassinout(1,1)
  stmass = stmassinout(1,1)
  wood = woodinout(1,1)
  stblcp = stblcpinout(1,1)
  fastcp = fastcpinout(1,1)
  lai = laiinout(1,1)
  sai = saiinout(1,1)
  tauss = taussinout(1,1)
  smoiseq(:) = smoiseqinout(1,:,1)
  smcwtd = smcwtdinout(1,1)
  deeprech = deeprechinout(1,1)
  rech = rechinout(1,1)
  grain = graininout(1,1)
  gdd = gddinout(1,1)
  pgs = pgsinout(1,1)
  gecros_state(:) = gecros_stateinout(1,:,1)
  t2mv = t2mvout(1,1)
  t2mb = t2mbout(1,1)
  q2mv = q2mvout(1,1)
  q2mb = q2mbout(1,1)
  trad = tradout(1,1)
  nee = neeout(1,1)
  gpp = gppout(1,1)
  npp = nppout(1,1)
  fveg = fvegout(1,1)
  runsf = runsfout(1,1)
  runsb = runsbout(1,1)
  ecan = ecanout(1,1)
  edir = edirout(1,1)
  etran = etranout(1,1)
  fsa = fsaout(1,1)
  fira = firaout(1,1)
  apar = aparout(1,1)
  psn = psnout(1,1)
  sav = savout(1,1)
  sag = sagout(1,1)
  rssun = rssunout(1,1)
  rssha = rsshaout(1,1)
  bgap = bgapout(1,1)
  wgap = wgapout(1,1)
  tgv = tgvout(1,1)
  tgb = tgbout(1,1)
  chv = chvout(1,1)
  chb = chbout(1,1)
  shg = shgout(1,1)
  shc = shcout(1,1)
  shb = shbout(1,1)
  evg = evgout(1,1)
  evb = evbout(1,1)
  ghv = ghvout(1,1)
  ghb = ghbout(1,1)
  irg = irgout(1,1)
  irc = ircout(1,1)
  irb = irbout(1,1)
  tr = trout(1,1)
  evc = evcout(1,1)
  chleaf = chleafout(1,1)
  chuc = chucout(1,1)
  chv2 = chv2out(1,1)
  chb2 = chb2out(1,1)
  rs = rsout(1,1)

  rainf = prcp * (1.0 - fpice)/dt  ! added by Shugong for LIS output 
  snowf = prcp * fpice/dt          ! added by Shugong for LIS output 

  deallocate(zsoil)
  deallocate(zsnso)

  deallocate(snice)
  deallocate(snliq)
  deallocate(tsnow)   ! Added by Zhuo Wang and Shugong Wang on 10/30/2018

end subroutine noahmp_driver_401

real function month_d_401(a12, nowdate) result (nowval)
  !
  ! Given a set of 12 values, taken to be valid on the fifteenth of each month (Jan through Dec)
  ! and a date in the form <YYYYMMDD[HHmmss]> ....
  ! 
  ! Return a value valid for the day given in <nowdate>, as an interpolation from the 12
  ! monthly values.
  !
  use kwm_date_utilities_401
  implicit none
  real, dimension(12), intent(in) :: a12 ! 12 monthly values, taken to be valid on the 15th of
  !                                      ! the month
  character(len=12), intent(in) :: nowdate ! Date, in the form <YYYYMMDD[HHmmss]>
  integer :: nowy, nowm, nowd
  integer :: prevm, postm
  real    :: factor
  integer, dimension(12) :: ndays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  !
  ! Handle leap year by setting the number of days in February for the year in question.
  !
  read(nowdate(1:8),'(I4,I2,I2)') nowy, nowm, nowd
  ndays(2) = nfeb(nowy)

  !
  ! Do interpolation between the fifteenth of two successive months.
  !
  if (nowd == 15) then
     nowval = a12(nowm)
     return
  else if (nowd < 15) then
     postm = nowm
     prevm = nowm - 1
     if (prevm == 0) prevm = 12
     factor = real(ndays(prevm)-15+nowd)/real(ndays(prevm))
  else if (nowd > 15) then
     prevm = nowm
     postm = nowm + 1
     if (postm == 13) postm = 1
     factor = real(nowd-15)/real(ndays(prevm))
  endif

  nowval = a12(prevm)*(1.0-factor) + a12(postm)*factor

end function month_d_401 

SUBROUTINE calc_declin_401 ( nowdate, latitude, longitude, cosz, yearlen, julian)
  use kwm_date_utilities_401 
!---------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
  character(len=19), intent(in)  :: nowdate    ! YYYY-MM-DD_HH:mm:ss
  real,              intent(in)  :: latitude
  real,              intent(in)  :: longitude
  real,              intent(out) :: cosz
  integer,           intent(out) :: yearlen
  real,              intent(out) :: JULIAN

  REAL                           :: hrang
  real                           :: DECLIN
  real                           :: GMT
  real                           :: tloctim
  REAL                           :: OBECL
  REAL                           :: SINOB
  REAL                           :: SXLONG
  REAL                           :: ARG
  integer                        :: iyear
  integer                        :: iday
  integer                        :: ihour
  integer                        :: iminute
  integer                        :: isecond

  REAL, PARAMETER :: DEGRAD = 3.14159265/180.
  REAL, PARAMETER :: DPD    = 360./365.

  !
  ! Determine the number of days in the year
  !

  read(nowdate(1:4), '(I4)') iyear
  yearlen = 365
  if (mod(iyear,4) == 0) then
     yearlen = 366
     if (mod(iyear,100) == 0) then
        yearlen = 365
        if (mod(iyear,400) == 0) then
           yearlen = 366
           if (mod(iyear,3600) == 0) then
              yearlen = 365
           endif
        endif
     endif
  endif

  !
  ! Determine the Julian time (floating-point day of year).
  !

  call geth_idts(nowdate(1:10), nowdate(1:4)//"-01-01", iday)
  read(nowdate(12:13), *) ihour
  read(nowdate(15:16), *) iminute
  read(nowdate(18:19), *) isecond
  GMT = REAL(IHOUR) + IMINUTE/60.0 + ISECOND/3600.0
  JULIAN = REAL(IDAY) + GMT/24.

! for short wave radiation

  DECLIN=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.

  OBECL=23.5*DEGRAD
  SINOB=SIN(OBECL)

!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

  IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
  IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
  ARG=SINOB*SIN(SXLONG)
  DECLIN=ASIN(ARG)

  TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! Local time in hours
  tloctim = AMOD(tloctim+24.0, 24.0)
  HRANG=15.*(TLOCTIM-12.)*DEGRAD
  COSZ=SIN(LATITUDE*DEGRAD)*SIN(DECLIN)+COS(LATITUDE*DEGRAD)*COS(DECLIN)*COS(HRANG)
  COSZ=MIN(COSZ,1.0);   !Added by kwH 3/1/16 to address floating point roundoff errors 
  COSZ=MAX(COSZ,-1.0);  !

!KWM   write(wrf_err_message,10)DECDEG/DEGRAD
!KWM10 FORMAT(1X,'*** SOLAR DECLINATION ANGLE = ',F6.2,' DEGREES.',' ***')
!KWM   CALL wrf_debug (50, wrf_err_message)

END SUBROUTINE calc_declin_401 

! Subroutine SNOW_INIT grabbed from NOAH-MP-WRF
SUBROUTINE SNOW_INIT_401 ( jts, jtf, its, itf, ims, ime, jms, jme, NSNOW, NSOIL, ZSOIL,  &
     SWE, tgxy, SNODEP, ZSNSOXY, TSNOXY, SNICEXY, SNLIQXY, ISNOWXY)

! ------------------------------------------------------------------------------------------
  IMPLICIT NONE
! ------------------------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: jts,jtf,its,itf,ims,ime, jms,jme,NSNOW,NSOIL
  REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: SWE 
  REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: SNODEP
  REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: tgxy
  REAL,    INTENT(IN), DIMENSION(1:NSOIL) :: ZSOIL

  INTEGER, INTENT(OUT), DIMENSION(ims:ime, jms:jme) :: ISNOWXY
  REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:NSOIL,jms:jme) :: ZSNSOXY
  REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: TSNOXY
  REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: SNICEXY
  REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: SNLIQXY

!local
  INTEGER :: I,J,IZ
  REAL,                 DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: DZSNOXY
  REAL,                 DIMENSION(ims:ime, -NSNOW+1:NSOIL,jms:jme) :: DZSNSOXY
! ------------------------------------------------------------------------------------------


  DO J = jts,jtf
     DO I = its,itf
        IF (SNODEP(I,J) < 0.025) THEN
           ISNOWXY(I,J) = 0
           DZSNOXY(I,-NSNOW+1:0,J) = 0.
        ELSE
           IF ((SNODEP(I,J) >= 0.025) .AND. (SNODEP(I,J) <= 0.05)) THEN
              ISNOWXY(I,J)    = -1
              DZSNOXY(I,0,J)  = SNODEP(I,J)
           ELSE IF ((SNODEP(I,J) > 0.05) .AND. (SNODEP(I,J) <= 0.10)) THEN
              ISNOWXY(I,J)    = -2
              DZSNOXY(I,-1,J) = SNODEP(I,J)/2.
              DZSNOXY(I, 0,J) = SNODEP(I,J)/2.
           ELSE IF ((SNODEP(I,J) > 0.10) .AND. (SNODEP(I,J) <= 0.25)) THEN
              ISNOWXY(I,J)    = -2
              DZSNOXY(I,-1,J) = 0.05
              DZSNOXY(I, 0,J) = SNODEP(I,J) - DZSNOXY(I,-1,J)
           ELSE IF ((SNODEP(I,J) > 0.25) .AND. (SNODEP(I,J) <= 0.35)) THEN
              ISNOWXY(I,J)    = -3
              DZSNOXY(I,-2,J) = 0.05
              DZSNOXY(I,-1,J) = 0.5*(SNODEP(I,J)-DZSNOXY(I,-2,J))
              DZSNOXY(I, 0,J) = 0.5*(SNODEP(I,J)-DZSNOXY(I,-2,J))
           ELSE IF (SNODEP(I,J) > 0.35) THEN
              ISNOWXY(I,J)     = -3
              DZSNOXY(I,-2,J) = 0.05
              DZSNOXY(I,-1,J) = 0.10
              DZSNOXY(I, 0,J) = SNODEP(I,J) - DZSNOXY(I,-1,J) - DZSNOXY(I,-2,J)
           END IF
        END IF
     ENDDO
  ENDDO

  DO J = jts,jtf
     DO I = its,itf
        TSNOXY( I,-NSNOW+1:0,J) = 0.
        SNICEXY(I,-NSNOW+1:0,J) = 0.
        SNLIQXY(I,-NSNOW+1:0,J) = 0.
        DO IZ = ISNOWXY(I,J)+1, 0
           TSNOXY(I,IZ,J)  = tgxy(I,J)  ! [k]
           SNLIQXY(I,IZ,J) = 0.00
           SNICEXY(I,IZ,J) = 1.00 * DZSNOXY(I,IZ,J) * (SWE(I,J)/SNODEP(I,J))  ! [kg/m3]
        END DO

        DO IZ = ISNOWXY(I,J)+1, 0
           DZSNSOXY(I,IZ,J) = -DZSNOXY(I,IZ,J)
        END DO

        DZSNSOXY(I,1,J) = ZSOIL(1)
        DO IZ = 2,NSOIL
           DZSNSOXY(I,IZ,J) = (ZSOIL(IZ) - ZSOIL(IZ-1))
        END DO

        ZSNSOXY(I,ISNOWXY(I,J)+1,J) = DZSNSOXY(I,ISNOWXY(I,J)+1,J)
        DO IZ = ISNOWXY(I,J)+2 ,NSOIL
           ZSNSOXY(I,IZ,J) = ZSNSOXY(I,IZ-1,J) + DZSNSOXY(I,IZ,J)
        ENDDO

     END DO
  END DO

END SUBROUTINE SNOW_INIT_401
