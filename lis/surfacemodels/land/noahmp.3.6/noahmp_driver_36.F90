!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#define LIS_NoahMP_TEST 0
! !INTERFACE
subroutine noahmp_driver_36(iloc, jloc, &
     landuse_tbl_name, soil_tbl_name,  gen_tbl_name,             & 
                            noahmp_tbl_name, landuse_scheme_name, soil_scheme_name,     &
                            dveg_opt, crs_opt, btr_opt, run_opt, sfc_opt, frz_opt,      &
                            inf_opt, rad_opt, alb_opt , snf_opt, tbot_opt, stc_opt,     &
                            nslcats, nlucats, nslpcats,                                 &
                            latitude, longitude,                                        &
                            year    , month   , day     , hour    , minute   ,          &
                            dt      , upd_alb_flag,albd_in,  albi_in, nsoil   , sldpth  , nsnow,  & ! in : model configuration 
                            shdfac_monthly,                                             &   
                            vegetype  , soiltype, slopetype, urban_vegetype,            &
                            ice_flag, st_flag, sc_idx, iz0tlnd , smceq,                 & ! in : user options
                            tair    , psurf   , wind_e  , wind_n  , qair    ,           & ! in : forcing
                            swdown  , lwdown  , prcp    ,                               &
                            tbot    , pblh    , zlvl    ,                               & ! in : forcing
                            p_csoil , p_bexp  , p_dksat , p_dwsat , p_psisat,           & ! SY: in : calibratable parameters for enabling OPTUE
                            p_quartz, p_smcmax, p_smcref, p_smcwlt,                     & ! SY: in : calibratable parameters for enabling OPTUE
                            p_czil  , p_frzk  , p_refdk , p_refkdt, p_slope ,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_topt  , p_rgl   , p_rsmax , p_rsmin , p_hs, p_nroot,     & ! SY: in : calibratable parameters for enabling OPTUE
                            p_CH2OP , p_DLEAF , p_Z0MVT , p_HVT   , p_HVB   ,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_RC    , p_RHOL1 , p_RHOL2 , p_RHOS1 , p_RHOS2 ,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_TAUL1 , p_TAUL2 , p_TAUS1 , p_TAUS2 , p_XL    ,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_CWPVT , p_C3PSN , p_KC25  , p_AKC   , p_KO25  ,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_AKO   , p_AVCMX , p_AQE   , p_LTOVRC, p_DILEFC,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_DILEFW, p_RMF25 , p_SLA   , p_FRAGR , p_TMIN  ,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_VCMX25, p_TDLEF , p_BP    , p_MP    , p_QE25  ,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_RMS25 , p_RMR25 , p_ARM   , p_FOLNMX, p_WDPOOL,          & ! SY: in : calibratable parameters for enabling OPTUE
                            p_WRRAT , p_MRP   ,                                        & ! SY: in : calibratable parameters for enabling OPTUE
                            albold  , sneqvo  ,                                         & ! in/out : 
                            sstc    , sh2o    , smc     , tah     , eah     , fwet    , & ! in/out : 
                            canliq  , canice  , tv      , tg      , qsnow   ,           & ! in/out : 
                            isnow   , zss     , snowh   , sneqv   , snowice , snowliq , & ! in/out : 
                            zwt     , wa      , wt      , wslake  , lfmass  , rtmass  , & ! in/out : 
                            stmass  , wood    , stblcp  , fastcp  , lai     , sai     , & ! in/out : 
                            cm      , ch      , tauss   ,                               & ! in/out : 
                            smcwtd  ,deeprech , rech    ,                               & ! in/out :
                            fsa     , fsr     , fira    , fsh     , ssoil   , fcev    , & ! out : 
                            fgev    , fctr    , ecan    , etran   , edir    , trad    , & ! out :
                            subsnow ,                                                   & ! out :
                            tgb     , tgv     , t2mv    , t2mb    , q2v     , q2b     , & ! out :
                            runsrf  , runsub  , apar    , psn     , sav     , sag     , & ! out :
                            fsno    , nee     , gpp     , npp     , fveg    , albedo  , & ! out :
                            qsnbot  , ponding , ponding1, ponding2, rssun   , rssha   , & ! out :
                            bgap    , wgap    , chv     , chb     , emissi  ,           & ! out :
                            shg     , shc     , shb     , evg     , evb     , ghv     , & ! out :
                            ghb     , irg     , irc     , irb     , tr      , evc     , & ! out :
                            chleaf  , chuc    , chv2    , chb2    , fpice   , &
                            !ag (12Sep2019)
                            rivsto, fldsto, fldfrc, &
                            sfcheadrt)  ! out 
  
  ! use LIS_FORC_AttributesMod 
  use module_sf_noahlsm_36, only: slcats, lucats, slpcats
  use noahmp_globals_36, only: dveg, opt_crs , opt_btr , opt_run ,      &
                               opt_sfc , opt_frz, opt_inf , opt_rad ,   &
                               opt_alb , opt_snf , opt_tbot, opt_stc,   & ! SY
                               CSOIL_DATA, BB, SATDK, SATDW, & ! SY
                               SATPSI, QTZ, MAXSMC, REFSMC, WLTSMC, & ! SY
                               CZIL_DATA, FRZK_DATA, REFDK_DATA, REFKDT_DATA, SLOPE_DATA, & ! SY
                               TOPT_DATA, RGLTBL, RSMAX_DATA, RSTBL, HSTBL, NROTBL, & ! SY
                               CH2OP, DLEAF, Z0MVT, HVT, HVB, RC, RHOL, RHOS, TAUL, TAUS, & ! SY
                               XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, & ! SY
                               LTOVRC,  DILEFC,  DILEFW,  RMF25,  SLA,  FRAGR,  TMIN, & ! SY
                               VCMX25,  TDLEF,  BP, MP, QE25, RMS25, RMR25, ARM, & ! SY
                               FOLNMX, WDPOOL, WRRAT, MRP ! SY    
  use noahmp_routines_36, only: noahmp_sflx_36, redprm 
  implicit none
  character(len=256), intent(in) :: landuse_tbl_name      ! Noah model landuse parameter table
  character(len=256), intent(in) :: soil_tbl_name         ! Noah model soil parameter table
  character(len=256), intent(in) :: gen_tbl_name          ! Noah model general parameter table
  character(len=256), intent(in) :: noahmp_tbl_name       ! NoahMP parameter table
  character(len=256), intent(in) :: landuse_scheme_name   ! Landuse classficiation scheme
  character(len=256), intent(in) :: soil_scheme_name      ! Soil classification scheme  
  integer, intent(in) :: dveg_opt             ! vegetation model ( 1->prescribed [table LAI, shdfac=FVEG]; 2->dynamic; 3->table LAI, calculate FVEG 4->table LAI, shdfac=maximum)    
  integer, intent(in) :: crs_opt              ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  integer, intent(in) :: btr_opt              ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  integer, intent(in) :: run_opt              ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  integer, intent(in) :: sfc_opt              ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  integer, intent(in) :: frz_opt              ! supercooled liquid water (1-> NY06; 2->Koren99)
  integer, intent(in) :: inf_opt              ! frozen soil permeability (1-> NY06; 2->Koren99)
  integer, intent(in) :: rad_opt              ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1--Fveg)
  integer, intent(in) :: alb_opt              ! snow surface albedo (1->BATS; 2->CLASS)
  integer, intent(in) :: snf_opt              ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
  integer, intent(in) :: tbot_opt             ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
  integer, intent(in) :: stc_opt              ! snow/soil temperature time scheme (1->semi-implicit; 2->fully implicit)                                                             
  integer, intent(in) :: nslcats              ! the number of total soil types in parameter tabel [-]
  integer, intent(in) :: nlucats              ! the number of total land cover types in parameter tabel [-]
  integer, intent(in) :: nslpcats             ! the number of total slope category for Noah baseflow [-] 
  real,    intent(in) :: latitude             ! latitude in decimal degree [-]
  real,    intent(in) :: longitude            ! longitude in decimal degree [-]
  integer, intent(in) :: year                 ! year of the current time step [-]
  integer, intent(in) :: month                ! month of the current time step [-] 
  integer, intent(in) :: day                  ! day of the current time step [-]
  integer, intent(in) :: hour                 ! hour of the current time step [-] 
  integer, intent(in) :: minute               ! minute of the current time step [-]
  real,    intent(in) :: dt                   ! time step in seconds [s]  
  logical             :: upd_alb_flag
  real                :: albd_in(2)
  real                :: albi_in(2)
  integer, intent(in) :: nsoil                ! number of soil layers [-] 
  real,    intent(in) :: sldpth(nsoil)        ! thickness of soil layers [-]  
  integer, intent(in) :: nsnow                ! maximum number of snow layers (e.g. 3) [-] 
  real,    intent(in) :: shdfac_monthly(12)   ! monthly values for green vegetation fraction [-]
  integer, intent(in) :: vegetype               ! land cover type index [-] 
  integer, intent(in) :: soiltype             ! soil type index [-] 
  integer, intent(in) :: slopetype            ! slope type for Noah baseflow [-] 
  integer, intent(in) :: urban_vegetype         ! urban land cover type index [-]   
  integer, intent(in) :: ice_flag             ! ice flag: 0 = no ice, 1 = ice  (default: 0) [-]
  integer, intent(in) :: st_flag              ! surface type 1=soil, 2=lake (default: 1) [-] 
  integer, intent(in) :: sc_idx               ! soil color type, an integer index from 1 to 8 (default: 4) [-]  
  integer, intent(in) :: iz0tlnd              ! option of Chen adjustment of Czil (default: 0-not use, other wise 1) [-]   
  real,    intent(in) :: smceq(nsoil)         ! equilibrium soil water content [m^3 m-3]
  real,    intent(in) :: tair                 ! air temperature [K] 
  real,    intent(in) :: psurf                ! air pressure [Pa] 
  real,    intent(in) :: wind_e               ! eastard wind speed [m s-1] 
  real,    intent(in) :: wind_n               ! northward wind speed [m s-1] 
  real,    intent(in) :: qair                 ! near Surface Specific Humidity [kg kg-1] 
  real,    intent(in) :: swdown               ! downward solar radiation [w/m2]
  real,    intent(in) :: lwdown               ! downward longwave radiation [w/m2]
  real,    intent(in) :: prcp                 ! total precip (rainfall+snowfall) Rate [kg m-2 s-1]   
  !real,    intent(in) :: rainf                ! rainfall Rate [kg m-2 s-1]   
  !real,    intent(in) :: snowf                ! snowfall Rate [kg m-2 s-1] 
  real,    intent(in) :: tbot                 ! deep-layer soil temperature [K]
  real,    intent(in) :: pblh                 ! planetary boundary layer height [m] 
  real, intent(inout) :: zlvl                 ! reference height of temperature and humidity [m] 
  ! SY: Begin for enabling OPTUE, added 02/08/18. 
  real,    intent(in) :: p_csoil                ! vol. soil heat capacity [j/m3/K]
  real,    intent(in) :: p_bexp                 ! B parameter
  real,    intent(in) :: p_dksat                ! saturated soil hydraulic conductivity
  real,    intent(in) :: p_dwsat                ! saturated soil hydraulic diffusivity
  real,    intent(in) :: p_psisat               ! saturated soil matric potential
  real,    intent(in) :: p_quartz               ! soil quartz content
  real,    intent(in) :: p_smcmax               ! porosity, saturated value of soil moisture (volumetric)
  real,    intent(in) :: p_smcref               ! reference soil moisture (field capacity)
  real,    intent(in) :: p_smcwlt               ! wilting point soil moisture (volumetric)
  real,    intent(in) :: p_czil                 ! Calculate roughness length of heat
  real,    intent(in) :: p_frzk                 ! frozen ground parameter
  real,    intent(in) :: p_refdk                ! parameters in the surface runoff parameteriz.
  real,    intent(in) :: p_refkdt               ! parameters in the surface runoff parameteriz.
  real,    intent(in) :: p_slope                ! slope index (0 - 1)
  real,    intent(in) :: p_topt                 ! optimum transpiration air temperature
  real,    intent(in) :: p_rgl                  ! parameter used in radiation stress function
  real,    intent(in) :: p_rsmax                ! maximum stomatal resistance
  real,    intent(in) :: p_rsmin                ! minimum Canopy Resistance [s/m]
  real,    intent(in) :: p_hs                   ! parameter used in vapor pressure deficit function
  real,    intent(in) :: p_nroot
  real,    intent(in) :: p_CH2OP                ! maximum intercepted h2o per unit lai+sai [mm]
  real,    intent(in) :: p_DLEAF                ! characteristic leaf dimension [m]
  real,    intent(in) :: p_Z0MVT                ! momentum roughness length [m]
  real,    intent(in) :: p_HVT                  ! top of canopy [m]
  real,    intent(in) :: p_HVB                  ! bottom of canopy [m]
  real,    intent(in) :: p_RC                   ! tree crown radius [m]
  real,    intent(in) :: p_RHOL1                ! leaf reflectance (1=vis)
  real,    intent(in) :: p_RHOL2                ! leaf reflectance (2=nir)
  real,    intent(in) :: p_RHOS1                ! stem reflectance (1=vis)
  real,    intent(in) :: p_RHOS2                ! stem reflectance (2=nir)
  real,    intent(in) :: p_TAUL1                ! leaf transmittance (1=vis)
  real,    intent(in) :: p_TAUL2                ! leaf transmittance (2=nir)
  real,    intent(in) :: p_TAUS1                ! stem transmittance (1=vis)
  real,    intent(in) :: p_TAUS2                ! stem transmittance (2=nir)
  real,    intent(in) :: p_XL                   ! leaf/stem orientation index
  real,    intent(in) :: p_CWPVT                ! empirical canopy wind parameter
  real,    intent(in) :: p_C3PSN                ! photosynthetic pathway (0. = c4, 1. = c3)
  real,    intent(in) :: p_KC25                 ! co2 michaelis-menten constant at 25c [pa]
  real,    intent(in) :: p_AKC                  ! q10 for kc25
  real,    intent(in) :: p_KO25                 ! o2 michaelis-menten constant at 25c [pa]
  real,    intent(in) :: p_AKO                  ! q10 for ko25
  real,    intent(in) :: p_AVCMX                ! q10 for vcmx25
  real,    intent(in) :: p_AQE                  ! q10 for qe25
  real,    intent(in) :: p_LTOVRC               ! leaf turnover [1/s]
  real,    intent(in) :: p_DILEFC               ! coeficient for leaf stress death [1/s]
  real,    intent(in) :: p_DILEFW               ! coeficient for leaf stress death [1/s]
  real,    intent(in) :: p_RMF25                ! leaf maintenance respiration at 25c [umol co2/m**2/s]
  real,    intent(in) :: p_SLA                  ! single-side leaf area per Kg [m2/kg]
  real,    intent(in) :: p_FRAGR                ! fraction of growth respiration. original was 0.3
  real,    intent(in) :: p_TMIN                 ! minimum temperature for photosynthesis [k]
  real,    intent(in) :: p_VCMX25               ! maximum rate of carboxylation at 25c. unit [umol co2/m**2/s]
  real,    intent(in) :: p_TDLEF                ! characteristic T for leaf freezing [K]
  real,    intent(in) :: p_BP                   ! minimum leaf conductance [umol/m**2/s]
  real,    intent(in) :: p_MP                   ! slope of conductance-to-photosynthesis relationship
  real,    intent(in) :: p_QE25                 ! quantum efficiency at 25c [umol co2 / umol photon]
  real,    intent(in) :: p_RMS25                ! stem maintenance respiration at 25c [umol co2/kg bio/s]
  real,    intent(in) :: p_RMR25                ! root maintenance respiration at 25c [umol co2/kg bio/s]
  real,    intent(in) :: p_ARM                  ! q10 for maintenance respiration
  real,    intent(in) :: p_FOLNMX               ! foliage nitrogen concentration when f(n)=1 [%]
  real,    intent(in) :: p_WDPOOL               ! wood pool (switch 1 or 0) depending on woody or not [-]
  real,    intent(in) :: p_WRRAT                ! wood to non-wood ratio
  real,    intent(in) :: p_MRP                  ! microbial respiration parameter [umol co2 /kg c/ s]
  ! SY: End for enabling OPTUE, added 02/08/18
  real, intent(inout) :: albold               ! snow albedo at last time step [-]
  real, intent(inout) :: sneqvo               ! snow mass at the last time step [mm] 
  real, intent(inout) :: sstc(nsnow+nsoil)    ! snow/soil temperature [K]
  real, intent(inout) :: sh2o(nsoil)          ! volumetric liquid soil moisture [m^3 m-3]
  real, intent(inout) :: smc(nsoil)           ! volumetric soil moisture, ice + liquid [m^3 m-3] 
  real, intent(inout) :: tah                  ! canopy air temperature [K]  
  real, intent(inout) :: eah                  ! canopy air vapor pressure [Pa]  
  real, intent(inout) :: fwet                 ! wetted or snowed fraction of canopy [-]  
  real, intent(inout) :: canliq               ! intercpted liquid water [mm]
  real, intent(inout) :: canice               ! intercepted ice mass [mm]
  real, intent(inout) :: tv                   ! vegetation temperature [K] 
  real, intent(inout) :: tg                   ! ground temperature (skin temperature) [K]
  real, intent(inout) :: qsnow                ! snowfall on the ground [mm s-1]
  integer, intent(inout) :: isnow             ! actual number of snow layers [-]
  real, intent(inout) :: zss(nsnow+nsoil)     ! snow/soil layer-bottom depth from snow surface [m]
  real, intent(inout) :: snowh                ! snow height [m] 
  real, intent(inout) :: sneqv                ! snow water equivalent [mm] 
  real, intent(inout) :: snowice(nsnow)       ! snow-layer ice [mm] 
  real, intent(inout) :: snowliq(nsnow)       ! snow-layer liquid water [mm] 
  real, intent(inout) :: zwt                  ! depth to water table [m]
  real, intent(inout) :: wa                   ! water storage in aquifer [mm]
  real, intent(inout) :: wt                   ! water in aquifer and saturated soil [mm]
  real, intent(inout) :: wslake               ! lake water storage (can be negative) [mm]
  real, intent(inout) :: lfmass               ! leaf mass (used only for dveg_opt=2)  [g/m2]
  real, intent(inout) :: rtmass               ! mass of fine roots [g/m2] 
  real, intent(inout) :: stmass               ! stem mass [g/m2] 
  real, intent(inout) :: wood                 ! mass of wood (including woody roots) [g/m2]
  real, intent(inout) :: stblcp               ! stable carbon in deep soil [g/m2]
  real, intent(inout) :: fastcp               ! short-lived carbon in shallow soil [g/m2] 
  real, intent(inout) :: lai                  ! leaf area index [-]
  real, intent(inout) :: sai                  ! stem area index [-] 
  real, intent(inout) :: cm                   ! momentum drag coefficient [s/m] 
  real, intent(inout) :: ch                   ! sensible heat exchange coefficient [s/m] 
  real, intent(inout) :: tauss                ! snow aging term [-] 
  real, intent(inout) :: smcwtd               ! soil water content between bottom of the soil and water table [m^3 m-3]
  real, intent(inout) :: deeprech             ! recharge to or from the water table when deep [m]
  real, intent(inout) :: rech                 ! recharge to or from the water table when shallow [m]
  real,   intent(out) :: fsa                  ! total absorbed solar radiation [W m-2]
  real,   intent(out) :: fsr                  ! total reflected solar radiation [W m-2] 
  real,   intent(out) :: fira                 ! total net longwave radiation to atmosphere [W m-2] 
  real,   intent(out) :: fsh                  ! total sensible heat to atmosphere [W m-2]
  real,   intent(out) :: ssoil                ! ground heat flux to soil [W m-2]
  real,   intent(out) :: fcev                 ! canopy evaporative heat to atmosphere [W m-2]
  real,   intent(out) :: fgev                 ! ground evaporative heat to atmosphere [W m-2]  
  real,   intent(out) :: fctr                 ! transpiration heat to atmosphere [W m-2]
  real,   intent(out) :: ecan                 ! evaporation rate of canopy water [kg m-2 s-1] 
  real,   intent(out) :: etran                ! transpiration rate [kg m-2 s-1] 
  real,   intent(out) :: edir                 ! direct evaporation rate from surface [kg m-2 s-1] 
  real,   intent(out) :: trad                 ! surface radiative temperature [K]  
  real,   intent(out) :: subsnow              ! snow sublimation rate [kg m-2 s-1]
  real,   intent(out) :: tgb                  ! ground temperature (K)
  real,   intent(out) :: tgv                  ! ground surface temperature [K]
  real,   intent(out) :: t2mv                 ! 2-m air temperature over vegetated part [K]
  real,   intent(out) :: t2mb                 ! 2-m height air temperature [K] 
  real,   intent(out) :: q2v                  ! 2-m specific humidity over vegetation [kg kg-1]
  real,   intent(out) :: q2b                  ! 2-m air specific humidity [kg kg-1]
  real,   intent(out) :: runsrf               ! surface runoff [kg m-2 s-1]
  real,   intent(out) :: runsub               ! baseflow (saturation excess) [kg m-2 s-1] 
  real,   intent(out) :: apar                 ! photosynthesis active energy by canopy [W m-2]
  real,   intent(out) :: psn                  ! total photosynthesis of CO2 [umol m-2 s-1] 
  real,   intent(out) :: sav                  ! solar radiation absorbed by vegetation [W m-2] 
  real,   intent(out) :: sag                  ! solar radiation absorbed by ground [W m-2] 
  real,   intent(out) :: fsno                 ! snow-cover fraction on the ground [-]
  real,   intent(out) :: nee                  ! net ecosystem exchange of CO2 [g/m2s] 
  real,   intent(out) :: gpp                  ! net instantaneous assimilation of carbon [g/m2s] 
  real,   intent(out) :: npp                  ! net primary productivity of carbon [g/m2s]
  real,   intent(out) :: fveg                 ! green vegetation fraction [-]
  real,   intent(out) :: albedo               ! surface albedo [-] 
  real,   intent(out) :: qsnbot               ! melting water out of snow bottom [kg m-2 s-1]
  real,   intent(out) :: ponding              ! surface ponding [mm]
  real,   intent(out) :: ponding1             ! surface ponding1 [mm]
  real,   intent(out) :: ponding2             ! surface ponding2 [mm] 
  real,   intent(out) :: rssun                ! sunlit stomatal resistance (s/m)
  real,   intent(out) :: rssha                ! shaded stomatal resistance (s/m)
  real,   intent(out) :: bgap                 ! between canopy gap fraction for beam [-]
  real,   intent(out) :: wgap                 ! within canopy gap fraction for beam [-] 
  real,   intent(out) :: chv                  ! sensible heat exchange coeffieient over vegetated fraction [s/m]
  real,   intent(out) :: chb                  ! sensible heat exchange coefficient over bare-ground fraction [s/m]
  real,   intent(out) :: emissi               ! surface emissivity [-] 
  real,   intent(out) :: shg                  ! ground sensible heat (+ to atm) [W m-2]     
  real,   intent(out) :: shc                  ! canopy sensible heat (+ to atm) [W m-2]   
  real,   intent(out) :: shb                  ! bare ground sensible heat (+ to atm) [W m-2]     
  real,   intent(out) :: evg                  ! ground evaporation heat (+ to atm) [W m-2]  
  real,   intent(out) :: evb                  ! bare ground evaporation heat (+ to atm) [W m-2]  
  real,   intent(out) :: ghv                  ! ground heat flux (+ to soil]) [W m-2] 
  real,   intent(out) :: ghb                  ! bare ground heat flux (+ to soil) [W m-2] 
  real,   intent(out) :: irg                  ! ground net long wave radiation (+ to atm) [W m-2] 
  real,   intent(out) :: irc                  ! canopy net long wave radiation (+ to atm) [W m-2] 
  real,   intent(out) :: irb                  ! bare ground net long wave radiation (+ to atm) [W m-2] 
  real,   intent(out) :: tr                   ! transpiration heat (+ to atm) [W m-2] 
  real,   intent(out) :: evc                  ! canopy evaporation heat (+ to atm) [W m-2] 
  real,   intent(out) :: chleaf               ! leaf exchange coefficient [-]  
  real,   intent(out) :: chuc                 ! under canopy exchange coefficient  [-]
  real,   intent(out) :: chv2                 ! sensible heat exchange coefficient over vegetated fraction [-] 
  real,   intent(out) :: chb2                 ! sensible heat exchange coefficient over bare-ground [-] 
  real,   intent(out) :: fpice                ! snow fraction in precipitation [-] 
  !ag (12Sep2019)
  real, intent(inout) :: rivsto               ! river storage
  real, intent(inout) :: fldsto               ! flood storage
  real, intent(inout) :: fldfrc               ! flood storage
  real, intent(inout) :: sfcheadrt            ! extra output for WRF-HYDRO [m] 

  ! external function
  real, external      :: month_d_36

  ! local variables 
  integer             :: iloc, jloc 
  real                :: dx
  integer             :: isurban
  character(len=12)   :: nowdate
  integer :: n 
  integer :: ice
  integer :: ist
  integer :: isc
  real    :: q2
  real    :: sfctmp
  real    :: uu
  real    :: vv
  real    :: soldn
  real    :: lwdn
  !real    :: prcp
  real    :: co2air
  real    :: o2air
  real    :: cosz
  real    :: foln
  real    :: sfcprs
  integer :: yearlen
  real    :: julian
  real    :: shdfac
  real    :: shdmax
  real    :: lat
  real    :: z0
  real, allocatable, dimension(:) :: zsoil
  real, allocatable, dimension(:) :: ficeold
  real, allocatable, dimension(:) :: zsnso
  real, allocatable, dimension(:) :: snice
  real, allocatable, dimension(:) :: snliq
  real, allocatable, dimension(:) :: stc

  real    :: qfx
  real    :: flxsum
  real    :: ir_sh_ev_gh
  real    :: fsa_fsr
  real    :: dz8w
  real    :: qsfc
  real    :: qc             ! cloud water mixing ratio, unit [-], not actually used? 
  real    :: psfc
  integer :: snl_idx
#if(LIS_NoahMP_TEST)
  write(*,*) " Noah model landuse parameter table: ",  landuse_tbl_name    
  write(*,*) " Noah model soil parameter table: ",     soil_tbl_name       
  write(*,*) " Noah model general parameter table: ",  gen_tbl_name        
  write(*,*) " NoahMP parameter table: ",              noahmp_tbl_name     
  write(*,*) " Landuse classficiation scheme: ",       landuse_scheme_name 
  write(*,*) " Soil classification scheme: ",          soil_scheme_name    
#endif

  !!!! kludge fixes by Shugong Wang to
  !!!!   prevent crashes (next three lines)
  if (wood .gt. 2500.0)  wood=2000.0
  if (wood .lt. 0.0) wood = 0.0 
  if (tv .gt. tg+20.0) tv = tg+20.0  
  
  slcats   = nslcats     
  lucats   = nlucats     
  slpcats  = nslpcats 
  
  ! 12 MP options   
  dveg     = dveg_opt 
  opt_crs  = crs_opt  
  opt_btr  = btr_opt  
  opt_run  = run_opt  
  opt_sfc  = sfc_opt  
  opt_frz  = frz_opt  
  opt_inf  = inf_opt  
  opt_rad  = rad_opt  
  opt_alb  = alb_opt  
  opt_snf  = snf_opt  
  opt_tbot = tbot_opt 
  opt_stc  = stc_opt  


  ! dummy arguments for Noah-MP 
!  iloc = 1
!  jloc = 1 
#ifndef WRF_HYDRO
  sfcheadrt = 0.0 
#endif

  if (opt_sfc == 3) then
     stop "(opt_sfc == 3) and (opt_sfc == 4) are not for offline use"
  endif

  ist     = st_flag   ! surface type:  ist=1 => soil;  ist=2 => lake
  isc     = sc_idx    ! soil color type
  ice     = ice_flag  ! surface type:  ice=0 => soil;  ice=1 => sea-ice

  ! set ZSOIL 
  allocate(zsoil(nsoil))
  ! zsoil is negative.
  zsoil(1) = -sldpth(1)
  do n = 2, nsoil
     zsoil(n) = zsoil(n-1) - sldpth(n)
  enddo

  !!!!!! convert latitude from decimal degree into radius  
  lat     = latitude * (3.1415926535/180.0)
  ! dx is horizontal grid spacing. dx should not be used. set it to a dummy value 1.0 first. 
  dx      = 1.0  
  
  isurban = urban_vegetype  
  
  allocate(ficeold(-nsnow+1:0))
  allocate(zsnso(-nsnow+1:nsoil))
  allocate(snice(-nsnow+1:0))
  allocate(snliq(-nsnow+1:0))
  allocate(stc(-nsnow+1:nsoil))


  zsnso = 0.0
  snice = 0.0
  snliq = 0.0
  stc   = 0.0
  ficeold = 0.0 

  ! state variables 
  zsnso = 0.0
  snice = 0.0
  snliq = 0.0
  stc   = 0.0
  ficeold = 0.0

  zsnso(-nsnow+1:nsoil) = zss(1:nsnow+nsoil) 
  snice(-nsnow+1:0)     = snowice(1:nsnow)
  snliq(-nsnow+1:0)     = snowliq(1:nsnow) 
  stc(-nsnow+1:nsoil)   = sstc(1:nsnow+nsoil) 
  ! ice fraction at the last timestep, add check for both snice and snliq are 0.0
  do snl_idx=isnow+1,0
    if(snice(snl_idx)+snliq(snl_idx)>0.0) then
      ficeold(snl_idx)  = snice(snl_idx) / (snice(snl_idx)+snliq(snl_idx))
    else 
      ficeold(snl_idx)  = 0.0
    endif
  enddo

  ! ficeold(isnow+1:0)    = snice(isnow+1:0) / ( snice(isnow+1:0)+snliq(isnow+1:0) ) ! ice fraction at the last timestep

  ! cosz, yearlen, and julian are calculated in subroutine calc_declin_36 
  ! be careful here!!!, LIS uses GMT. the date for calc_declin_36 should be local time. Longitude is need to cnvert GMT 
  ! into local time!!!. To be implemented after bondville test 
  write(nowdate,'(I4.4,4I2.2)') year, month, day, hour, minute
  call calc_declin_36(nowdate(1:4)//"-"//nowdate(5:6)//"-"//nowdate(7:8)//"_"//nowdate(9:10)//":"//nowdate(11:12)//":00", &
      latitude, longitude, cosz, yearlen, julian)

  ! SY: Begin for enabling OPTUE
  ! SY: Begin lines following those in REDPRM and NoahMP36_setup
  ! SY: Begin SOIL PARAMETERS
  CSOIL_DATA = p_csoil
  BB(soiltype) = p_bexp
  SATDK(soiltype) = p_dksat
  SATDW(soiltype) = p_dwsat
  SATPSI(soiltype) = p_psisat
  QTZ(soiltype) = p_quartz
  MAXSMC(soiltype) = p_smcmax
  REFSMC(soiltype) = p_smcref
  WLTSMC(soiltype) = p_smcwlt
  ! SY: End SOIL PARAMETERS
  ! SY: Begin UNIVERSAL PARAMETERS
  CZIL_DATA = p_czil
  FRZK_DATA = p_frzk
  REFDK_DATA = p_refdk
  REFKDT_DATA = p_refkdt
  SLOPE_DATA(slopetype) = p_slope
  ! SY: End UNIVERSAL PARAMETERS
  ! SY: Begin VEGETATION PARAMETERS
  TOPT_DATA = p_topt
  RGLTBL(vegetype) = p_rgl   
  RSMAX_DATA = p_rsmax
  RSTBL(vegetype) = p_rsmin   
  HSTBL(vegetype) = p_hs  
  NROTBL(vegetype) = p_nroot  
  ! SY: End VEGETATION PARAMETERS  
  ! SY: End lines following those in REDPRM and NoahMP36_setup
  ! SY: Begin lines following those in read_mp_veg_parameters and NoahMP36_setup
  CH2OP(vegetype) = p_CH2OP
  DLEAF(vegetype) = p_DLEAF
  Z0MVT(vegetype) = p_Z0MVT
  HVT(vegetype) = p_HVT
  HVB(vegetype) = p_HVB
  RC(vegetype) = p_RC
  RHOL(vegetype,1) = p_RHOL1
  RHOL(vegetype,2) = p_RHOL2
  RHOS(vegetype,1) = p_RHOS1
  RHOS(vegetype,2) = p_RHOS2
  TAUL(vegetype,1) = p_TAUL1
  TAUL(vegetype,2) = p_TAUL2
  TAUS(vegetype,1) = p_TAUS1
  TAUS(vegetype,2) = p_TAUS2
  XL(vegetype) = p_XL
  CWPVT(vegetype) = p_CWPVT
  C3PSN(vegetype) = p_C3PSN
  KC25(vegetype) = p_KC25
  AKC(vegetype) = p_AKC
  KO25(vegetype) = p_KO25
  AKO(vegetype) = p_AKO
  AVCMX(vegetype) = p_AVCMX
  AQE(vegetype) = p_AQE
  LTOVRC(vegetype) = p_LTOVRC
  DILEFC(vegetype) = p_DILEFC
  DILEFW(vegetype) = p_DILEFW
  RMF25(vegetype) = p_RMF25
  SLA(vegetype) = p_SLA
  FRAGR(vegetype) = p_FRAGR
  TMIN(vegetype) = p_TMIN
  VCMX25(vegetype) = p_VCMX25
  TDLEF(vegetype) = p_TDLEF
  BP(vegetype) = p_BP
  MP(vegetype) = p_MP
  QE25(vegetype) = p_QE25
  RMS25(vegetype) = p_RMS25
  RMR25(vegetype) = p_RMR25
  ARM(vegetype) = p_ARM
  FOLNMX(vegetype) = p_FOLNMX
  WDPOOL(vegetype) = p_WDPOOL
  WRRAT(vegetype) = p_WRRAT
  MRP(vegetype) = p_MRP
  ! SY: End lines following those in read_mp_veg_parameters and NoahMP36_setup
  ! SY: End for enabling OPTUE

  call redprm (vegetype, soiltype, slopetype, zsoil, nsoil, isurban)


  if ( dveg == 1 ) then
    ! with dveg==1, shdfac is fed directly to fveg
    shdfac = month_d_36(shdfac_monthly, nowdate)
  else
    ! with dveg==2, fveg is computed from lai and sai, and shdfac is unused
    shdfac = -1.e38
  endif
  shdmax = maxval(shdfac_monthly)

  ! assign forcing variables 
  sfctmp = tair 
  sfcprs = psurf 
  q2     = qair 
  soldn  = swdown
  lwdn   = lwdown
  dz8w = -1.e36 ! not used
  qsfc = q2 !?
  psfc = sfcprs !?
  uu = wind_e
  vv = wind_n 
  co2air  = 355.e-6 * sfcprs ! partial pressure of co2 (pa) ! from noah-mp-wrf
  o2air   = 0.209   * sfcprs ! partial pressure of o2 (pa)  ! from noah-mp-wrf
  foln    = 1.0              ! foliage nitrogen, (fraction) a kind of forcing from WRF? 
                             ! in simple driver, it is set to 1.0. In NoahMP code, it 
                             ! is just input (check intent) 
  call  noahmp_sflx_36 (                                                   &
               iloc    , jloc    , lat     , yearlen , julian  , cosz    , & ! in : time/space-related
               dt      , dx      , upd_alb_flag,albd_in, albi_in, &
               dz8w    , nsoil   , zsoil   , nsnow   , & ! in : model configuration 
               shdfac  , shdmax  , vegetype  , isurban , ice   , ist     , & ! in : vegetation/soil characteristics
               isc     , smceq   ,                                         & ! in : vegetation/soil characteristics
               iz0tlnd ,                                                   & ! in : user options
               sfctmp  , sfcprs  , psfc    , uu      , vv      , q2      , & ! in : forcing
               qc      , soldn   , lwdn    , prcp    , tbot    , co2air  , & ! in : forcing
               o2air   , foln    , ficeold , pblh    , zlvl,               & ! in : forcing
               albold  , sneqvo  ,                                         & ! in/out : 
               stc     , sh2o    , smc     , tah     , eah     , fwet    , & ! in/out : 
               canliq  , canice  , tv      , tg      , qsfc    , qsnow   , & ! in/out : 
               isnow   , zsnso   , snowh   , sneqv   , snice   , snliq   , & ! in/out : 
               zwt     , wa      , wt      , wslake  , lfmass  , rtmass  , & ! in/out : 
               stmass  , wood    , stblcp  , fastcp  , lai     , sai     , & ! in/out : 
               cm      , ch      , tauss   ,                               & ! in/out : 
               smcwtd  , deeprech, rech    ,                               & ! in/out :
               fsa     , fsr     , fira    , fsh     , ssoil   , fcev    , & ! out : 
               fgev    , fctr    , ecan    , etran   , edir    , trad    , & ! out :
               subsnow ,                                                   & ! out :
               tgb     , tgv     , t2mv    , t2mb    , q2v     , q2b     , & ! out :
               runsrf  , runsub  , apar    , psn     , sav     , sag     , & ! out :
               fsno    , nee     , gpp     , npp     , fveg    , albedo  , & ! out :
               qsnbot  , ponding , ponding1, ponding2, rssun   , rssha   , & ! out :
               bgap    , wgap    , chv     , chb     , emissi  ,           & ! out :
               shg     , shc     , shb     , evg     , evb     , ghv     , & ! out :
               ghb     , irg     , irc     , irb     , tr      , evc     , & ! out :
               chleaf  , chuc    , chv2    , chb2    , fpice   , &
               !ag (12Sep2019)
               rivsto  , fldsto, fldfrc            &
#ifdef WRF_HYDRO
               ,sfcheadrt                                                  & ! in/out :
#endif
               )


  qfx = fgev + fcev + fctr
  ir_sh_ev_gh = fira + fsh + qfx + ssoil
  flxsum = (-fsa) + fira + fsh + qfx + ssoil
  fsa_fsr = fsa + fsr
  
  ! state variables 
  sstc(1:nsnow+nsoil) = stc(-nsnow+1:nsoil) 
  zss(1:nsnow+nsoil)  = zsnso(-nsnow+1:nsoil) 
  snowice(1:nsnow)    = snice(-nsnow+1:0)     
  snowliq(1:nsnow)    = snliq(-nsnow+1:0)     
  
  deallocate(zsoil)
  deallocate(ficeold)
  deallocate(zsnso)
  deallocate(snice)
  deallocate(snliq)
  deallocate(stc)
end subroutine noahmp_driver_36

real function month_d_36(a12, nowdate) result (nowval)
  !
  ! Given a set of 12 values, taken to be valid on the fifteenth of each month (Jan through Dec)
  ! and a date in the form <YYYYMMDD[HHmmss]> ....
  ! 
  ! Return a value valid for the day given in <nowdate>, as an interpolation from the 12
  ! monthly values.
  !
  use kwm_date_utilities_36
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

end function month_d_36 

SUBROUTINE calc_declin_36 ( nowdate, latitude, longitude, cosz, yearlen, julian)
  use kwm_date_utilities_36 
!---------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------------------

  REAL, PARAMETER :: DEGRAD = 3.14159265/180.
  REAL, PARAMETER :: DPD    = 360./365.
! !ARGUMENTS:
  character(len=19), intent(in)  :: nowdate    ! YYYY-MM-DD_HH:mm:ss
  real,              intent(in)  :: latitude
  real,              intent(in)  :: longitude
  real,              intent(out) :: cosz
  integer,           intent(out) :: yearlen
  real,              intent(out) :: JULIAN

  REAL                           :: hrang
  real                           :: DECLIN
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
  julian = real(iday) + real(ihour)/24.

!
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

END SUBROUTINE calc_declin_36 

! Subroutine SNOW_INIT grabbed from NOAH-MP-WRF
SUBROUTINE SNOW_INIT_36 ( jts, jtf, its, itf, ims, ime, jms, jme, NSNOW, NSOIL, ZSOIL,  &
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

END SUBROUTINE SNOW_INIT_36
