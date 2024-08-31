!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: crocus_driver
! \label{crocus_driver}
!
! !REVISION HISTORY:
! July 17 2019 Mahdi Navari started for implementing Crocus
! !INTERFACE:!
!
SUBROUTINE crocus_driver(n, &
                         nsnow, &
                         nimpur, &
                         !ttile,                  &
                         !latitude,longitude,     &
                         year, &
                         month, &
                         day, &
                         hour, &
                         minute, &
                         SNOWRES_opt, &
                         !TPTIME,                 &
                         OMEB_BOOL, &
                         GLACIER_BOOL, &
                         HIMPLICIT_WIND_opt, &
                         SNOWSWE, &
                         SNOWRHO, &
                         SNOWHEAT, &
                         SNOWALB, &
                         SNOWGRAN1, &
                         SNOWGRAN2, &
                         SNOWHIST, &
                         SNOWAGE, &
                         PTSTEP, &
                         PPS, &
                         SRSNOW, &
                         RRSNOW, &
                         TA, &
                         TG, & 
                         SW_RAD, &
                         QA, &
                         Wind_E, & ! Eastward Wind
                         Wind_N, & ! Northward Wind
                         LW_RAD, &
                         UREF, &
                         SLOPE, & ! replace ZP_DIRCOSZW with SLOPE and compute the cosine in the driver
                         ZREF, &
                         Z0NAT, &
                         Z0EFF, &
                         Z0HNAT, &
                         ALB, &
                         D_G, &
                         SNOWLIQ, &
                         SNOWTEMP, &
                         SNOWDZ, &
                         THRUFAL, &
                         GRNDFLUX, &
                         SNDRIFT, &! snowdrift only modifies the  surface snow properties
                         !  (density and microstructure) but without any mass change
                         !  unless you use the  SNOWDRIFT_SUBLIM option
                         !  which can activate sublimation in case of snow drift
                         !  but without any erosion or accumulation fluxes between
                         !  grid points (METEO FRANCE ticket # 1592).
                         RI_n, &
                         EMISNOW, &
                         CDSNOW, &
                         USTARSNOW, &
                         CHSNOW, &
                         SNOWHMASS, &
                         QS, &
                         PERMSNOWFRAC, & ! ZP_VEGTYPE
                         LAT, &
                         LON, &
                         ! ZP_BLOWSNW,          & ! We can not use snow_sytron it is designed for the
                         !  METEO FRANCEinternal use .(METEO FRANCE ticket # 1592).
                         !  Therefore it will be defined as a local variable.
                         SNOWDRIFT_opt, & !IO%CSNOWDRIFT
                         SNOWDRIFT_SUBLIM_BOOL, & ! IO%LSNOWDRIFT_SUBLIM
                         SNOW_ABS_ZENITH_BOOL, &! IO%LSNOW_ABS_ZENITH
                         SNOWMETAMO_opt, & !IO%CSNOWMETAMO
                         SNOWRAD_opt, &! IO%CSNOWRAD
                         ATMORAD_BOOL, &!IO%LATMORAD
                         IMPWET, &
                         IMPDRY, &
                         SNOWFALL_opt, &!IO%CSNOWFALL
                         SNOWCOND_opt, &! IO%CSNOWCOND
                         SNOWHOLD_opt, &! IO%CSNOWHOLD
                         SNOWCOMP_opt, &! IO%CSNOWCOMP
                         SNOWZREF_opt, &! IO%CSNOWZREF
                         SNOWMAK_dz, &!
                         SNOWCOMPACT_BOOL, &! IO%LSNOWCOMPACT_BOOL
                         SNOWMAK_BOOL, &! IO%LSNOWMAK_BOOL
                         SNOWTILLER_BOOL, &! IO%LSNOWTILLER
                         SELF_PROD_BOOL, &! IO%LSELF_PROD
                         SNOWMAK_PROP_BOOL, &! IO%LSNOWMAK_PROP
                         PRODSNOWMAK_BOOL, &! IO%LPRODSNOWMAK)
                         SLOPE_DIR, &  ! IN    - !Typical slope aspect in the grid  (from N clockwise) [Radians]
                         SAND              , & ! IN    - Soil sand fraction (-) [-]
                         SILT              , & ! IN    - Soil silt fraction (-) [-]
                         CLAY              , & ! IN    - Soil clay fraction (-) [-]
                         POROSITY          , & ! IN    - Soil porosity (m3 m-3) [m3/m3]
                         XWGI              ,&  ! IN    - soil volumetric frozen water content
                         XWG               ,&  ! IN    - soil volumetric liquid water content
                         usemonalb )

   USE LIS_coreMod, only: LIS_rc
   !USE LIS_logMod,  only  : LIS_logunit, LIS_endrun
   !USEe LIS_timeMgrMod, only : LIS_date2time, LIS_tick
   USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME! MN: declaration of temporal types.  DATE_TIME%.... (YEAR, MONTH, DAY)
   USE MODD_CSTS, ONLY: XTT, XPI, XDAY, XLMTT, XLSTT, &
                        XG, XBOLTZ, XAVOGADRO, XMD, &
                        XMV, XRD, XRV, XCPD, XP00 ! need these to compute air density
   USE modi_surface_cd ! for drag coefficient for momentum
   USE modi_surface_aero_cond  ! for drag coefficient for heat
   USE modd_snow_metamo ! use XVVISC3= 0.023
   USE MODD_SNOW_PAR, ONLY: XRHOSMAX_ES, XSNOWDMIN, XRHOSMIN_ES, XEMISSN, &
      XRHO_SNOWMAK, XPSR_SNOWMAK, XPTA_SEUIL, &
      XPROD_SCHEME, XPROD_COUNT, XTIMESNOWMAK
   USE MODD_SURF_ATM, ONLY: XRIMAX
   USE MODI_INI_SURF_CSTS
   USE MODD_SURF_PAR, ONLY: XUNDEF
   USE MODD_CONST_ATM, ONLY: JPNBANDS_ATM
   USE MODE_SNOW3L
   USE MODE_TARTES, ONLY: SNOWCRO_TARTES
   USE MODE_THERMOS
   USE MODE_CRODEBUG
   USE MODD_REPROD_OPER
   USE MODD_ISBA_PAR,   ONLY : XDRYWGHT, XSPHSOIL, XCONDQRTZ, XCONDOTH1, XCONDOTH2 ! for soil thermal conductivity  

   implicit none

   TYPE(DATE_TIME)     :: TPTIME           ! current date and time
!*  0.1  declarations of arguments
   integer, intent(in) :: n                ! nest id
   integer, intent(in) :: nsnow            ! number of snow layers [-]
   integer, intent(in) :: nimpur           ! number of impurities [-]
   integer, intent(in) :: year             ! year of the current time step [-]
   integer, intent(in) :: month            ! month of the current time step [-]
   integer, intent(in) :: day              ! day of the current time step [-]
   integer, intent(in) :: hour             ! hour of the current time step [-]
   integer, intent(in) :: minute           ! minute of the current time step [-]

   CHARACTER(LEN=3), INTENT(IN)   :: SNOWRES_opt
!                         SNOWRES_opt  = ISBA-SNOW3L turbulant exchange option
!                         'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!                         'RIL' = Limit Richarson number under very stable
!                                  conditions (currently testing)
!                         'M98'  = Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus
   LOGICAL, INTENT(IN)            :: OMEB_BOOL
!                         True = coupled to MEB. This means surface fluxes ae IMPOSED
!                         as an upper boundary condition to the explicit snow schemes.
!                         If = False, then energy
!                         budget and fluxes are computed herein.
   LOGICAL, INTENT(IN)            :: GLACIER_BOOL
!                         True = Over permanent snow and ice,
!                         initialise WGI=WSAT,
!                         Hsnow>=10m and allow 0.8<SNOALB<0.85
!                         False = No specific treatment
   CHARACTER(LEN=3), INTENT(IN)  :: HIMPLICIT_WIND_opt
!                         wind implicitation option
!                         'OLD' = direct
!                         'NEW' = Taylor serie, order 1
   REAL*8, INTENT(INOUT)         :: SNOWSWE(nsnow)
!                         Z_PSNOWSWE  = Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
   REAL*8, INTENT(INOUT)         :: SNOWRHO(nsnow)
!                         SNOWRHO  = Snow layer(s) averaged density (kg/m3)
   REAL*8, INTENT(INOUT)         :: SNOWHEAT(nsnow)
!                         SNOWHEAT = Snow layer(s) heat content (J/m2)
   REAL*8, INTENT(INOUT)         :: SNOWALB
!                         SNOWALB = Prognostic surface snow albedo
!                               (does not include anything but
!                               the actual snow cover)
   REAL*8, INTENT(INOUT)         :: SNOWGRAN1(nsnow) !SNOWGRAN1 = Snow layers grain feature 1
   REAL*8, INTENT(INOUT)         :: SNOWGRAN2(nsnow) !SNOWGRAN2 = Snow layer grain feature 2
   REAL*8, INTENT(INOUT)         :: SNOWHIST(nsnow)  ! SNOWHIST  = Snow layer grain historical
!                                 parameter (only for non dendritic snow)
   REAL*8, INTENT(INOUT)         :: SNOWAGE(nsnow)  ! Snow grain age
! REAL,  INTENT(INOUT)         :: ZP_SNOWIMPUR
!                               Snow impurity content (g) (LOCATION,LAYER,NIMPUR)) Impur type :1/BC 2/Dust
!                               CSNOWRAD='B92'  In that case, the direct-diffuse partition is not used, a 3 band
!                               spectral fixed repartition is used from the global radiation and the impact of the
!                               deposition of light absorbing impurities is parameterized from the age of snow
!                               as detailed in Vionnet et al 2012.
   REAL*8, INTENT(IN)  :: PTSTEP ! PTSTEP  = time step of the integration
   REAL*8, INTENT(IN)  :: PPS ! ZP_PS = surface pressure
   REAL*8, INTENT(IN) :: RRSNOW  !  RRSNOW  = rain rate [kg/(m2 s)]
   REAL*8, INTENT(IN) :: SRSNOW  !  SRSNOW = snow rate (SWE) [kg/(m2 s)]
   REAL*8, INTENT(IN) :: TA      ! TA = atmospheric temperature at level za (K)
   REAL*8, INTENT(IN) :: TG      !  TG = Surface soil temperature (effective
!                               temperature the of layer lying below snow)
   REAL*8, INTENT(IN) :: SW_RAD!   SW_RAD = incoming solar radiation (W/m2)
   REAL*8, INTENT(IN) :: QA ! QA = atmospheric specific humidity at level za
   REAL*8, INTENT(IN) :: Wind_E   !  Eastward Wind
   REAL*8, INTENT(IN) :: Wind_N   !  Northward Wind
   REAL*8, INTENT(IN) :: LW_RAD   !  LW_RAD = atmospheric infrared radiation (W/m2)!
   REAL, INTENT(IN)   :: UREF     !  UREF  = reference height of the wind
   REAL, INTENT(IN)   :: SLOPE    !  SLOPE = Slope
   REAL, INTENT(IN)   :: ZREF     !  ZREF  = reference height of the first
!                                  atmospheric level
   REAL, INTENT(IN)    :: Z0NAT    !  Z0NAT (PZ0) = grid box average roughness length
   REAL, INTENT(IN)   :: Z0EFF    !  Z0EFF = roughness length for momentum
   REAL, INTENT(IN)   :: Z0HNAT   !  Z0HNAT (PZOH)  = grid box average roughness length for heat
   REAL*8, DIMENSION (12) :: ALB      !  ALB = soil/vegetation albedo
   REAL, INTENT(IN)   :: D_G      !  D_G  = Assumed first soil layer thickness (m)
!                                  Used to calculate ground/snow heat flux
   REAL*8 , INTENT(IN) :: SAND ! Soil SAND fraction [-]
   REAL*8 , INTENT(IN) :: SILT ! Soil SILT fraction [-]
   REAL*8 , INTENT(IN) :: CLAY ! Soil CLAY fraction [-]
   REAL*8 , INTENT(IN) :: POROSITY ! Soil porosity (m3 m-3) [m3/m3]
   REAL*8 , INTENT(IN) :: XWGI ! soil volumetric frozen water content (m3/m3)
   REAL*8 , INTENT(IN) :: XWG  ! soil volumetric liquid water content (m3/m3)
   REAL*8, INTENT(INOUT)  :: SNOWLIQ(nsnow) !  SNOWLIQ  = Snow layer(s) liquid water content (m)
   REAL*8, INTENT(INOUT)  :: SNOWTEMP(nsnow)!  SNOWTEMP = Snow layer(s) temperature (m)
   REAL*8, INTENT(INOUT)  :: SNOWDZ(nsnow)  !  SNOWDZ = Snow layer(s) thickness (m)
   REAL*8, INTENT(OUT)    :: THRUFAL
!                          THRUFAL  = rate that liquid water leaves snow pack:
!                                paritioned into soil infiltration/runoff by ISBA [kg/(m2 s)]
   REAL*8, INTENT(INOUT)  :: GRNDFLUX
!                         GRNDFLUX = soil/snow interface heat flux (W/m2)
   REAL, INTENT(OUT)      :: SNDRIFT    !   SNDRIFT  = blowing snow sublimation (kg/m2/s)
   REAL, INTENT(OUT)      :: RI_n
!                          RI_n = Ridcharson number (If not OMED initalized to undefined in the snow3L_isba.F90)
   REAL*8, INTENT(OUT)    :: EMISNOW    !   EMISNOW  = snow surface emissivity
   REAL, INTENT(OUT)      :: CDSNOW     !   CDSNOW = drag coefficient for momentum over snow
   REAL, INTENT(OUT)      :: USTARSNOW  !   USTARSNOW  = friction velocity over snow (m/s)
   REAL, INTENT(OUT)      :: CHSNOW     !   CHSNOW = drag coefficient for heat over snow
   REAL*8, INTENT(OUT)    :: SNOWHMASS
!                         SNOWHMASS  = heat content change due to mass
!                              changes in snowpack (J/m2): for budget calculations only.
   REAL*8, INTENT(OUT)     :: QS
!                          QS = surface humidity
   REAL, INTENT(IN)        :: PERMSNOWFRAC
!                          PPERMSNOWFRAC  = fraction of permanet snow/ice
   REAL, INTENT(IN)         :: LAT
   REAL, INTENT(IN)         :: LON
   CHARACTER(4), INTENT(IN):: SNOWDRIFT_opt  ! Snowdrift scheme :
!                          Mechanical transformation of snow grain and compaction + effect of wind
!                          on falling snow properties
!                              'NONE': No snowdrift scheme
!                               'DFLT': falling snow falls as purely dendritic
!                              'GA01': Gallee et al 2001
!                              'VI13': Vionnet et al 2013
   LOGICAL, INTENT(IN)     :: SNOWDRIFT_SUBLIM_BOOL  !   activate sublimation during drift
   LOGICAL, INTENT(IN)     :: SNOW_ABS_ZENITH_BOOL   !   activate parametrization of solar absorption for polar regions
   CHARACTER(LEN=3), INTENT(IN)  :: SNOWMETAMO_opt
!                                B92 (historical version, Brun et al 92), C13, T07, F06 (see Carmagnola et al 2014)
   CHARACTER(LEN=3), INTENT(IN)  :: SNOWRAD_opt
!                               Radiative transfer scheme. HSNOWRAD=B92 Brun et al 1992.
!                               HSNOWRAD=T17 (Tuzet et al. 2017) (Libois et al. 2013)
!                               TARTES with impurities content scheme
   LOGICAL, INTENT(IN)            :: ATMORAD_BOOL !   activate atmotartes scheme
   REAL, INTENT(IN)               :: IMPWET(nimpur)
   REAL, INTENT(IN)               :: IMPDRY(nimpur)
!                                 Dry and wet deposit coefficient from Forcing File(g/m²/s)
   CHARACTER(len=3), INTENT(IN)  :: SNOWFALL_opt
!                                 New options for multiphysics version (Cluzet et al 2016)
!                                 Falling snow scheme
!                                 SNOWFALL_opt=V12 Vionnet et al. 2012 from Brun et al. 1989
!                                 SNOWFALL_opt=A76 Anderson et al. 1976
!                                 SNOWFALL_opt=S02 Lehning el al. 2002
!                                 SNOWFALL_opt=P75 Pahaut 1975
!                                 SNOWFALL_opt=NZE Constant density 200 kg/m3 (who knows ?)
   CHARACTER(len=3), INTENT(IN)  :: SNOWCOND_opt
!                                 Thermal conductivity scheme
!                                 SNOWCOND_opt=Y81 default Crocus from Yen et al. 1981
!                                 SNOWCOND_opt=I02 ISBA_ES snow conductivity parametrization (Boone et al. 2002)
!                                 SNOWCOND_opt=C11 Calonne et al. 2011 snow conductivity parametrization
   CHARACTER(len=3), INTENT(IN)  :: SNOWHOLD_opt
!                                 liquid water content scheme
!                                 SNOWHOLD_opt=B92 default Crocus from Brun et al. 1992 or Vionnet et al. 2012
!                                 SNOWHOLD_opt=B02 ISBA_ES  parametrization (Boone et al. 2002)
!                                 SNOWHOLD_opt=O04 CLM parametrization (Oleson et al 2004)
!                                 SNOWHOLD_opt=S02 SNOWPACK aprametrization (Lehning et al 2002)
   CHARACTER(len=3), INTENT(IN)  :: SNOWCOMP_opt

   CHARACTER(len=3), INTENT(IN)  :: SNOWZREF_opt
!                                 reference height is constant or variable from the snow surface
!                                 SNOWZREF_opt='CST' constant reference height from the snow surface
!                                  SNOWZREF_opt='VAR' variable reference height from the snow
!                                 surface (i.e. constant  from the ground)
   REAL*8, INTENT(IN)              :: SNOWMAK_dz
!                         Snowmaking thickness (m)
   LOGICAL, INTENT(IN)             :: SNOWCOMPACT_BOOL
   LOGICAL, INTENT(IN)             :: SNOWMAK_BOOL
   LOGICAL, INTENT(IN)          :: SNOWTILLER_BOOL
   LOGICAL, INTENT(IN)          :: SELF_PROD_BOOL
   LOGICAL, INTENT(IN)          :: SNOWMAK_PROP_BOOL
   LOGICAL, INTENT(INOUT)          :: PRODSNOWMAK_BOOL
   REAL, INTENT(IN)               :: SLOPE_DIR  ! typical slope aspect in the grid
   REAL*8 :: ZP_SWNETSNOW
   REAL*8 :: ZP_SWNETSNOWS
   REAL*8 :: ZP_LWNETSNOW
   REAL*8 :: ZP_RNSNOW
   REAL*8 :: ZP_HSNOW
   REAL*8 :: ZP_GFLUXSNOW   ! it has not been initialized  in the CALL_MODEL
   REAL*8 :: ZP_HPSNOW
   REAL*8 :: ZP_LES3L
   REAL*8 :: ZP_LEL3L
   REAL*8 :: ZP_EVAP
   REAL*8 :: ZP_EVAPCOR  ! it has not been initialized  in the CALL_MODEL

! ***************************************************************************
! Local variable
! ***************************************************************************
   REAL*8 :: ZP_GSFCSNOW  !heat flux between the surface and sub-surface
!                                           snow layers (for energy budget diagnostics) (W/m2)
   REAL*8 :: SOILCOND !  SOILCOND = soil thermal conductivity [W/(m K)]


! === Only call snow model when there is snow on the surface
!              exceeding a minimum threshold OR if the equivalent
!              snow depth falling during the current time step exceeds
!              this limit.
   INTEGER   :: JWRK
   REAL       :: ZSNOW, ZSNOWFALL
!                                      ZSNOW        = snow depth (m)
!                                      ZSNOWFALL    = minimum equivalent snow depth
!                                                     for snow falling during the
!                                                     current time step (m)
   logical, intent(in) :: usemonalb  ! if usemonalb == .true., then the alb value passed to 
               !LDT will be used as the background snow-free albedo term.  
               ! if usemonalb == .false., then alb will be sett to 0.2  

   ZP_LES3L = 0.0  ! 1
   ZP_LEL3L = 0.0  ! 1
   ZP_EVAP = 0.0  ! 1
   ZP_RNSNOW = 0 ! 1
   ZP_HSNOW = 0 ! 1
   ZP_HPSNOW = 0 ! 1
   THRUFAL = 0.0 ! 1 PTHRUFAL(:)    = 0.0
   ZP_EVAPCOR = 0.0 ! 1 PEVAPCOR(:)  = 0.0
   QS =  LIS_rc%udef !XUNDEF ! 1 QS(:)  = XUNDEF
   RI_n = LIS_rc%udef !XUNDEF ! 1 PRI(:)  = XUNDEF
   EMISNOW = 0.99 ! NOTE:  snow makeing is not active, so we can use EMISNOWout(:) = 0.99
   ZP_SWNETSNOW = 0.0   ! 1    ZSWNET_N(:)       = 0.0
   ZP_SWNETSNOWS = 0.0 ! 1    ZSWNET_NS(:)     = 0.0
   ZP_LWNETSNOW = 0.0   ! 1    ZLWNET_N(:)       = 0.0

! it is not initialized in the SURFEX-Crocus. Initialized to zero here, otherwise the value in the snowcro.F90 would be -9.255963134931783E+061
   ZP_GSFCSNOW = 0

! ===========================CALL MODEL ====================================
   ZSNOW = 0.
   ZSNOWFALL = 0.0

   DO JWRK = 1, SIZE(SNOWSWE)
      ZSNOW = ZSNOW + SNOWSWE(JWRK)/SNOWRHO(JWRK)
   END DO
   ZSNOWFALL = SRSNOW*PTSTEP/XRHOSMAX_ES    ! maximum possible snowfall depth (m)
   IF (ZSNOW >= XSNOWDMIN .OR. ZSNOWFALL >= XSNOWDMIN) THEN

      CALL CALL_MODEL(1, nsnow, 1)
   ELSE

! Remove trace amounts of snow and reinitialize snow prognostic variables
! if snow cover is ablated.
! MN: TODO For now set these variables to ZERO and UNDIFF
! TODO look at snow3L_isba.F90 for details

      ZP_LES3L = 0.0
      ZP_LEL3L = 0.0
      ZP_EVAP = 0.0
      !THRUFAL = 0.0
      THRUFAL = MAX(0.0, sum(SNOWSWE)/PTSTEP + SRSNOW + RRSNOW) ! kg m-2 s-1   ! - PEVAP(:)*ZPSN(:)
      SNOWALB = LIS_rc%udef !XUNDEF
      ZP_GSFCSNOW = 0.0
      ZP_EVAPCOR = 0.0
      ZP_GFLUXSNOW = 0 ! TO DO check snow3L_isba.F90!  NOTE: it has not been initialized in the CALL_MODEL but in the outpit it is zero
      !SNOWHMASS = 0.0  ! TO DO check snow3L_isba.F90
      GRNDFLUX = 0.0
      SNOWSWE(:) = 0.0
!Prognostic variables forced to XUNDEF for correct outputs
      SNOWHEAT(:) = LIS_rc%udef ! XUNDEF !1.0E+020 ! LIS_rc%udef
      SNOWRHO(:) = LIS_rc%udef !  XUNDEF ! 10000000000000 !LIS_rc%udef
      SNOWAGE(:) = LIS_rc%udef !  XUNDEF !LIS_rc%udef

      SNOWTEMP(:) = 273.16 !LIS_rc%udef
      SNOWLIQ(:) = LIS_rc%udef !  XUNDEF ! LIS_rc%udef
      SNOWDZ(:) = 0.0
      !ZP_SNOWIMPUR (:,:,:) = 0
      SNOWGRAN1(:) = LIS_rc%udef ! XUNDEF ! LIS_rc%udef
      SNOWGRAN2(:) = LIS_rc%udef ! XUNDEF ! LIS_rc%udef
      SNOWHIST(:) = LIS_rc%udef !  XUNDEF ! LIS_rc%udef
      ! MN added
      !EMISNOW = 0.99 !
      SNOWHMASS = LIS_rc%udef !  XUNDEF ! LIS_rc%udef

   ENDIF
! ===============================================================

!================================================================
CONTAINS
!
!================================================================
   SUBROUTINE CALL_MODEL(KSIZE1, KSIZE2, KSIZE4)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: KSIZE1
      INTEGER, INTENT(IN) :: KSIZE2
      INTEGER, INTENT(IN) :: KSIZE4

      REAL*8, DIMENSION(1:KSIZE1)  :: ZP_PEW_A_COEF
      REAL*8, DIMENSION(1:KSIZE1)  :: ZP_PEW_B_COEF
      REAL*8, DIMENSION(1:KSIZE1)  :: ZP_PET_A_COEF
      REAL*8, DIMENSION(1:KSIZE1)  :: ZP_PEQ_A_COEF
      REAL*8, DIMENSION(1:KSIZE1)  :: ZP_PET_B_COEF
      REAL*8, DIMENSION(1:KSIZE1)  :: ZP_PEQ_B_COEF

      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2)  :: SNOWSWEinout
!  Z_PSNOWSWE  = Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2)  :: SNOWRHOinout
!  SNOWRHO  = Snow layer(s) averaged density (kg/m3)
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2)  :: SNOWHEATinout
!  SNOWHEAT = Snow layer(s) heat content (J/m2)
      REAL*8, DIMENSION(1:KSIZE1)  :: SNOWALBinout
!  SNOWALB = Prognostic surface snow albedo
! (does not include anything but
! the actual snow cover)
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2)  :: SNOWGRAN1inout
!  SNOWGRAN1 = Snow layers grain feature 1
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2)  :: SNOWGRAN2inout
!  SNOWGRAN2 = Snow layer grain feature 2
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2)  :: SNOWHISTinout
!  SNOWHIST  = Snow layer grain historical
! parameter (only for non
! dendritic snow)
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2)  :: SNOWAGEinout  ! Snow grain age
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2, 1:nimpur)  :: ZP_SNOWIMPURinout
!                          Snow impurity content (g) (LOCATION,LAYER,NIMPUR)) Impur type :1/BC 2/Dust
!REAL, INTENT(IN) :: PTSTEP           ! PTSTEP  = time step of the integration
      REAL*8, DIMENSION(1:KSIZE1) :: PPSin               ! ZP_PS = surface pressure
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_PAin          ! pressure at lowest atmos. level
      REAL*8, DIMENSION(1:KSIZE1) :: SRSNOWin    ! ZP_SRSNOW  = snow rate (SWE) [kg/(m2 s)]
      REAL*8, DIMENSION(1:KSIZE1) :: RRSNOWin   ! ZP_RRSNOW = rain rate [kg/(m2 s)]
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_PSN3Lin          ! ZP_PSN3L = snow fraction
      REAL*8, DIMENSION(1:KSIZE1) :: TAin               ! ZP_TA  = atmospheric temperature at level za (K)
      REAL*8, DIMENSION(1:KSIZE1) :: TGin               ! ZP_TG  = Surface soil temperature (effective
! temperature the of layer lying below snow)
      REAL*8, DIMENSION(1:KSIZE1) :: SW_RADin   ! ZP_SW_RAD = incoming solar radiation (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: QAin          !Z_PQA = atmospheric specific humidity at level za
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_VMODin   ! ZP_VMOD = modulus of the wind parallel to the orography (m/s)
      REAL*8, DIMENSION(1:KSIZE1) :: Wind_Ein
      REAL*8, DIMENSION(1:KSIZE1) :: Wind_Nin
      REAL*8, DIMENSION(1:KSIZE1) :: LW_RADin   ! ZP_LW_RAD = atmospheric infrared radiation (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_RHOAin          ! ZP_RHOA  = air density
      REAL*8, DIMENSION(1:KSIZE1) :: UREFin           ! ZP_UREF  = reference height of the wind
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_EXNSin          ! ZP_EXNS  = Exner function at surface
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_EXNAin          ! ZP_EXNA  = Exner function at lowest atmos level
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_DIRCOSZWin ! ZP_DIRCOSZW = Cosinus of the angle between the
!  normal to the surface and the vertical
      REAL*8, DIMENSION(1:KSIZE1) :: SLOPEin !    ZP_SLOPE  = Slope (degree, from LDT )
      REAL*8, DIMENSION(1:KSIZE1) :: ZREFin          ! ZP_ZREF  = reference height of the first atmospheric level
      REAL*8, DIMENSION(1:KSIZE1) :: Z0NATin           ! Z0NAT (PZ0) = grid box average roughness length
      REAL*8, DIMENSION(1:KSIZE1) :: Z0EFFin          ! Z0EFF  = roughness length for momentum
      REAL*8, DIMENSION(1:KSIZE1) :: Z0HNATin    ! Z0HNAT (PZOH)  = grid box average roughness length for heat
      REAL*8, DIMENSION(1:KSIZE1) :: ALBin           ! ALB  = soil/vegetation albedo
      REAL*8, DIMENSION(1:KSIZE1) :: SOILCONDin! SOILCOND = soil thermal conductivity [W/(m K)]
      REAL*8, DIMENSION(1:KSIZE1) :: D_Gin           ! D_G  = Assumed first soil layer thickness (m)
!  Used to calculate ground/snow heat flux
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2) :: SNOWLIQout   ! SNOWLIQ  = Snow layer(s) liquid water content (m)
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2) :: SNOWDZout           ! SNOWDZ = Snow layer(s) thickness (m)
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE2) :: SNOWTEMPinout ! SNOWTEMP = Snow layer(s) temperature (m)
      REAL*8, DIMENSION(1:KSIZE1) :: THRUFALout   ! THRUFAL  = rate that liquid water leaves snow pack:
!  paritioned into soil infiltration/runoff by ISBA [kg/(m2 s)]
      REAL*8, DIMENSION(1:KSIZE1) :: GRNDFLUXinout !REAL, DIMENSION(1:KSIZE1) :: GRNDFLUXinout
!                         GRNDFLUX = soil/snow interface heat flux (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_GSFCSNOWout ! SURFEX Feb 21 20
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_EVAPCORout   ! ZP_EVAPCOR  = evaporation/sublimation correction term:
!  extract any evaporation exceeding the
!  actual snow cover (as snow vanishes)
!  and apply it as a surface soil water
!  sink. [kg/(m2 s)]
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_GFLXCORout
!                          ZP_GFLXCOR  = flux correction to underlying soil for vanishing snowpack
!  (to put any energy excess from snow to soil) (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_SWNETSNOWout
!                         ZP_SWNETSNOW = net shortwave radiation entering top of snowpack
!  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_SWNETSNOWSout
!                         ZP_SWNETSNOWS= net shortwave radiation in uppermost layer of snowpack
!  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
! Used for surface energy budget diagnostics
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_LWNETSNOWout
!                          ZP_LWNETSNOW = net longwave radiation entering top of snowpack
!  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_RNSNOWout   ! ZP_RNSNOW = net radiative flux from snow (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_HSNOWout! REAL, DIMENSION(1:KSIZE1) :: ZP_HSNOWout   ! ZP_HSNOW  = sensible heat flux from snow (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_GFLUXSNOWout   ! ZP_GFLUXSNOW  = net heat flux from snow (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_HPSNOWout   ! ZP_HPSNOW = heat release from rainfall (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_LES3Lout           ! ZP_LES3L  = sublimation (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_LEL3Lout   ! ZP_LEL3L  = evaporation heat flux from snow (W/m2)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_EVAPout   ! ZP_EVAP = total evaporative flux (kg/m2/s)
      REAL*8, DIMENSION(1:KSIZE1) :: SNDRIFTout   ! SNDRIFT  = blowing snow sublimation (kg/m2/s)
      REAL*8, DIMENSION(1:KSIZE1) :: RI_nout          ! RI_n = Ridcharson number
      REAL*8, DIMENSION(1:KSIZE1) :: EMISNOWout   ! EMISNOW  = snow surface emissivity
      REAL*8, DIMENSION(1:KSIZE1) :: CDSNOWout   ! CDSNOW = drag coefficient for momentum over snow
      REAL*8, DIMENSION(1:KSIZE1) :: USTARSNOWout   ! USTARSNOW  = friction velocity over snow (m/s)
      REAL*8, DIMENSION(1:KSIZE1) :: CHSNOWout   ! CHSNOW = drag coefficient for heat over snow
      REAL*8, DIMENSION(1:KSIZE1) :: SNOWHMASSout
!                         SNOWHMASS  = heat content change due to mass
!  changes in snowpack (J/m2): for budget
!  calculations only.
      REAL*8, DIMENSION(1:KSIZE1) :: QSout          ! QS = surface humidity
      REAL*8, DIMENSION(1:KSIZE1) :: PERMSNOWFRACin
!                         PPERMSNOWFRAC  = fraction of permanet snow/ice
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_ZENITHin   ! solar zenith angle
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_ANGL_ILLUMin
!                         Effective illumination angle, Angle between the sun and the
!                          normal to the ground (=zenith if no slope) used in TARTES
      REAL*8, DIMENSION(1:KSIZE1) :: LATin
      REAL*8, DIMENSION(1:KSIZE1) :: LONin
      REAL*8, DIMENSION(1:KSIZE1, 1:4)   :: ZP_BLOWSNWin
!                         Properties of deposited blowing snow (from Sytron or Meso-NH/Crocus)
!                          1 : Deposition flux (kg/m2/s)
!                          2 : Density of deposited snow (kg/m3)
!                          3 : SGRA1 of deposited snow
!                          4 : SGRA2 of deposited snow


!
! see SURFEX ticket # 1599: You should prefer CSNOWRAD='B92' at the moment for the first tests.
! In that case, the direct-diffuse partition is not used, a 3 band spectral fixed repartition
!  is used from the global radiation as detailed in Vionnet et al 2012.
! 'T17' option comes with a number of other problems (impurities management, memory issues, etc.).
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE4) :: ZP_DIR_SWin  ! dimension (spectral band: 186 bands? )
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE4) :: ZP_SCA_SWin ! dimension (spectral band: 186 bands? )
!                                        direct and diffuse spectral irradiance (W/m2/um)
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE4) :: ZP_SPEC_ALBout ! dimension (spectral band: 186 bands? )
      REAL*8, DIMENSION(1:KSIZE1, 1:KSIZE4) :: ZP_DIFF_RATIOout ! dimension (spectral band: 186 bands? )
!                                        spectral albedo and diffuse to total irradiance ratio

      REAL*8, DIMENSION(1:KSIZE1, 1:nimpur) :: IMPWETin  ! dimension
      REAL*8, DIMENSION(1:KSIZE1, 1:nimpur) :: IMPDRYin  ! dimension
!                                        Dry and wet deposit coefficient from Forcing File(g/m²/s)

      REAL*8, DIMENSION(1:KSIZE1) :: SNOWMAK_dzin !Snowmaking thickness (m)
      LOGICAL, DIMENSION(1:KSIZE1) :: PRODSNOWMAK_BOOLinout
      REAL*8, DIMENSION(1) :: ZP_TIME ! current time
!REAL, DIMENSION(1:KSIZE1) :: ZPZENITH     ! Solar zenithal angle
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_AZIMSOL  ! Solar azimuthal angle
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_TSUN  ! Solar time
      REAL*8, DIMENSION(1:KSIZE1) :: SLOPE_DIRin  ! typical slope aspect in the grid
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_CDN  ! neutral drag coefficient for momentum (not used)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_AC  ! aerodynamical resistance(not used)
      REAL*8, DIMENSION(1:KSIZE1) :: ZP_RA  ! aerodynamical resistance(not used)
      REAL*8 ::PCONDDRY ! soil dry thermal conductivity (W m-1 K-1) (MN: will not be used)  
      REAL*8 ::ZCONDSLDZ  ! soil solids thermal conductivity (W m-1 K-1)
      REAL*8 ::ZQUARTZ ! local variables for thermal conductivity  
      REAL*8 ::ZGAMMAD ! local variables for thermal conductivity  
      REAL*8 ::watsat ! local variables for thermal conductivity
      REAL*8 ::tkm ! local variables for thermal conductivity
      REAL*8 ::tkmg ! local variables for thermal conductivity
      REAL*8 :: XCONDI !    = 2.22
      REAL*8 :: XCONDWTR ! = 0.57   ! W/(m K)  Water thermal conductivity
      ! local variables for thermal conductivity
      !REAL*8 :: XWGI          ! soil volumetric frozen water content (m3/m3)
      !REAL*8 :: XWG           ! soil volumetric liquid water content (m3/m3)
      REAL*8 :: XWGMIN  ! = 0.001   ! (m3 m-3)
      REAL*8 :: ZFROZEN2DF
      REAL*8 :: ZUNFROZEN2DF
      REAL*8 :: ZWORK1, ZWORK2 , ZWORK3
      REAL*8 :: ZCONDSATDF
      REAL*8 :: ZSATDEGDF
      REAL*8 :: ZKERSTENDF
      REAL*8 :: ZLOG_CONDI  ! = LOG(XCONDI)
      REAL*8 :: ZLOG_CONDWTR ! = LOG(XCONDWTR)
      REAL*8 ::  ZCONDDRYZ ! soil dry thermal conductivity (W m-1 K-1) 
      character(len=12)  :: nowdate ! the date of each time step, ( yyyymmddhhmm )
      real, external     ::  crocus81_month_d  ! external function (follows this main program):  given an array (dimension 12)
      !                                      ! representing monthly values for some parameter, return a value for 
      !                                      ! a specified date.

!_______________________________________________
! test : initialize to zero
! ______________________________________________
      ZP_PEW_A_COEF(:) = 0
      ZP_PEW_B_COEF(:) = 0
      ZP_PET_A_COEF(:) = 0
      ZP_PEQ_A_COEF(:) = 0
      ZP_PET_B_COEF(:) = 0
      ZP_PEQ_B_COEF(:) = 0
      SNOWSWEinout(:, :) = 0
      SNOWRHOinout(:, :) = 0
      SNOWHEATinout(:, :) = 0
      SNOWALBinout(:) = 0
      SNOWGRAN1inout(:, :) = 0
      SNOWGRAN2inout(:, :) = 0
      SNOWHISTinout(:, :) = 0
      SNOWAGEinout(:, :) = 0
      ZP_SNOWIMPURinout(:, :, :) = 0
      PPSin(:) = 0
      SRSNOWin(:) = 0.
      RRSNOWin(:) = 0.
      ZP_PSN3Lin(:) = 0.
      TAin(:) = 0
      TGin(:) = 0
      SW_RADin(:) = 0
      QAin(:) = 0
      ZP_VMODin(:) = 0
      Wind_Ein(:) = 0
      Wind_Nin(:) = 0
      LW_RADin(:) = 0
      ZP_RHOAin(:) = 0
      UREFin(:) = 0
      ZP_EXNSin(:) = 0
      ZP_EXNAin(:) = 0
      ZP_DIRCOSZWin(:) = 0
      SLOPEin(:) = 0
      ZREFin(:) = 0
      Z0NATin(:) = 0
      Z0EFFin(:) = 0
      Z0HNATin(:) = 0
      ALBin(:) = 0
      SOILCONDin(:) = 0
      D_Gin(:) = 0
      SNOWLIQout(:, :) = 0
      SNOWDZout(:, :) = 0
      SNOWTEMPinout(:, :) = 0
      THRUFALout(:) = 0
      ZP_GSFCSNOWout(:) = 0.0
      ZP_EVAPCORout(:) = 0
      !ZP_GFLXCORout(:) = 0
      GRNDFLUXinout(:) = 0
      ZP_SWNETSNOWout(:) = 0
      ZP_SWNETSNOWSout(:) = 0
      ZP_LWNETSNOWout(:) = 0
      ZP_RNSNOWout(:) = 0
      ZP_HSNOWout(:) = 0
      ZP_GFLUXSNOWout(:) = 0
      ZP_HPSNOWout(:) = 0
      ZP_LES3Lout(:) = 0
      ZP_LEL3Lout(:) = 0
      ZP_EVAPout(:) = 0
      SNDRIFTout(:) = 0
      RI_nout(:) = 0.2            ! snow3L_isba   --> PRI(:)  = XUNDEF  ,   it has not been initialized in the CALL_MODEL
      EMISNOWout(:) = 0.99  ! see ini_surf_csts.F90  NOTE:  snow makeing is not active, so we can use EMISNOWout(:) = 0.99
      CDSNOWout(:) = 0 ! will be computed
      USTARSNOWout(:) = 0
      CHSNOWout(:) = 0 ! will be computed
      SNOWHMASSout(:) = 0
      QSout(:) = 0
      PERMSNOWFRACin(:) = 0
      ZP_ZENITHin(:) = 0
      ZP_ANGL_ILLUMin(:) = 0
      LATin(:) = 0
      LONin(:) = 0
      ZP_DIR_SWin(:, :) = 0
      ZP_SCA_SWin(:, :) = 0
      ZP_SPEC_ALBout(:, :) = 0
      ZP_DIFF_RATIOout(:, :) = 0
      IMPWETin(:, :) = 0
      IMPDRYin(:, :) = 0
      SNOWMAK_dzin(:) = 0
      SLOPE_DIRin(:) = 0
      ZP_BLOWSNWin(:, :) = 0
      PRODSNOWMAK_BOOLinout(:) = .FALSE.
!_______________________________________________
! end test : initialize to zero
! ______________________________________________
      USTARSNOW = 0 ! it has not been initialized in the CALL_MODEL
      !CHSNOW = 0 ! it has not been initialized in the CALL_MODEL
      !SNOWHMASS = 0 ! it has not been initialized in the CALL_MODEL
      XRIMAX = 0.2 ! in surfex it is defined in the namelist

      CALL ini_csts ! routine to initialize the module MODD_CST

      CALL INI_SURF_CSTS_SUB ! routine to initialize all surface parameter
! ***************************************************************************
! Set up local variables for dimension match
! ***************************************************************************
      SNOWSWEinout(1, :) = SNOWSWE(:)
      SNOWRHOinout(1, :) = SNOWRHO(:)
      SNOWHEATinout(1, :) = SNOWHEAT(:)
      SNOWALBinout(1) = SNOWALB
      SNOWGRAN1inout(1, :) = SNOWGRAN1(:)
      SNOWGRAN2inout(1, :) = SNOWGRAN2(:)
      SNOWHISTinout(1, :) = SNOWHIST(:)
      SNOWAGEinout(1, :) = SNOWAGE(:)
      ZP_SNOWIMPURinout(1, nsnow, nimpur) = 0 ! ZP_SNOWIMPUR(:,:) ! B92--> snowcro.F90 uses
!                                                 age to compute the effect of impurity Vionnet et al 2012
      PPSin(1) = PPS
      SRSNOWin(1) = SRSNOW
      RRSNOWin(1) = RRSNOW
!ZP_PSN3Lin(1,1) = ZP_PSN3L
      TAin(1) = TA
      TGin(1) = TG
      SW_RADin(1) = SW_RAD
      QAin(1) = QA
!ZP_VMODin(1,1) = ZP_VMOD
      Wind_Ein(1) = Wind_E
      Wind_Nin(1) = Wind_N
      LW_RADin(1) = LW_RAD
!ZP_RHOAin(1,1) = ZP_RHOA
      UREFin(1) = UREF
!ZP_EXNSin(1,1) = ZP_EXNS
!ZP_EXNAin(1,1) = ZP_EXNA
!ZP_DIRCOSZWin(1,1) = ZP_DIRCOSZW
      SLOPEin(1) = SLOPE
      ZREFin(1) = ZREF
      Z0NATin(1) = Z0NAT
      Z0EFFin(1) = Z0EFF
      Z0HNATin(1) = Z0HNAT
      !ALBin(1) = 0.2 ! ALB  !  soil/vegetation albedo   set to 0.2 in the SURFEX-Crocus
      SOILCONDin(1) = SOILCOND
      D_Gin(1) = D_G
      SNOWLIQout(1, :) = SNOWLIQ
      SNOWDZout(1, :) = SNOWDZ
      SNOWTEMPinout(1, :) = SNOWTEMP
      THRUFALout(1) = THRUFAL
      ZP_GSFCSNOWout(1) = ZP_GSFCSNOW
      ZP_EVAPCORout(1) = ZP_EVAPCOR
      ZP_GFLXCORout(1) = 0 ! this is a local variable and set to zero after CALL to 'snowcor'
      GRNDFLUXinout(1) = 0 ! GRNDFLUX   , update: it is 0 in the SURFEX-Crocus
      ZP_SWNETSNOWout(1) = 0 ! ZP_SWNETSNOW        update: it is 0 in the SURFEX-Crocus
      ZP_SWNETSNOWSout(1) = 0 ! ZP_SWNETSNOWS     update: it is 0 in the SURFEX-Crocus
      ZP_LWNETSNOWout(1) = 0 ! ZP_LWNETSNOW     update: it is 0 in the SURFEX-Crocus
      ZP_RNSNOWout(1) = 0 ! ZP_RNSNOW   update: it is 0 in the SURFEX-Crocus
      ZP_HSNOWout(1) = 0 ! ZP_HSNOW   update: it is 0 in the SURFEX-Crocus
      ZP_GFLUXSNOWout(1) = 0 ! ZP_GFLUXSNOW      , update: it is 0 in the SURFEX-Crocus
      ZP_HPSNOWout(1) = 0 ! ZP_HPSNOW
      ZP_LES3Lout(1) = 0 ! ZP_LES3L    update: it is 0 in the SURFEX-Crocus
      ZP_LEL3Lout(1) = 0 !ZP_LEL3L
      ZP_EVAPout(1) = 0 ! ZP_EVAP
      SNDRIFTout(1) = 0 ! SNDRIFT  , update: it is not initialized in the scow3L_isba.F90 but, it is initialized to 0 in the snowcro.F90
      RI_nout(1) = 0 !                    , update: it is 0 in the SURFEX-Crocus
      EMISNOWout(1) = EMISNOW
      CDSNOWout(1) = 0 ! CDSNOW ! will be computed  , update: it is 0 in the SURFEX-Crocus  --> the call to subroutine   MODI_SURFACE_CD_SUB to compute CDSNOWout was commented out
      USTARSNOWout(1) = 0 ! USTARSNOW , update: it is not initialized in the scow3L_isba.F90
      CHSNOWout(1) = 0 !CHSNOW ! will be computed  , update: it is 0 in the SURFEX-Crocus  --> the call to subroutine   MODI_SURFACE_AERO_COND_SUB to compute CHSNOWout was commented out
      SNOWHMASSout(1) = 0 ! SNOWHMASS     , update: it is not initialized in the scow3L_isba.F90 but it is initialized to 0 in the snowcro.F90
      QSout(1) = 0 ! QS
      PERMSNOWFRACin(1) = PERMSNOWFRAC
!ZP_ZENITHin (1,1)                = ZP_ZENITH
!ZP_ANGL_ILLUMin (1,1)           = ZP_ANGL_ILLUM
      LATin(1) = LAT
!LATin(1)=45.17
      LONin(1) = LON

! T17 --> 186 band , B92 --> the direct-diffuse partition is not used, a 3 band spectral fixed repartition
! is used from the global radiation as detailed in Vionnet et al 2012.
      ZP_DIR_SWin(1, :) = 0 ! ZP_DIR_SW(:)
      ZP_SCA_SWin(1, :) = 0 ! ZP_SCA_SW(:)
!ZP_SPEC_ALBout (1,:) = ZP_SPEC_ALB(:) ! (spectral band: T17 --> 186 bands )
!ZP_DIFF_RATIOout(1,:) = ZP_DIFF_RATIO(:) ! dimension (spectral band: T17--> 186 bands )
!                                                 spectral albedo and diffuse to total irradiance ratio
      IMPWETin(1, :) = IMPWET(:) ! dimension (1,nimpur,1)
      IMPDRYin(1, :) = IMPDRY(:) ! dimension (1,nimpur,1)
!                                                Dry and wet deposit coefficient from Forcing File(g/m²/s)
      SNOWMAK_dzin(1) = SNOWMAK_dz
      PRODSNOWMAK_BOOLinout(1) = PRODSNOWMAK_BOOL
      SLOPE_DIRin(1) = SLOPE_DIR

      ZP_BLOWSNWin(1, 1:4) = 0 ! No snow drift; set the value to zero

  ! if the usemonalb flag is .true., we want to provide alb from the user-specified
  ! trend through the year, rather than let lsmruc calculate it for us.
  !
  write(nowdate, "(i0.4,i0.2,i0.2,i0.2,i0.2)"), year, month, day, hour, minute ! yyyymmddhhmm
  if (usemonalb .eqv. .true. ) then
     ALBin(1) = crocus81_month_d(ALB, nowdate)
  else
     ALBin(1) = 0.2 ! soil/vegetation albedo (ALB) set to 0.2 in the SURFEX-Crocus
  endif
! ***************************************************************************
! Compute variables
! ***************************************************************************

! Compute wind speed
! --------------------------------------------
      ZP_VMODin(1) = sqrt(Wind_Ein(1)*Wind_Ein(1) + Wind_Nin(1)*Wind_Nin(1))

! Compute cosine of Slope
! ---------------------------------------------------
      ZP_DIRCOSZWin(1) = COS(SLOPEin(1))


! Compute the snow fraction (is it a fraction of total precip? or fraction of grid cell covered with snow?) : I think it is fraction of grid cell
! -------------------------------------------------------
!if (SRSNOWin(1) /= 0. ) then
! ZP_PSN3Lin(1)= SRSNOWin(1) /(SRSNOWin(1) +RRSNOWin(1))
!endif
! Compute air density
! -----------------------------------------
      XRD = XAVOGADRO*XBOLTZ/XMD
      XRV = XAVOGADRO*XBOLTZ/XMV

      ZP_RHOAin(1) = PPSin(1)/(XRD*TAin(1)*(1.+((XRV/XRD) - 1.)*QAin(1)) + XG*ZREFin(1))
! specific humidity (conversion from kg/m3 to kg/kg)  ! coupling_isban.F90

! Compute the Exner functions
! -------------------------------------------------------------
      XCPD = 7.*XRD/2.
      ZP_EXNSin(1) = (PPSin(1)/XP00)**(XRD/XCPD) ! Exner function at surface

! For ZP_EXNAin we need pressure at lowest atmos. level
! Assistance SURFEX tiket #1599 --> XPA = XPS - XRHOA * XZREF * XG
      ZP_PAin(1) = PPSin(1) - ZP_RHOAin(1)*ZREFin(1)*XG
      ZP_EXNAin(1) = (ZP_PAin(1)/XP00)**(XRD/XCPD) ! Exner function at lowest atmos. level

! Compute ZP_ZENITH, ZP_AZIMSOL
! ------------------------------------------------------------------------------
      ZP_TIME = LIS_rc%hr*3600 + LIS_rc%mn*60 + LIS_rc%ss ! current time
      CALL SUNPOS(YEAR, MONTH, DAY, ZP_TIME, &
                  LONin, LATin, ZP_TSUN, ZP_ZENITHin, ZP_AZIMSOL, 1)

! Compute the effective illumination angle
! ----------------------------------------------------------------------------------
! ZP_ZENITH & ZP_AZIMSOL from SUBROUTINE SUNPOS
! SLOPE_DIR ! Note: in SURFEX is defined as the direction of S.S.O. (deg from N clockwise) 
! Here we use LDT output and ASPECT unit is in radians. So we need to edit the following equation

      ZP_ANGL_ILLUMin(1) = ACOS((COS(ZP_ZENITHin(1))*COS(ACOS(ZP_DIRCOSZWin(1)))) + &
                                (SIN(ZP_ZENITHin(1))*SIN(ACOS(ZP_DIRCOSZWin(1))*COS(ZP_AZIMSOL(1) - (SLOPE_DIRin(1)))))) ! SLOPE_DIRin(1)*XPI/180 

!  Initialize / compute implicit coefficients:
! Comes form SUBROUTINE COUPLING_ISBA_n
      ZP_PEW_A_COEF(1) = 0.
      ZP_PEW_B_COEF(1) = SQRT(Wind_Ein(1)**2 + Wind_Nin(1)**2)
! Comes form MODULE MODE_COUPLING_CANOPY
      ZP_PET_A_COEF(1) = 0.
      ZP_PET_B_COEF(1) = TAin(1)/(ZP_PAin(1)/XP00)**(XRD/XCPD)
      ZP_PEQ_A_COEF(1) = 0.
      ZP_PEQ_B_COEF(1) = QAin(1)
! Compute (Recharson number)
! ---------------------------------------------------------------
! RI (Recharson number) >0.2-0.25 (a threshold beyond which the turbulence is suppressed)
! Crocus : Ri,crit, is set to 0.2. user manual Page = 170
! I the CALL_MODEL the Ri has not been initialized and uninitalized variable passed into the snowcro.F90
! If Ri needed to be computed we can use the following subroutine
! CALL SURFACE_RI(PTS, PQSAT, PEXNS, PEXNA, PTA, PQA, &
! PZREF, PUREF, PDIRCOSZW, PVMOD, ZRI )

! Compute drag coefficient for momentum
! -----------------------------------------------------------------------------------
! see modi_surface_cd.F90
! NOTE: we need to compute RI_nout(1,1), ZREFin(1,1), UREFin(1,1), Z0EFFin(1,1), Z0HNATin(1,1)
! PUREF (UREFin) ! reference height of the wind
! PZ0EFF (Z0EFFin) ! roughness length for momentum with subgrid-scale orography
! PZ0H (Z0HNATin) ! roughness length for heat
! NOTE: MODI_SURFACE_CD_SUB.F90 need 1D variable.

! CALL MODI_SURFACE_CD_SUB(RI_nout, ZREFin, UREFin, Z0EFFin, Z0HNATin, &
! CDSNOWout , ZP_CDN)

! Compute drag coefficient for heat
! ---------------------------------------------------------------------
! see modi_surface_aero_cond.F90
! PZREF (ZREFin) ! reference height of the first atmospheric level
! ZP_AC (out) ! aerodynamical conductance
! ZP_RA (out) ! aerodynamical resistance

! CALL MODI_SURFACE_AERO_COND_SUB(RI_nout, ZREFin, UREFin, ZP_VMODin, Z0NATin ,&
! Z0HNATin, ZP_AC, ZP_RA, CHSNOWout,SNOWRES_opt ,1)

!XCPV   = 4.* XRV
!XGAMW  = (XCL - XCPV) / XRV
!XBETAW = (XLVTT/XRV) + (XGAMW * XTT)
!XALPW  = LOG(XESTT) + (XBETAW /XTT) + (XGAMW *LOG(XTT))

!XSTEFAN = ( 2.* XPI**5 / 15. ) * ( (XBOLTZ / XPLANCK)* XBOLTZ ) * (XBOLTZ/(XLIGHTSPEED*XPLANCK))**2  ! use ini_csts.F90

! ___________________________________________ Thermal conductivity _________________________________________________
! First method (thrmcondz.F90 from SURFEX-Crocus)
! Compute thermal conductivity for dry soil (NOTE: for wet soil, we need to use SURFEX-Crocus soil.F90 in which needs LSM variables)
! ---------------------------------------------------------------------
! SAND(PSANDZ)     ! soil sand fraction (-)
! POROSITY(PWSATZ)     ! soil porosity (m3 m-3)
! PCONDDRY  ! soil dry thermal conductivity     (W m-1 K-1)
! SOILCOND (PCONDSLD)  ! soil solids thermal  conductivity (W m-1 K-1)

! Part of thrmcondz.F90 from SURFEX-Crocus
ZQUARTZ  = LIS_rc%udef ! XUNDEF
ZGAMMAD   = LIS_rc%udef !XUNDEF
! SOILCOND  = XUNDEF ! A dummy argument with the INTENT(IN) attribute shall not be defined nor become undefined
PCONDDRY  = LIS_rc%udef !XUNDEF
! Quartz content estimated from sand fraction:
IF (SAND/=LIS_rc%udef) THEN  !  WHERE(SAND/=XUNDEF)
   ZQUARTZ   = 0.038 + 0.95*SAND
! Note, ZGAMMAD (soil dry density) can be supplied from obs, but
! for mesoscale modeling, we use the following approximation
! from Peters-Lidard et al. 1998:
   ZGAMMAD   = (1.0-POROSITY)*XDRYWGHT
END IF ! WHERE
! Soil solids conductivity:
IF (ZQUARTZ >  0.20 .AND. SAND/=LIS_rc%udef) THEN  ! WHERE(ZQUARTZ >  0.20 .AND. SAND/=XUNDEF)
   SOILCOND  = (XCONDQRTZ**ZQUARTZ)*                        &
                    (XCONDOTH1**(1.0-ZQUARTZ))
! END  WHERE
ELSE IF (ZQUARTZ <= 0.20 .AND. SAND/=LIS_rc%udef) THEN ! WHERE(ZQUARTZ <= 0.20 .AND. SAND/=XUNDEF)
   SOILCOND  = (XCONDQRTZ**ZQUARTZ)*                        &
                    (XCONDOTH2**(1.0-ZQUARTZ))
END IF !WHERE
! Soil dry conductivity:
IF (SAND/=LIS_rc%udef) THEN  ! WHERE(SAND/=XUNDEF)
   PCONDDRY     = (0.135*ZGAMMAD + 64.7)/                   &
                         (XDRYWGHT - 0.947*ZGAMMAD)
END IF !WHERE



ZCONDSLDZ = SOILCOND ! soil solids thermal  conductivity (W m-1 K-1)  --> use this thrid method 
ZCONDDRYZ = PCONDDRY !  

! ---------------------------------------------------------------------             
! Second method (iniTimeConst.F90 from CLM2 )
! ---------------------------------------------------------------------
! tkmg   !thermal conductivity, soil minerals  [W/m-K]  (new)
! tkdry  !thermal conductivity, dry soil       (W/m/Kelvin)
! SOILCOND (tksatu) !thermal conductivity, saturated soil [W/m-K]  (new)     
! watsat !volumetric soil water at saturation (porosity)
              watsat = 0.489 - 0.00126*SAND*100.0
              !if(SAND.eq.0.and.CLAY.eq.0) then
              !   SAND = 0.00001
              !   CLAY = 0.00001
              !endif
              tkm              = (8.80*SAND*100.0+2.92*&
                   CLAY*100.0)/(SAND*100.0+CLAY*100.0) ! W/(m K)
              !bd               = (1.-watsat)*2.7e3
              tkmg   = tkm ** (1.- watsat)
              SOILCOND = tkmg*0.57**watsat

! ---------------------------------------------------------------------
! Third method using soildif.F90 with some assumptions
! NOTE: In future when Crocus conected to a LSM we can get the soil thermal conductivity from the LSM
! or we can get soil parameters from the LSM and compute the soil thermal conductivity using following 
! equations from the soildif.F90 (SURFEX-Crocus)
! I have made two assumptions:   
! 1-  The volumetric soil water content of the snow-covered ground is 80% of POROSITY
! this assumption is based on the volumetric soil water content of the Col de Porte site for one year 
! (the observed values are between 72%-85%). 
! 2-  Soil does not freeze during the cold season. This assumption is based on the frozen and unfrozen
! part of the soil matrix in the Col de Porte site. The maximum frozen fraction is less than 0.1 and 
! the effect of the frozen thermal conductivity is negligible. 
! ---------------------------------------------------------------------
!
XCONDI    = 2.22
XCONDWTR  = 0.57   ! W/(m K)  Water thermal conductivity

!XWGI          ! soil liquid water equivalent volumetric ice content profile (m3/m3)
!XWG           ! soil volumetric water content profile (m3/m3)
XWGMIN   = 0.001   ! (m3 m-3)
!POROSITY (XWSAT)          ! porosity profile                        (m3/m3)

ZLOG_CONDI   = LOG(XCONDI)
ZLOG_CONDWTR = LOG(XCONDWTR)
! 
! For the stand-alone version, we need to provide some default values for CWGI and XWG
!      XWGI = 0.0 ! MN set to zero 
!      XWG = POROSITY * 0.8 ! MN assume volumetric soil water content of the snow covered ground is 80% of POROSITY (For Col de Porte it is between 72-85)
      ZFROZEN2DF   = XWGI/( XWGI + MAX(XWG,XWGMIN))
      ZUNFROZEN2DF = (1.0-ZFROZEN2DF)* POROSITY
!
      ZWORK1      = LOG(ZCONDSLDZ)*(1.0- POROSITY)
      ZWORK2      = ZLOG_CONDI*( POROSITY-ZUNFROZEN2DF)
      ZWORK3      = ZLOG_CONDWTR*ZUNFROZEN2DF
      ZCONDSATDF  = EXP(ZWORK1+ZWORK2+ZWORK3)
!
      ZSATDEGDF   = MAX(0.1, (XWGI+XWG)/POROSITY)
      ZSATDEGDF   = MIN(1.0,ZSATDEGDF)
      ZKERSTENDF  = LOG10(ZSATDEGDF) + 1.0
      ZKERSTENDF  = (1.0-ZFROZEN2DF)*ZKERSTENDF + ZFROZEN2DF *ZSATDEGDF
!
! Thermal conductivity of soil:
!
      SOILCOND = ZKERSTENDF*(ZCONDSATDF-ZCONDDRYZ) + ZCONDDRYZ

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

      TPTIME%TDATE%YEAR = year
      TPTIME%TDATE%MONTH = month
      TPTIME%TDATE%DAY = day
      TPTIME%TIME = hour*3600 + minute*60

      SOILCONDin = SOILCOND
      ZP_PSN3Lin(1) = 1   !TODO  assume fraction is 1  (In SURFEX-Crocus is start from zero and in several time step it became 1.  I was not able to find out where is it computed)

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! call model physics here
      CALL SNOWCRO(SNOWRES_opt, &
                   TPTIME, &
                   OMEB_BOOL, &
                   GLACIER_BOOL, &
                   HIMPLICIT_WIND_opt, &
                   ZP_PEW_A_COEF, &
                   ZP_PEW_B_COEF, &
                   ZP_PET_A_COEF, &
                   ZP_PEQ_A_COEF, &
                   ZP_PET_B_COEF, &
                   ZP_PEQ_B_COEF, &
                   SNOWSWEinout, &
                   SNOWRHOinout, &
                   SNOWHEATinout, &
                   SNOWALBinout, &
                   SNOWGRAN1inout, &
                   SNOWGRAN2inout, &
                   SNOWHISTinout, &
                   SNOWAGEinout, &
                   ZP_SNOWIMPURinout, &
                   PTSTEP, &
                   PPSin, &
                   SRSNOWin, &
                   RRSNOWin, &
                   ZP_PSN3Lin, &
                   TAin, TGin, &
                   SW_RADin, &
                   QAin, &
                   ZP_VMODin, &
                   LW_RADin, &
                   ZP_RHOAin, &
                   UREFin, &
                   ZP_EXNSin, &
                   ZP_EXNAin, &
                   ZP_DIRCOSZWin, &
                   ZREFin, &
                   Z0NATin, &
                   Z0EFFin, &
                   Z0HNATin, &
                   ALBin, &
                   SOILCONDin, &
                   D_Gin, &
                   SNOWLIQout, &
                   SNOWTEMPinout, &
                   SNOWDZout, &
                   THRUFALout, &
                   GRNDFLUXinout, &
                   ZP_EVAPCORout, &
                   ZP_GFLXCORout, &
                   ZP_SWNETSNOWout, &
                   ZP_SWNETSNOWSout, &
                   ZP_LWNETSNOWout, &
                   ZP_RNSNOWout, &
                   ZP_HSNOWout, &
                   ZP_GFLUXSNOWout, &
                   ZP_HPSNOWout, &
                   ZP_LES3Lout, &
                   ZP_LEL3Lout, &
                   ZP_EVAPout, &
                   SNDRIFTout, &
                   RI_nout, &
                   EMISNOWout, &
                   CDSNOWout, &
                   USTARSNOWout, &
                   CHSNOWout, &
                   SNOWHMASSout, &
                   QSout, &
                   PERMSNOWFRACin, &
                   ZP_ZENITHin, &
                   ZP_ANGL_ILLUMin, &
                   LATin, &
                   LONin, &
                   ZP_BLOWSNWin, &
                   SNOWDRIFT_opt, &
                   SNOWDRIFT_SUBLIM_BOOL, &
                   SNOW_ABS_ZENITH_BOOL, &
                   SNOWMETAMO_opt, &
                   SNOWRAD_opt, &
                   ATMORAD_BOOL, &
                   ZP_DIR_SWin, &
                   ZP_SCA_SWin, &
                   ZP_SPEC_ALBout, &
                   ZP_DIFF_RATIOout, &
                   ZP_GSFCSNOWout, &
                   IMPWETin, &
                   IMPDRYin, &
                   SNOWFALL_opt, &
                   SNOWCOND_opt, &
                   SNOWHOLD_opt, &
                   SNOWCOMP_opt, &
                   SNOWZREF_opt, &
                   SNOWMAK_dzin, &
                   SNOWCOMPACT_BOOL, &
                   SNOWMAK_BOOL, &
                   SNOWTILLER_BOOL, &
                   SELF_PROD_BOOL, &
                   SNOWMAK_PROP_BOOL, &
                   PRODSNOWMAK_BOOLinout, &
                   KSIZE1, KSIZE2, KSIZE4)
!
      ZP_GFLXCORout = 0.0 ! see snow3L_isba.F90


! --------------------------------------------------------------------------------------------------------------
      SNOWHEAT(:) = SNOWHEATinout(1, :)
      SNOWRHO(:) = SNOWRHOinout(1, :)
      SNOWSWE(:) = SNOWSWEinout(1, :)
      SNOWALB = SNOWALBinout(1)
      SNOWGRAN1(:) = SNOWGRAN1inout(1, :)
      SNOWGRAN2(:) = SNOWGRAN2inout(1, :)
      SNOWHIST(:) = SNOWHISTinout(1, :)
      SNOWAGE(:) = SNOWAGEinout(1, :)
      !ZP_SNOWIMPUR   = ZP_SNOWIMPURinout(1,nsnow,nimpur)
      !ZP_PSN3L                = ZP_PSN3Lin(1,1)
      !ZP_RHOA                = ZP_RHOAin(1,1)
      !SOILCOND           = SOILCONDin(1,1)
      SNOWLIQ(:) = SNOWLIQout(1, :)
      SNOWDZ(:) = SNOWDZout(1, :)
      SNOWTEMP(:) = SNOWTEMPinout(1, :)
      THRUFAL = THRUFALout(1)
      ZP_GSFCSNOW = ZP_GSFCSNOWout(1)
      ZP_EVAPCOR = ZP_EVAPCORout(1)
      ! ZP_GFLXCOR = ZP_GFLXCORout(1)  ! it is a local variable and set to zero after a CALL to 'snowcro'
      GRNDFLUX = GRNDFLUXinout(1)
      ZP_SWNETSNOW = ZP_SWNETSNOWout(1)
      ZP_SWNETSNOWS = ZP_SWNETSNOWSout(1)
      ZP_LWNETSNOW = ZP_LWNETSNOWout(1)
      ZP_RNSNOW = ZP_RNSNOWout(1)
      ZP_HSNOW = ZP_HSNOWout(1)
      ZP_GFLUXSNOW = ZP_GFLUXSNOWout(1)
      ZP_HPSNOW = ZP_HPSNOWout(1)
      ZP_LES3L = ZP_LES3Lout(1)
      ZP_LEL3L = ZP_LEL3Lout(1)
      ZP_EVAP = ZP_EVAPout(1)
      SNDRIFT = SNDRIFTout(1)
      RI_n = RI_nout(1)
      EMISNOW = EMISNOWout(1)
      CDSNOW = CDSNOWout(1)
      USTARSNOW = USTARSNOWout(1)
      CHSNOW = CHSNOWout(1)
      SNOWHMASS = SNOWHMASSout(1)
      QS = QSout(1)
      PRODSNOWMAK_BOOL = PRODSNOWMAK_BOOLinout(1)

   END SUBROUTINE CALL_MODEL

END SUBROUTINE crocus_driver


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
real function crocus81_month_d(a12, nowdate) result (nowval)
  !
  ! given a set of 12 values, taken to be valid on the fifteenth of each month (jan through dec)
  ! and a date in the form <yyyymmdd[hhmmss]> ....
  ! 
  ! return a value valid for the day given in <nowdate>, as an interpolation from the 12
  ! monthly values.
  !
  use kwm_date_utilities_crocus81
  implicit none
  real*8, dimension(12), intent(in) :: a12 ! 12 monthly values, taken to be valid on the 15th of
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

end function crocus81_month_d
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------





