module noahmp_globals

  ! Maybe most of these can be moved to a REDPRM use statement?
  use module_sf_noahlsm, only: &
       &                       SLCATS,     &
       &                       LUCATS,     &
       &                       CSOIL_DATA, & 
       &                       BB,         &
       &                       SATDK,      &
       &                       SATDW,      &
       &                       F11,        &
       &                       SATPSI,     &
       &                       QTZ,        &
       &                       DRYSMC,     &
       &                       MAXSMC,     &
       &                       REFSMC,     &
       &                       WLTSMC,     &
       &                       RSTBL,      &
       &                       RGLTBL,     &
       &                       HSTBL,      &
       &                       NROTBL,     &
       &                       TOPT_DATA,  &
       &                       RSMAX_DATA, &
       &                       ZBOT_DATA,  &
       &                       CZIL_DATA,  &
       &                       FRZK_DATA,  &
       &                       SLOPE_DATA, &
       &                       REFDK_DATA, &
       &                       REFKDT_DATA
       
  implicit none

! ==================================================================================================
!------------------------------------------------------------------------------------------!
! Physical Constants:                                                                      !
!------------------------------------------------------------------------------------------!

  REAL, PARAMETER :: GRAV   = 9.80616   !acceleration due to gravity (m/s2)
  REAL, PARAMETER :: SB     = 5.67E-08  !Stefan-Boltzmann constant (w/m2/k4)
  REAL, PARAMETER :: VKC    = 0.40      !von Karman constant
  REAL, PARAMETER :: TFRZ   = 273.16    !freezing/melting point (k)
  REAL, PARAMETER :: HSUB   = 2.8440E06 !latent heat of sublimation (j/kg)
  REAL, PARAMETER :: HVAP   = 2.5104E06 !latent heat of vaporization (j/kg)
  REAL, PARAMETER :: HFUS   = 0.3336E06 !latent heat of fusion (j/kg)
  REAL, PARAMETER :: CWAT   = 4.188E06  !specific heat capacity of water (j/m3/k)
  REAL, PARAMETER :: CICE   = 2.094E06  !specific heat capacity of ice (j/m3/k)
  REAL, PARAMETER :: CPAIR  = 1004.64   !heat capacity dry air at const pres (j/kg/k)
  REAL, PARAMETER :: TKWAT  = 0.6       !thermal conductivity of water (w/m/k)
  REAL, PARAMETER :: TKICE  = 2.2       !thermal conductivity of ice (w/m/k)
  REAL, PARAMETER :: TKAIR  = 0.023     !thermal conductivity of air (w/m/k)
  REAL, PARAMETER :: RAIR   = 287.04    !gas constant for dry air (j/kg/k)
  REAL, PARAMETER :: RW     = 461.269   !gas constant for  water vapor (j/kg/k)
  REAL, PARAMETER :: DENH2O = 1000.     !density of water (kg/m3)
  REAL, PARAMETER :: DENICE = 917.      !density of ice (kg/m3)

!------------------------------------------------------------------------------------------!
! From the VEGPARM.TBL tables, as functions of vegetation category.
!------------------------------------------------------------------------------------------!
  INTEGER :: NROOT        !rooting depth [as the number of layers] ( Assigned in REDPRM )
  REAL    :: RGL          !parameter used in radiation stress function ( Assigned in REDPRM )
  REAL    :: RSMIN        !minimum Canopy Resistance [s/m] ( Assigned in REDPRM )
  REAL    :: HS           !parameter used in vapor pressure deficit function ( Assigned in REDPRM )
  REAL    :: RSMAX        !maximum stomatal resistance ( Assigned in REDPRM )
  REAL    :: TOPT         !optimum transpiration air temperature.

!KWM  CHARACTER(LEN=256) ::  LUTYPE
!KWM  INTEGER                    :: LUCATS, BARE
!KWM  INTEGER, PARAMETER         :: NLUS=50
!KWM  INTEGER, DIMENSION(1:NLUS) :: NROTBL
!KWM  REAL,    DIMENSION(1:NLUS) :: RSTBL, RGLTBL, HSTBL
!KWM  REAL                       :: TOPT_DATA,RSMAX_DATA

! not further used in this version (niu):

!KWM  REAL,    DIMENSION(1:NLUS) ::  SNUPTBL, LAITBL,        &
!KWM                                 ALBTBL, SHDTBL, MAXALB
!KWM  REAL                       ::  CMCMAX_DATA,CFACTR_DATA,SBETA_DATA,&
!KWM                                 SALP_DATA  ,SMLOW_DATA ,SMHIGH_DATA

!KWM  REAL,    DIMENSION(NLUS)   ::  LAIMINTBL    !KWM
!KWM  REAL,    DIMENSION(NLUS)   ::  LAIMAXTBL    !KWM
!KWM  REAL,    DIMENSION(NLUS)   ::  EMISSMINTBL  !KWM
!KWM  REAL,    DIMENSION(NLUS)   ::  EMISSMAXTBL  !KWM
!KWM  REAL,    DIMENSION(NLUS)   ::  ALBEDOMINTBL !KWM
!KWM  REAL,    DIMENSION(NLUS)   ::  ALBEDOMAXTBL !KWM
!KWM  REAL,    DIMENSION(NLUS)   ::  Z0MINTBL     !KWM
!KWM  REAL,    DIMENSION(NLUS)   ::  Z0MAXTBL     !KWM


!------------------------------------------------------------------------------------------!
! From the SOILPARM.TBL tables, as functions of soil category.
!------------------------------------------------------------------------------------------!
  REAL    :: BEXP         !B parameter ( Assigned in REDPRM )
  REAL    :: SMCDRY       !dry soil moisture threshold where direct evap from top
                          !layer ends (volumetric) ( Assigned in REDPRM )
  REAL    :: F1           !soil thermal diffusivity/conductivity coef ( Assigned in REDPRM )
  REAL    :: SMCMAX       !porosity, saturated value of soil moisture (volumetric)
  REAL    :: SMCREF       !reference soil moisture (field capacity) (volumetric) ( Assigned in REDPRM )
  REAL    :: PSISAT       !saturated soil matric potential ( Assigned in REDPRM )
  REAL    :: DKSAT        !saturated soil hydraulic conductivity ( Assigned in REDPRM )
  REAL    :: DWSAT        !saturated soil hydraulic diffusivity ( Assigned in REDPRM )
  REAL    :: SMCWLT       !wilting point soil moisture (volumetric) ( Assigned in REDPRM )
  REAL    :: QUARTZ       !soil quartz content ( Assigned in REDPRM )

!KWM  CHARACTER*4 SLTYPE
!KWM  INTEGER                     :: SLCATS
!KWM  INTEGER, PARAMETER          :: NSLTYPE=30
!KWM  REAL, DIMENSION (1:NSLTYPE) :: BB,DRYSMC,F11,                           &
!KWM        MAXSMC, REFSMC,SATPSI,SATDK,SATDW, WLTSMC,QTZ
!------------------------------------------------------------------------------------------!
! From the GENPARM.TBL file
!------------------------------------------------------------------------------------------!
  REAL    :: SLOPE       !slope index (0 - 1) ( Assigned in REDPRM )
  REAL    :: CSOIL       !vol. soil heat capacity [j/m3/K] ( Assigned in REDPRM )
  REAL    :: ZBOT        !Depth (m) of lower boundary soil temperature ( Assigned in REDPRM )
  REAL    :: CZIL        !Calculate roughness length of heat ( Assigned in REDPRM )

  REAL    :: KDT         !used in compute maximum infiltration rate (in INFIL) ( Assigned in REDPRM )
  REAL    :: FRZX        !used in compute maximum infiltration rate (in INFIL) ( Assigned in REDPRM )

! LSM GENERAL PARAMETERS

!KWM  INTEGER :: SLPCATS
!KWM  INTEGER, PARAMETER :: NSLOPE=30
!KWM  REAL, DIMENSION (1:NSLOPE) :: SLOPE_DATA
!KWM  REAL ::  FXEXP_DATA,CSOIL_DATA,REFDK_DATA ,         &
!KWM           REFKDT_DATA,FRZK_DATA ,ZBOT_DATA ,CZIL_DATA

! =====================================options for different schemes================================
! options for dynamic vegetation: 
! 1 -> off (use table LAI; use FVEG = SHDFAC from input)
! 2 -> on (together with OPT_CRS = 1)
! 3 -> off (use table LAI; calculate FVEG)
! 4 -> off (use table LAI; use maximum vegetation fraction)

  INTEGER :: DVEG    != 2   !

! options for canopy stomatal resistance
! 1-> Ball-Berry; 2->Jarvis

  INTEGER :: OPT_CRS != 1    !(must 1 when DVEG = 2)

! options for soil moisture factor for stomatal resistance
! 1-> Noah (soil moisture) 
! 2-> CLM  (matric potential)
! 3-> SSiB (matric potential)

  INTEGER :: OPT_BTR != 1    !(suggested 1)

! options for runoff and groundwater
! 1 -> TOPMODEL with groundwater (Niu et al. 2007 JGR) ;
! 2 -> TOPMODEL with an equilibrium water table (Niu et al. 2005 JGR) ;
! 3 -> original surface and subsurface runoff (free drainage)
! 4 -> BATS surface and subsurface runoff (free drainage)

  INTEGER :: OPT_RUN != 1    !(suggested 1)

! options for surface layer drag coeff (CH & CM)
! 1->M-O ; 2->original Noah (Chen97); 3->MYJ consistent; 4->YSU consistent. 

  INTEGER :: OPT_SFC != 1    !(1 or 2 or 3 or 4)

! options for supercooled liquid water (or ice fraction)
! 1-> no iteration (Niu and Yang, 2006 JHM); 2: Koren's iteration 

  INTEGER :: OPT_FRZ != 1    !(1 or 2)

! options for frozen soil permeability
! 1 -> linear effects, more permeable (Niu and Yang, 2006, JHM)
! 2 -> nonlinear effects, less permeable (old)

  INTEGER :: OPT_INF != 1    !(suggested 1)

! options for radiation transfer
! 1 -> modified two-stream (gap = F(solar angle, 3D structure ...)<1-FVEG)
! 2 -> two-stream applied to grid-cell (gap = 0)
! 3 -> two-stream applied to vegetated fraction (gap=1-FVEG)

  INTEGER :: OPT_RAD != 1    !(suggested 1)

! options for ground snow surface albedo
! 1-> BATS; 2 -> CLASS

  INTEGER :: OPT_ALB != 2    !(suggested 2)

! options for partitioning  precipitation into rainfall & snowfall
! 1 -> Jordan (1991); 2 -> BATS: when SFCTMP<TFRZ+2.2 ; 3-> SFCTMP<TFRZ

  INTEGER :: OPT_SNF != 1    !(suggested 1)

! options for lower boundary condition of soil temperature
! 1 -> zero heat flux from bottom (ZBOT and TBOT not used)
! 2 -> TBOT at ZBOT (8m) read from a file (original Noah)

  INTEGER :: OPT_TBOT != 2   !(suggested 2)

! options for snow/soil temperature time scheme (only layer 1)
! 1 -> semi-implicit; 2 -> full implicit (original Noah)

  INTEGER :: OPT_STC != 1    !(suggested 1)
! ==================================================================================================
! runoff parameters used for SIMTOP and SIMGM:
  REAL, PARAMETER :: TIMEAN = 10.5   !gridcell mean topgraphic index (global mean)
  REAL, PARAMETER :: FSATMX = 0.38   !maximum surface saturated fraction (global mean)

! adjustable parameters for snow processes

  REAL, PARAMETER :: M      = 1.0 ! 2.50   !melting factor (-) 
  REAL, PARAMETER :: Z0SNO  = 0.002  !snow surface roughness length (m) (0.002)
  REAL, PARAMETER :: SSI    = 0.03   !liquid water holding capacity for snowpack (m3/m3) (0.03)
  REAL, PARAMETER :: SWEMX  = 1.00   !new snow mass to fully cover old snow (mm)
                                     !equivalent to 10mm depth (density = 100 kg/m3)

! NOTES: things to add or improve
! 1. lake model: explicit representation of lake water storage, sunlight through lake
!    with different purity, turbulent mixing of surface laker water, snow on frozen lake, etc.
! 2. shallow snow wihtout a layer: melting energy
! 3. urban model to be added.
! 4. irrigation
!------------------------------------------------------------------------------------------!
END MODULE NOAHMP_GLOBALS
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
#if 0 
MODULE NOAHMP_VEG_PARAMETERS

    IMPLICIT NONE

    INTEGER, PARAMETER :: MAX_VEG_PARAMS = 33
    INTEGER, PARAMETER :: MVT   = 27
    INTEGER, PARAMETER :: MBAND = 2

    INTEGER, PRIVATE :: ISURBAN
    INTEGER :: ISWATER
    INTEGER :: ISBARREN
    INTEGER :: ISSNOW
    INTEGER :: EBLFOREST

    REAL :: CH2OP(MVT)       !maximum intercepted h2o per unit lai+sai (mm)
    REAL :: DLEAF(MVT)       !characteristic leaf dimension (m)
    REAL :: Z0MVT(MVT)       !momentum roughness length (m)
    REAL :: HVT(MVT)         !top of canopy (m)
    REAL :: HVB(MVT)         !bottom of canopy (m)
    REAL :: DEN(MVT)         !tree density (no. of trunks per m2)
    REAL :: RC(MVT)          !tree crown radius (m)
    REAL :: SAIM(MVT,12)     !monthly stem area index, one-sided
    REAL :: LAIM(MVT,12)     !monthly leaf area index, one-sided
    REAL :: SLA(MVT)         !single-side leaf area per Kg [m2/kg]
    REAL :: DILEFC(MVT)      !coeficient for leaf stress death [1/s]
    REAL :: DILEFW(MVT)      !coeficient for leaf stress death [1/s]
    REAL :: FRAGR(MVT)       !fraction of growth respiration  !original was 0.3 
    REAL :: LTOVRC(MVT)      !leaf turnover [1/s]

    REAL :: C3PSN(MVT)       !photosynthetic pathway: 0. = c4, 1. = c3
    REAL :: KC25(MVT)        !co2 michaelis-menten constant at 25c (pa)
    REAL :: AKC(MVT)         !q10 for kc25
    REAL :: KO25(MVT)        !o2 michaelis-menten constant at 25c (pa)
    REAL :: AKO(MVT)         !q10 for ko25
    REAL :: VCMX25(MVT)      !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    REAL :: AVCMX(MVT)       !q10 for vcmx25
    REAL :: BP(MVT)          !minimum leaf conductance (umol/m**2/s)
    REAL :: MP(MVT)          !slope of conductance-to-photosynthesis relationship
    REAL :: QE25(MVT)        !quantum efficiency at 25c (umol co2 / umol photon)
    REAL :: AQE(MVT)         !q10 for qe25
    REAL :: RMF25(MVT)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
    REAL :: RMS25(MVT)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
    REAL :: RMR25(MVT)       !root maintenance respiration at 25c (umol co2/kg bio/s)
    REAL :: ARM(MVT)         !q10 for maintenance respiration
    REAL :: FOLNMX(MVT)      !foliage nitrogen concentration when f(n)=1 (%)
    REAL :: TMIN(MVT)        !minimum temperature for photosynthesis (k)

    REAL :: XL(MVT)          !leaf/stem orientation index
    REAL :: RHOL(MVT,MBAND)  !leaf reflectance: 1=vis, 2=nir
    REAL :: RHOS(MVT,MBAND)  !stem reflectance: 1=vis, 2=nir
    REAL :: TAUL(MVT,MBAND)  !leaf transmittance: 1=vis, 2=nir
    REAL :: TAUS(MVT,MBAND)  !stem transmittance: 1=vis, 2=nir

    REAL :: MRP(MVT)         !microbial respiration parameter (umol co2 /kg c/ s)
    REAL :: CWPVT(MVT)       !empirical canopy wind parameter

    REAL :: WRRAT(MVT)       !wood to non-wood ratio
    REAL :: WDPOOL(MVT)      !wood pool (switch 1 or 0) depending on woody or not [-]
    REAL :: TDLEF(MVT)       !characteristic T for leaf freezing [K]

    INTEGER :: IK,IM
    REAL :: TMP10(MVT*MBAND)
    REAL :: TMP11(MVT*MBAND)
    REAL :: TMP12(MVT*MBAND)
    REAL :: TMP13(MVT*MBAND)
    REAL :: TMP14(MVT*12)
    REAL :: TMP15(MVT*12)
    REAL :: TMP16(MVT*5)

    real slarea(MVT)
    real eps(MVT,5)

CONTAINS
  subroutine read_mp_veg_parameters(FILENAME_VEGTABLE,DATASET_IDENTIFIER)
    implicit none
    character(len=*), intent(in) :: FILENAME_VEGTABLE
    character(len=*), intent(in) :: DATASET_IDENTIFIER
    integer :: ierr

    ! Temporary arrays used in reshaping namelist arrays
    REAL :: TMP10(MVT*MBAND)
    REAL :: TMP11(MVT*MBAND)
    REAL :: TMP12(MVT*MBAND)
    REAL :: TMP13(MVT*MBAND)
    REAL :: TMP14(MVT*12)
    REAL :: TMP15(MVT*12)
    REAL :: TMP16(MVT*5)

    integer :: NVEG
    character(len=256) :: VEG_DATASET_DESCRIPTION

    NAMELIST / noah_mp_usgs_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_usgs_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_modis_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_modis_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    ! MPC change: enable use of alternative veg tables
    !  - in this case, tables using attributes from other models used in the PLUMBER experiment

    NAMELIST / noah_mp_plumberCABLE_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberCABLE_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_plumberCHTESSEL_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberCHTESSEL_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_plumberSUMMA_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberSUMMA_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
    CH2OP  = -1.E36
    DLEAF  = -1.E36
    Z0MVT  = -1.E36
    HVT    = -1.E36
    HVB    = -1.E36
    DEN    = -1.E36
    RC     = -1.E36
    RHOL   = -1.E36
    RHOS   = -1.E36
    TAUL   = -1.E36
    TAUS   = -1.E36
    XL     = -1.E36
    CWPVT  = -1.E36
    C3PSN  = -1.E36
    KC25   = -1.E36
    AKC    = -1.E36
    KO25   = -1.E36
    AKO    = -1.E36
    AVCMX  = -1.E36
    AQE    = -1.E36
    LTOVRC = -1.E36
    DILEFC = -1.E36
    DILEFW = -1.E36
    RMF25  = -1.E36
    SLA    = -1.E36
    FRAGR  = -1.E36
    TMIN   = -1.E36
    VCMX25 = -1.E36
    TDLEF  = -1.E36
    BP     = -1.E36
    MP     = -1.E36
    QE25   = -1.E36
    RMS25  = -1.E36
    RMR25  = -1.E36
    ARM    = -1.E36
    FOLNMX = -1.E36
    WDPOOL = -1.E36
    WRRAT  = -1.E36
    MRP    = -1.E36
    SAIM   = -1.E36
    LAIM   = -1.E36
    SLAREA = -1.E36
    EPS    = -1.E36

    open(15, file=trim(FILENAME_VEGTABLE), status='old', form='formatted', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,'("****** Error ******************************************************")')
       write(*,'("Cannot find file MPTABLE.TBL")')
       write(*,'("STOP")')
       write(*,'("*******************************************************************")')
       call wrf_error_fatal("STOP in Noah-MP read_mp_veg_parameters")
    endif

    if ( trim(DATASET_IDENTIFIER) == "USGS" ) then
       read(15,noah_mp_usgs_veg_categories)
       read(15,noah_mp_usgs_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "MODIFIED_IGBP_MODIS_NOAH" ) then
       read(15,noah_mp_modis_veg_categories)
       read(15,noah_mp_modis_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberCABLE" ) then
       read(15,noah_mp_plumberCABLE_veg_categories)
       read(15,noah_mp_plumberCABLE_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberCHTESSEL" ) then
       read(15,noah_mp_plumberCHTESSEL_veg_categories)
       read(15,noah_mp_plumberCHTESSEL_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberSUMMA" ) then
       read(15,noah_mp_plumberSUMMA_veg_categories)
       read(15,noah_mp_plumberSUMMA_parameters)
    else
       write(*,'("Unrecognized DATASET_IDENTIFIER in subroutine READ_MP_VEG_PARAMETERS")')
       write(*,'("DATASET_IDENTIFIER = ''", A, "''")') trim(DATASET_IDENTIFIER)
       call wrf_error_fatal("STOP in Noah-MP read_mp_veg_parameters")
    endif
    close(15)

    ! Problem.  Namelist reading of 2-d arrays doesn't work well when the arrays are declared with larger dimension than the
    ! variables in the provided namelist.  So we need to reshape the 2-d arrays after we've read them.

    if ( MVT > NVEG ) then

       ! 
       ! Reshape the 2-d arrays:
       ! 

       TMP10 = reshape( RHOL, (/ MVT*size(RHOL,2) /))
       TMP11 = reshape( RHOS, (/ MVT*size(RHOS,2) /))
       TMP12 = reshape( TAUL, (/ MVT*size(TAUL,2) /))
       TMP13 = reshape( TAUS, (/ MVT*size(TAUS,2) /))
       TMP14 = reshape( SAIM, (/ MVT*size(SAIM,2) /))
       TMP15 = reshape( LAIM, (/ MVT*size(LAIM,2) /))
       TMP16 = reshape( EPS,  (/ MVT*size(EPS ,2) /))

       RHOL(1:NVEG,:) = reshape( TMP10, (/ NVEG, size(RHOL,2) /))
       RHOS(1:NVEG,:) = reshape( TMP11, (/ NVEG, size(RHOS,2) /))
       TAUL(1:NVEG,:) = reshape( TMP12, (/ NVEG, size(TAUL,2) /))
       TAUS(1:NVEG,:) = reshape( TMP13, (/ NVEG, size(TAUS,2) /))
       SAIM(1:NVEG,:) = reshape( TMP14, (/ NVEG, size(SAIM,2) /))
       LAIM(1:NVEG,:) = reshape( TMP15, (/ NVEG, size(LAIM,2) /))
       EPS(1:NVEG,:)  = reshape( TMP16, (/ NVEG, size(EPS,2)  /))

       RHOL(NVEG+1:MVT,:) = -1.E36
       RHOS(NVEG+1:MVT,:) = -1.E36
       TAUL(NVEG+1:MVT,:) = -1.E36
       TAUS(NVEG+1:MVT,:) = -1.E36
       SAIM(NVEG+1:MVT,:) = -1.E36
       LAIM(NVEG+1:MVT,:) = -1.E36
       EPS( NVEG+1:MVT,:) = -1.E36
    endif

  end subroutine read_mp_veg_parameters

END MODULE NOAHMP_VEG_PARAMETERS
#endif
! ==================================================================================================
! ==================================================================================================
#if 0 
MODULE NOAHMP_RAD_PARAMETERS

    IMPLICIT NONE
 
    INTEGER I                ! loop index
    INTEGER, PARAMETER :: MSC   = 9
    INTEGER, PARAMETER :: MBAND = 2

    REAL :: ALBSAT(MSC,MBAND)   !saturated soil albedos: 1=vis, 2=nir
    REAL :: ALBDRY(MSC,MBAND)   !dry soil albedos: 1=vis, 2=nir
    REAL :: ALBICE(MBAND)       !albedo land ice: 1=vis, 2=nir
    REAL :: ALBLAK(MBAND)       !albedo frozen lakes: 1=vis, 2=nir
    REAL :: OMEGAS(MBAND)       !two-stream parameter omega for snow
    REAL :: BETADS              !two-stream parameter betad for snow
    REAL :: BETAIS              !two-stream parameter betad for snow
    REAL :: EG(2)               !emissivity

! saturated soil albedos: 1=vis, 2=nir
    DATA(ALBSAT(I,1),I=1,8)/0.15,0.11,0.10,0.09,0.08,0.07,0.06,0.05/
    DATA(ALBSAT(I,2),I=1,8)/0.30,0.22,0.20,0.18,0.16,0.14,0.12,0.10/

! dry soil albedos: 1=vis, 2=nir
    DATA(ALBDRY(I,1),I=1,8)/0.27,0.22,0.20,0.18,0.16,0.14,0.12,0.10/
    DATA(ALBDRY(I,2),I=1,8)/0.54,0.44,0.40,0.36,0.32,0.28,0.24,0.20/

! albedo land ice: 1=vis, 2=nir
    DATA (ALBICE(I),I=1,MBAND) /0.80, 0.55/

! albedo frozen lakes: 1=vis, 2=nir
    DATA (ALBLAK(I),I=1,MBAND) /0.60, 0.40/

! omega,betad,betai for snow
    DATA (OMEGAS(I),I=1,MBAND) /0.8, 0.4/
    DATA BETADS, BETAIS /0.5, 0.5/

! emissivity ground surface    
      DATA EG /0.97, 0.98/ ! 1-soil;2-lake

END MODULE NOAHMP_RAD_PARAMETERS
#endif
! ==================================================================================================


! ==================================================================================================

!MODULE MODULE_SF_NOAHMPLSM
!
!  USE NOAHMP_ROUTINES
!  USE NOAHMP_GLOBALS
!  USE NOAHMP_VEG_PARAMETERS
!
!END MODULE MODULE_SF_NOAHMPLSM
