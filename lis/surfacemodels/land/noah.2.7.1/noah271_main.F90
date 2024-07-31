!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: noah271_main
! \label{noah271_main}
! 
! !REVISION HISTORY:
!  27 Oct 2010: David Mocko, changes for Noah2.7.1 in LIS6.1
!  17 Jan 2011: David Mocko, fixes to output diagnostics for energy
!                            & water balance and fix to conversion of
!                            U,V forcing height to same level as T,q
!  27 Oct 2011: Chris Franks, added output of near sfc RH as well as
!                             heat & momemtum exchange coefficients
!
! !INTERFACE:
subroutine noah271_main(n)
! !USES:
  use LIS_coreMod
  use LIS_vegDataMod,    only : LIS_gfrac
  use LIS_timeMgrMod,    only : LIS_isAlarmRinging
  use LIS_albedoMod,     only : LIS_alb
  use LIS_constantsMod,  only : LIS_CONST_RHOFW, LIS_CONST_TKFRZ, &
                                LIS_CONST_LATVAP
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use LIS_histDataMod
  use module_sf_noah271lsm, only : SFLX, SFCDIF, SNFRAC
  use noah271_lsmMod 
  use LIS_tbotAdjustMod, only: LIS_tbotTimeUtil,LIS_updateTbot

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!  This is the entry point for calling the Noah2.7.1 LSM physics.
!  This routine calls the {\tt SFLX} routine that performs the
!  land surface computations, to solve for water and energy equations.
!  For documentation of the {\tt SFLX} and other Noah2.7.1 routines,
!  please see the "NOAH LSM USERS GUIDE" maintained at the public
!  server at NOAA NCEP:
!  http://www.emc.ncep.noaa.gov/mmb/gcp/noahlsm/Noah\_LSM\_USERGUIDE\_2.7.1.htm
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  INTEGER, PARAMETER:: NSOLD=20   ! Maximum number of soil layers
  REAL,    PARAMETER :: SALP = 4.0
  INTEGER tid
  INTEGER K          ! Loop integer
  
  INTEGER ICE         ! flag to invoke sea-ice physics (0=land)
  INTEGER NROOT       ! Number of root layers, dep. on veg. type 
  INTEGER NSOIL       ! Number of soil layers (used in noah271_main.f)
  INTEGER SOILTYP     ! Soil type (integer index)
  INTEGER T         ! Current tile NOAH model is called on, can be removed
  REAL SLOPE          ! Slope Estimate of linear reservoir coefficient  
  
  REAL BETA           ! Ratio of actual/potential EVAP (dimensionless) 
  REAL DRIP           ! Excess canopy moisture (m)
  REAL EC             ! Canopy evaporation (W m-2)
  REAL EDIR           ! Direct soil evaporation (W m-2)
  REAL ET(NSOLD)      ! Plant transp. from each root/soil layer(W m-2)
  REAL ETT            ! Accum plant transpiration (W m-2)
  REAL ESNOW          ! Sublimation from snowpack (kg m-2 s-1) (or dep.)
  REAL FLX1           ! Precip-snow Surface  (W M-2)
  REAL FLX2           ! Freezing rain latent heat flux (W M-2)
  REAL FLX3           ! Phase-change heat flux from snowmelt (W M-2)
  REAL DEW            ! Dewfall amount  (m s-1)
!  REAL RUNOFF1        ! Ground surface runoff (m s-1)
!  REAL RUNOFF2        ! Underground runoff (m s-1)
  REAL RUNOFF3        ! Runoff within soil layers (m s-1)
  REAL Q1         
  REAL ALB            ! Quarterly Snow-free albedo (fraction)
  REAL ALBEDO         ! Total Sfc. albedo fraction (incl. snow-max alb)
  REAL DQSDT          ! Slope of sat. specific hum. curve for Penman
  REAL REFKDT         ! Scalar surface runoff parameter 
  REAL REFDK          ! A reference value for KDT
  REAL KDT            ! Based on REFKDT, DKSAT, and REFDK 
  REAL FRZK           ! Ice content threshold in soil 
  REAL FRZX           ! Adjust FRZK Parameter to actual soil type 
  REAL FRZFACT        ! Used in ice content in soil 
  REAL FFROZP         ! Flag used for precipitation type  (NEW)
  REAL DQSDT2         ! Slope of sat. specific hum. curve for Penman
  REAL DT             ! time step length (sec) of physical integration
  REAL ESAT           ! Saturation vapor pressure for water (Pa)
  REAL E              ! Saturation vapor pressure
  
  REAL ETA            ! Actual latent heat flux (W m-2) (NEG: up from sfc)
  REAL ETA_KINEMATIC
  REAL EVP            ! Total Evaporation (kg m-2 s-1) (=EC+ETT+EDIR)
  REAL ETP            ! Final Potential Evapotransp. (W m-2)
  REAL FUP            ! Upward ground LW radiation (W m-2)
  REAL SHTFLX         ! Sensible heat flux (W m-2)
  REAL SOLDN          ! Solar downward radiation (W m-2; positive,
                      !   not net shortwave)
  REAL SOLNET         !  (NEW) SOLDN*(1.0-ALBEDO) - INPUT to SFLX
  REAL LWDN           ! Downward Longwave radiation (W m-2)
  REAL PRCP           ! Precipitation rate conversion (kg m-2 s-1)
  REAL CPCP           ! Convective Precipitation (kg m-2 s-1)
  REAL PTU            ! Phota Thermal Unit (dep. on air temp and rad)
  
  REAL Q2             ! Mixing ratio at 1st middle level above skin
  REAL Q2SAT          ! Sat. Mix. ratio at 1st middle level above skin
  REAL Q2SATI
  REAL SFCSPD         ! Wind speed, sqrt(u*u+v*v)
  REAL UWIND          ! U-Wind component (m s-1)
  REAL VWIND          ! V-Wind component (m s-1)
  REAL SFCPRS         ! Surface pressure (Pascals)
  REAL SFCTMP         ! Surface temperature (K)
  REAL SHDFAC         ! 12-month green vegetation fraction
  REAL SHDMIN         ! Fixed minimum green veg fraction

!      REAL MLAI(13)       ! Monthly leaf area index values

  REAL SNOMLT         ! Snow melt (m) (water equivalent)
  REAL SNOALB         ! Maximum albedo expected over deep snow
  REAL MSTAVRZ        ! Avail root zone soil moisture (unitless fraction)
  REAL MSTAVTOT       ! Total Column Soil Moisture Availability (%)
  REAL SOILM          ! Total soil column water content (m->kg m-2)
  REAL SOILRZ         ! Root zone soil column water content (kg m-2)
  REAL SOIL1M         ! Top 1-m soil column water content (kg m-2)
  REAL GFLX           ! Soil heat flux (W m-2)
!      REAL T1             ! Initial skin temperature (K) 
  REAL T14            ! Ground sfc. temp. to the 4th power (K^+4)
  REAL T1V            ! Virtual temperature at ground (sub 1)
  REAL T2V            ! Virtual temp. at 1st mid. lev. above grnd (2)
! LW    REAL TAK            ! Air temperature in K
  REAL TBOT           ! Annually-fixed, soil-temp condition at ZBOT
  REAL TH2            ! Potential temperature (K) at hgt z above grnd
  REAL TH2V           ! Virtual potential temperature (K) at z
  REAL Z              ! Height (meters) above ground in atmos. forcing
  REAL ZLVL           ! Height (meters) above ground in atmos. forcing
  REAL Z0             ! Roughness length parameter
  REAL ZFAC           ! Factor to convert U,V to same level as T,q  
  
  REAL, allocatable :: SLDPTH(:)  ! Thickness values for NSOIL layers (meters)
  REAL SNCOVR         ! FRACTIONAL SNOW COVER (UNITLESS FRACTION, 0-1)
  
  REAL RSMIN          ! MINIMUM CANOPY RESISTANCE (S M-1)
  REAL RGL            ! FROM SOLAR RAD TERM OF CANOPY RESISTANCE FUNCTION
  REAL HS             ! USED IN VAPOR PRESS. DEF. TERM OF CAN. RES. FUNCTION
  REAL SNUP           ! THRESHOLD SNOW DEPTH (IN WATER EQUIVALENT M) THAT
  !  IMPLIES 100% SNOW COVER
  REAL XLAI           ! LEAF AREA INDEX (DIMENSIONLESS)
  REAL RC             ! CANOPY RESISTANCE (S M-1)
  REAL PC             ! PLANT COEFFICIENT (UNITLESS FRACTION, 0-1) WHERE
  !   PC*ETP = ACTUAL TRANSP
  REAL RCS            ! INCOMING SOLAR RC FACTOR (DIMENSIONLESS)
  REAL RCT            ! AIR TEMPERATURE RC FACTOR (DIMENSIONLESS)
  REAL RCQ            ! ATMOS VAPOR PRESSURE DEFICIT RC FACTOR (DIMENSIONLESS)
  REAL RCSOIL         ! SOIL MOISTURE RC FACTOR (DIMENSIONLESS)
  
  REAL SMCMAX         ! MAX SOIL MOISTURE CONTENT (POROSITY)
  REAL SMCREF         ! REFERENCE SOIL MOISTURE (ONSET OF SOIL MOISTURE
  !  STRESS IN TRANSPIRATION) - Volumetric
  REAL SMCREF1
  REAL SMCWLT         ! WILTING Point (VOLUMETRIC)
  REAL SMCWLT1
  REAL SMCDRY         ! DRY SOIL MOIST THRESHOLD - DIRECT EVAP 
  !  FRM TOP LAYER ENDS (VOLUMETRIC)
  REAL PSISAT         ! SATURATED SOIL POTENTIAL
  REAL DKSAT          ! SATURATED SOIL HYDRAULIC CONDUCTIVITY
  REAL BEXP           ! THE 'B' PARAMETER
  REAL DWSAT          ! SATURATED SOIL DIFFUSIVITY
  REAL F1             ! USED TO COMPUTE SOIL DIFFUSIVITY/CONDUCTIVITY
  REAL QUARTZ         ! SOIL QUARTZ CONTENT
  REAL SMLOW          ! Soil moisture wilt and reference parameter
  REAL SMHIGH         ! Soil moisture wilt and reference parameter
  REAL TSOIL          ! SOIL SURFACE TEMPERATURE (K)
  REAL ETANRG
  REAL CH, CM
  REAL SFCTSNO
  REAL E2SAT
  REAL EMISSI_OUT

  REAL, PARAMETER:: R = 287.04
  REAL, PARAMETER:: CP = 1004.5
!  REAL, PARAMETER:: CH2O = 4.2E6
!  REAL, PARAMETER:: CSOIL = 2.00E+6
!  REAL, PARAMETER:: CAIR = 1004.0
!  REAL, PARAMETER:: CICE = 2.106E6
  REAL, PARAMETER:: LVH2O = 2.501000E+6 ! Latent heat for evapo for water  
  REAL, PARAMETER:: EPS = 0.622 ! Water/(dry air) molec mass ratio (epsilon)
  REAL, PARAMETER  :: TRESH=.95E0, A2=17.67,A3=273.15,A4=29.65,   &
       T0=273.16E0, ELWV=2.50E6,  A23M4=A2*(A3-A4)
  !--  PARAMETER USED TO CALCULATE ROUGHNESS LENGTH OF HEAT.
  !      PARAMETER(CZIL = 0.2)
  ! - Changed for version 2.6 June 2nd 2003 *
  !      PARAMETER(CZIL = 0.075)
  ! - Changed for version 2.7.1 Jesse 20041225
  !      PARAMETER(CZIL = 0.1)
  ! - CZIL now set in soil parameter table
  REAL :: roottemp
  integer :: i
  real    :: ustar, vp, rho
  real :: temp
  real :: soilhtc, soilmtc
  real :: startht, startsm, startswe, startint
  logical             :: alarmCheck
  integer             :: iret
  real :: julian_in
  real :: bdsno
  integer :: yr

  real, parameter :: CONST_LATSUB=2.83E+6 ! to be consistent with Noah SFLX

  ! SY: Begin for FLDAS
  REAL :: WRSI_TimeStep
  REAL :: WR_TimeStep
  REAL :: AET_TimeStep
  real :: soilrzmax      
  REAL, parameter  :: SUMWR_EXITWRSIBELOWTHISQTY = 0
  REAL, parameter  :: SUMWR_EXITWRSIABOVETHISQTY = 10000
  REAL, parameter  :: SUMET_EXITWRSIABOVETHISQTY = 10000
  ! SY: End for FLDAS
  character*3   :: fnest
  write(fnest,'(i3.3)') n

  alarmCheck = LIS_isAlarmRinging(LIS_rc,"Noah271 model alarm "//trim(fnest))
  if(alarmCheck) then 
     ! Get Julian day of year
     call LIS_tbotTimeUtil(julian_in,yr)
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!        dt = LIS_rc%ts       
        dt = noah271_struc(n)%ts

        tbot = noah271_struc(n)%noah(t)%tempbot
        if(LIS_rc%plevel.eq.3) then 
           write(LIS_logunit,*) 'tbot ',t,tbot
        endif

!=== STATIC VEGETATION PARAMETERS
! ----------------------------------------------------------------------
! VEGETATION PARAMETERS ARE DEPENDENT ON VEGETATION TYPE (INDEX)
!   SHDFAC: VEGETATION GREENNESS FRACTION
!   NROOT:  NUMBER OF ROOT LAYERS
!   RSMIN:  MIMIMUM STOMATAL RESISTANCE
!   RGL:    PARAMETER USED IN SOLAR RAD TERM OF
!           CANOPY RESISTANCE FUNCTION
!   HS:     PARAMETER USED IN VAPOR PRESSURE DEFICIT TERM OF
!           CANOPY RESISTANCE FUNCTION
!   SNUP:   THRESHOLD SNOW DEPTH (IN WATER EQUIVALENT M) THAT
!           IMPLIES 100% SNOW COVER
!   Z0:     Roughness Length
!   XLAI:   Leaf Area Index
! ----------------------------------------------------------------------
! SSIB VEGETATION TYPES (DORMAN AND SELLERS, 1989; JAM)
!  1:  BROADLEAF-EVERGREEN TREES  (TROPICAL FOREST)
!  2:  BROADLEAF-DECIDUOUS TREES
!  3:  BROADLEAF AND NEEDLELEAF TREES (MIXED FOREST)
!  4:  NEEDLELEAF-EVERGREEN TREES
!  5:  NEEDLELEAF-DECIDUOUS TREES (LARCH)
!  6:  BROADLEAF TREES WITH GROUNDCOVER (SAVANNA)
!  7:  GROUNDCOVER ONLY (PERENNIAL)
!  8:  BROADLEAF SHRUBS WITH PERENNIAL GROUNDCOVER
!  9:  BROADLEAF SHRUBS WITH BARE SOIL
! 10:  DWARF TREES AND SHRUBS WITH GROUNDCOVER (TUNDRA)
! 11:  BARE SOIL
! 12:  CULTIVATIONS (THE SAME PARAMETERS AS FOR TYPE 7)
! 13:  GLACIAL (THE SAME PARAMETERS AS FOR TYPE 11)
! ----------------------------------------------------------------------

        NROOT = NOAH271_STRUC(N)%NOAH(T)%nroot
        RSMIN = NOAH271_STRUC(N)%NOAH(T)%rsmin
        RGL   = NOAH271_STRUC(N)%NOAH(T)%rgl
        HS    = NOAH271_STRUC(N)%NOAH(T)%hs
        SNUP  = NOAH271_STRUC(N)%NOAH(T)%snup
        Z0    = NOAH271_STRUC(N)%NOAH(T)%z0
        XLAI  = NOAH271_STRUC(N)%NOAH(T)%lai
        if(LIS_rc%plevel.eq.3) then 
           write(LIS_logunit,*) 'veg parameters ',t,noah271_struc(n)%noah(t)%vegt
           write(LIS_logunit,*) 'nroot ',nroot
           write(LIS_logunit,*) 'rgl ',rgl
           write(LIS_logunit,*) 'hs ',hs
           write(LIS_logunit,*) 'snup ',snup
           write(LIS_logunit,*) 'z0 ',z0
           write(LIS_logunit,*) 'xlai ',xlai
        endif
        !=== MONTHLY VEGETATION PARAMETERS
        tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
        SHDFAC =  LIS_gfrac(n)%greenness(tid)
        noah271_struc(n)%noah(t)%vegip = shdfac
        !    Minimum greenness fraction
        SHDMIN=0.0

        !    If ground surface is bare:
        IF (noah271_struc(n)%noah(t)%vegt .EQ. LIS_rc%bareclass) SHDFAC = 0.0
        IF (noah271_struc(n)%noah(t)%vegt .EQ. LIS_rc%snowclass) SHDFAC = 0.0
        IF (noah271_struc(n)%noah(t)%vegt .EQ. LIS_rc%waterclass) SHDFAC = 0.0

        if(LIS_rc%plevel.eq.3) then 
           write(LIS_logunit,*) 'gfrac ',t,shdfac
        endif
!=== STATIC SOIL PARAMETERS
! ----------------------------------------------------------------------
!  SOIL PARAMETERS ARE DEPENDENT ON SOIL TYPE (INDEX)
!    SMCMAX: MAX SOIL MOISTURE CONTENT (POROSITY)
!    SMCREF: REFERENCE SOIL MOISTURE (ONSET OF SOIL MOISTURE
!             STRESS IN TRANSPIRATION)
!    SMCWLT: WILTING PT SOIL MOISTURE CONTENT
!    SMCDRY: AIR DRY SOIL MOIST CONTENT LIMITS
!    PSISAT: SATURATED SOIL POTENTIAL
!    DKSAT:  SATURATED SOIL HYDRAULIC CONDUCTIVITY
!    BEXP:   THE 'B' PARAMETER
!    DWSAT:  SATURATED SOIL DIFFUSIVITY
!    F1:     USED TO COMPUTE SOIL DIFFUSIVITY/CONDUCTIVITY
!    QUARTZ: SOIL QUARTZ CONTENT
! ----------------------------------------------------------------------
! SOIL TYPES   ZOBLER (1986)      COSBY ET AL (1984) (quartz cont.(1))
!  1        COARSE            LOAMY SAND         (0.82)
!  2        MEDIUM            SILTY CLAY LOAM    (0.10)
!  3        FINE              LIGHT CLAY         (0.25)
!  4        COARSE-MEDIUM     SANDY LOAM         (0.60)
!  5        COARSE-FINE       SANDY CLAY         (0.52)
!  6        MEDIUM-FINE       CLAY LOAM          (0.35)
!  7        COARSE-MED-FINE   SANDY CLAY LOAM    (0.60)
!  8        ORGANIC           LOAM               (0.40)
!  9        GLACIAL LAND ICE  LOAMY SAND         (NA using 0.82)
! ----------------------------------------------------------------------
       
!     SLDPTH(1) = 0.1         ! Soil layer thicknesses (m)
!     SLDPTH(1) = 0.05         ! Soil layer thicknesses (m)
!     SLDPTH(2) = 0.35
!     SLDPTH(1) = 0.1         ! Soil layer thicknesses (m)
!     SLDPTH(2) = 0.3
!     SLDPTH(3) = 0.6
!     SLDPTH(4) = 1.0

!-------------------------------------------------------------------------
! Soil layer thicknesses (m)
!-------------------------------------------------------------------------
        nsoil = noah271_struc(n)%nslay
        allocate(SLDPTH(nsoil))
        do i = 1,nsoil
           SLDPTH(i) = noah271_struc(n)%lyrthk(i)
        enddo

        SOILTYP = NOAH271_STRUC(N)%NOAH(T)%soiltype
        SMCMAX =  NOAH271_STRUC(N)%NOAH(T)%smcmax
        PSISAT =  NOAH271_STRUC(N)%NOAH(T)%psisat
        DKSAT  =  NOAH271_STRUC(N)%NOAH(T)%dksat
        BEXP   =  NOAH271_STRUC(N)%NOAH(T)%bexp
        QUARTZ =  NOAH271_STRUC(N)%NOAH(T)%quartz
        if(LIS_rc%plevel.eq.3) then 
           write(LIS_logunit,*) 'soil params '
           write(LIS_logunit,*) 'soiltype ',soiltyp
           write(LIS_logunit,*) 'smcmax ',smcmax
           write(LIS_logunit,*) 'psisat ',psisat
           write(LIS_logunit,*) 'dksat ',dksat
           write(LIS_logunit,*) 'bexp ',bexp
           write(LIS_logunit,*) 'quartz ',quartz
        endif
        !   The following 5 parameters are just given here for reference
        !    and to force static storage allocation.

        !       SMCREF =  NOAH%SOILP(6)
        !       SMCWLT =  NOAH%SOILP(7)
        !       SMCDRY =  NOAH%SOILP(8)
        !       DWSAT  =  NOAH%SOILP(9)
        !       F1     =  NOAH%SOILP(10)

        !-- Here is where the above five parameters are actually derived.
        !    SET TWO SOIL MOISTURE WILT, SOIL MOISTURE REFERENCE PARAMETERS
        SMLOW = 0.5
        !     Changed in 2.6 from 3 to 6 on June 2nd 2003
        !       SMHIGH = 3.0
        SMHIGH = 6.0

        if((smcmax.ne.0).and.(bexp.ne.0).and.(dksat.ne.0)) then 
           DWSAT  = BEXP * DKSAT * (PSISAT/SMCMAX)
           F1     = ALOG10(PSISAT) + BEXP*ALOG10(SMCMAX) + 2.0
           SMCREF1 = SMCMAX*(5.79E-9/DKSAT)**(1.0/(2.0*BEXP+3.0))
           SMCREF = SMCREF1 + (SMCMAX-SMCREF1) / SMHIGH
           SMCWLT1 = SMCMAX * (200.0/PSISAT)**(-1.0/BEXP)
           SMCWLT = SMCWLT1 - SMLOW * SMCWLT1
           !    Current version SMCDRY values equate to SMCWLT
           !        smcwlt = 0.02
           SMCDRY = SMCWLT
        else
           dwsat = 0.0
           f1 = 0.0
           smcref = 0.0
           smcwlt = 0.0
           smcdry = 0.0
        endif


! ----------------------------------------------------------------------
! KDT IS DEFINED BY REFERENCE REFKDT AND DKSAT; REFDK=2.E-6 IS THE SAT.
! DK. VALUE FOR THE SOIL TYPE 2
! ----------------------------------------------------------------------
        REFDK=2.0E-6
        REFKDT=3.0
        KDT = REFKDT * DKSAT/REFDK

! ----------------------------------------------------------------------
! FROZEN GROUND PARAMETER, FRZK, DEFINITION: ICE CONTENT THRESHOLD ABOVE
! WHICH FROZEN SOIL IS IMPERMEABLE REFERENCE VALUE OF THIS PARAMETER FOR
! THE LIGHT CLAY SOIL (TYPE=3) FRZK = 0.15 M.
! ----------------------------------------------------------------------
        FRZK=0.15
! ----------------------------------------------------------------------
! TO ADJUST FRZK PARAMETER TO ACTUAL SOIL TYPE: FRZK * FRZFACT
! ----------------------------------------------------------------------
        if(smcref.eq.0) then 
           frzx = 0.0
        else
           FRZFACT = (SMCMAX / SMCREF) * (0.412 / 0.468)
           FRZX = FRZK * FRZFACT
        endif
     
!=== SLOPE TYPE
! ----------------------------------------------------------------------
! CLASS PARAMETER 'SLOPETYP' WAS INCLUDED TO ESTIMATE LINEAR RESERVOIR
! COEFFICIENT 'SLOPE' TO THE BASEFLOW RUNOFF OUT OF THE BOTTOM LAYER.
! LOWEST CLASS (SLOPETYP=0) MEANS HIGHEST SLOPE PARAMETER = 1.
! DEFINITION OF SLOPETYP FROM 'ZOBLER' SLOPE TYPE:
! SLOPE CLASS  PERCENT SLOPE
! 1            0-8
! 2            8-30
! 3            > 30
! 4            0-30
! 5            0-8 & > 30
! 6            8-30 & > 30
! 7            0-8, 8-30, > 30
! 8            GLACIAL ICE
! 9            OCEAN/SEA
! ----------------------------------------------------------------------

!--  SLOPETYP = 3
        SLOPE = 1.0

!=== MONTHLY (QUARTERLY, for now) ALBEDO (SNOW-FREE)
        ALB =  LIS_alb(n)%ALBSF(tid)
        
        !    Maximum Albedo over very Deep Snow

        SNOALB = LIS_alb(n)%MXSNALB(tid)

        if(LIS_rc%plevel.eq.3) then 
           write(LIS_logunit,*) 'alb ',alb
           write(LIS_logunit,*) 'snoalb ',snoalb
        endif
!=== THE FOLLOWING BREAKS DOWN THE FORCING VARIABLES
        SFCTMP = noah271_struc(n)%noah(t)%tair/noah271_struc(n)%forc_count
        Q2     = noah271_struc(n)%noah(t)%qair/noah271_struc(n)%forc_count
        SOLDN  = noah271_struc(n)%noah(t)%swdown/noah271_struc(n)%forc_count
! NUWRF EMK...WRF now passes full LWDOWN, so this must be multiplied by
! emissivity.
!        LWDN   = noah271_struc(n)%noah(t)%lwdown/noah271_struc(n)%forc_count
        LWDN   = noah271_struc(n)%noah(t)%lwdown * &
             noah271_struc(n)%noah(t)%emiss / noah271_struc(n)%forc_count

        UWIND  = (noah271_struc(n)%noah(t)%uwind)*(noah271_struc(n)%noah(t)%uwind)
        VWIND  = (noah271_struc(n)%noah(t)%vwind)*(noah271_struc(n)%noah(t)%vwind)
        SFCSPD = SQRT( UWIND + VWIND )/noah271_struc(n)%forc_count
        SFCPRS = noah271_struc(n)%noah(t)%psurf/noah271_struc(n)%forc_count
        PRCP   = noah271_struc(n)%noah(t)%rainf/noah271_struc(n)%forc_count
   !Convective Precipitation (kg m-2 s-1ec)
        CPCP   = noah271_struc(n)%noah(t)%rainf_c/noah271_struc(n)%forc_count
        
        if((SFCTMP.eq.LIS_rc%udef).or.&
             (Q2.eq.LIS_rc%udef).or.&
             (SOLDN.eq.LIS_rc%udef).or.&
             (LWDN.eq.LIS_rc%udef).or.&
             (UWIND.eq.LIS_rc%udef).or.&
             (VWIND.eq.LIS_rc%udef).or.&
             (SFCPRS.eq.LIS_rc%udef).or.&
             (PRCP.eq.LIS_rc%udef)) then 
           print*, 'PROBLEM in FORCING ',sfctmp,q2,soldn,lwdn,uwind,vwind,&
                sfcspd,sfcprs,prcp,cpcp
           stop
        endif
        if(sfctmp.lt.200) then
           print*, 'problem in forcing ',sfctmp, q2,soldn,lwdn,uwind,vwind,&
                sfcspd,sfcprs,prcp,cpcp
           stop
        endif

        
        if(LIS_rc%plevel.eq.3) then 
           write(LIS_logunit,*) 'forcing ',sfctmp,q2,soldn,lwdn,uwind,vwind,&
                sfcspd,sfcprs,prcp,cpcp
        endif
        !-- Prevent Numerical Instability for Wind Speed
        if(SFCSPD.le.0.01) SFCSPD=0.01

        !-- Prevent Numerical Instability with HUMIDITY

        IF (Q2 .LT. 0.1E-5)  Q2 = 0.1E-5

! Calculate Saturation Specific Humidity (Kg/Kg) and
!    Saturation vapor pressure for water (Pa) based on Specific
!    Humidity(Kg/Kg), Temperature(K), and Pressure (Pa)
!
!      FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989: 'A
!      SHORT COURSE IN CLOUD PHYSICS', PERGAMON PRESS, 3rd ED.
!        Pablo J. Grunmann, 3/6/98.
!
!      QSAT  = Saturation Specific humidity (Kg/Kg)
!      ESAT  = Saturation vapor pressure for water (Pa)
!      EPS   = Water/(dry air) molecular mass ratio (epsilon)
!      E     = Saturation vapor pressure
!
!-- Function E(SFCTMP) = Sat. vapor pressure (in Pascal) at
!                   temperature T (uses Clausius-Clapeyron).
        ESAT = E(SFCTMP)

!-- CALCULATE SATURATION MIXING RATIO (PABLO GRUNMANN, 05/28/98)

#if (!defined COUPLED)
        Q2SAT = 0.622 * ESAT /(SFCPRS - (1.-0.622)*ESAT)
        IF (Q2 .GE.  Q2SAT)  Q2 = Q2SAT*0.99

        ! Set the heights of the forcing variables (D. Mocko; J. Santanello)
        if ((noah271_struc(n)%zh.lt.0.0).or.(noah271_struc(n)%zm.lt.0.0)) then
           noah271_struc(n)%zh = noah271_struc(n)%noah(t)%fheight
           noah271_struc(n)%zm = noah271_struc(n)%noah(t)%fheight
           if (noah271_struc(n)%noah(t)%fheight.le.0.0) then
              write(LIS_logunit,*) 'Forcing height less than or'
              write(LIS_logunit,*) 'equal to zero!  Stopping run'
              call LIS_endrun()
           endif
        else
           if (noah271_struc(n)%zh.lt.noah271_struc(n)%zm) then
              ZFAC = LOG((noah271_struc(n)%zh + noah271_struc(n)%noah(t)%z0)&
                   /noah271_struc(n)%noah(t)%z0) /                    &
                   LOG((noah271_struc(n)%zm + noah271_struc(n)%noah(t)%z0)&
                   /noah271_struc(n)%noah(t)%z0)
              SFCSPD = ZFAC*SFCSPD
              !-- Prevent Numerical Instability for Wind Speed
              if (SFCSPD.le.0.01) SFCSPD=0.01
              UWIND = (ZFAC*SQRT(UWIND))*(ZFAC*SQRT(UWIND))
              VWIND = (ZFAC*SQRT(VWIND))*(ZFAC*SQRT(VWIND))
           endif
        endif
        ZLVL = noah271_struc(n)%zh
#else
        Q2SAT = noah271_struc(n)%noah(t)%q2sat
        ZLVL = 0.5*noah271_struc(n)%noah(t)%z
        !     print*, 'noah271_main ',q2sat, z
#endif

        !-- CALCULATE SLOPE OF SAT. SPECIFIC HUMIDITY CURVE FOR PENMAN: DQSDT2

        DQSDT2 = DQSDT (SFCTMP, SFCPRS)

        !-- CALC VIRTUAL TEMPS AND POTENTIAL TEMPS AT GRND (SUB 1) AND AT
        !    THE 1ST MDL LVL ABV THE GRND (SUB 2). EXPON IS CP DIVD BY R.

        TH2 = SFCTMP + ( 0.0098 * ZLVL )
        T2V = SFCTMP * (1.0 + 0.61 * Q2 )

        T1V  = noah271_struc(n)%noah(t)%T1 * (1.0 + 0.61 * Q2 )
        TH2V = TH2 * (1.0 + 0.61 * Q2 )

!==  OPTIONAL SUBROUTINE:  Calculate LW Radiation (Down)  =================
!      CALL OBTLWDN(SFCTMP,LWDN)

!--  Phota Thermal Unit (PTU)
        
        PTU    = 0.10
     
!--  Initialize SOILM for 1st timestep water balance
!       SOILM = 0.0
!--  Initialize ROOT ZONE COLUMN SOIL MOISTURE IN METERS (SOILRZ)
        SOILRZ = 0.0
        SOILRZMAX = 0.0 ! SY
!--  Initialize TOP 1-METER COLUMN SOIL MOISTURE IN METERS (SOIL1M)
        SOIL1M = 0.0

!==  CALCULATE CH (EXCHANGE COEFFICIENT)  =================================
!
!   CH IS THE SFC EXCHANGE COEFFICIENT FOR HEAT/MOISTURE
!   CM IS THE SFC EXCHANGE COEFFICIENT FOR MOMENTUM
!
! IMPORTANT NOTE: TO CALCULATE THE SFC EXCHANGE COEF (CH) FOR HEAT AND
!                 MOISTURE, SUBROUTINE SFCDIF BELOW CAN:
!
!    A)  BE CALLED HERE FROM THE DRIVER, THUS CH IS INPUT TO SFLX
!           (AS IS TYPICAL IN A COUPLED ATMOSPHERE/LAND MODEL), OR
!  * B)  BE CALLED INTERNALLY IN ROUTINE SFLX (THUS CH IS OUTPUT FROM SFLX),
!          BEFORE THE CALL TO ROUTINE "PENMAN"
!
!   OPTION B IS THE DEFAULT HERE.  THAT IS, IN THE UNCOUPLED, OFF-LINE LSM
!   REPRESENTED HEREIN BY THIS DRIVER, WE CALL SFCDIF LATER IN ROUTINE SFLX.
!
!   THE ROUTINE SFCDIF REPRESENTS THE SO-CALLED "SURFACE LAYER" OR THE
!   "CONSTANT FLUX LAYER" (THE LOWEST 20-100 M OF THE ATMOSPHERE).
!   HENCE ROUTINE SFCDIF EMBODIES THE "ATMOSPHERIC AERODYNAMIC RESISTANCE".
!
!   TO ENABLE THE FLEXIBILITY OF EITHER OPTION A OR B, WE PASS
!   THE ARGUMENTS "CH", "CM", AND "SFCSPD" (WIND SPEED:JUST CALCULATED ABOVE)
!   TO ROUTINE SFLX TO SUPPORT OPTION B -- THAT IS, FOR INPUT TO THE CALL TO
!   ROUTINE SFCDIF THEREIN.  IN OPTION A, THE ARGUMENTS "SFCSPD" AND "CM"
!   ARE NEITHER NEEDED IN ROUTINE SFLX, NOR ALTERED BY ROUTINE SFLX.
!
!     IF ONE CHOOSES OPTION A, THEN ONE MUST
!      1 - ACTIVATE (UNCOMMENT) THE CALL TO SFCDIF BELOW,
!      2 - ACTIVATE (UNCOMMENT) THE ASSIGNMENT OF "Z0" AND "CZIL" NEXT BELOW
!      3 - DE-ACTIVATE (COMMENT OUT) THE CALL TO SFCDIF IN ROUTINE SFLX.
!
!   Z0 and CZIL:
!
!   THE ROUGHNESS LENGTH PARAMETERS "Z0" AND "CZIL" MUST BE SET HERE IN THE
!   DRIVER TO SUPPORT THE "OPTION-A", I.E. THE CALL TO SFCDIF BELOW.  IN SO
!   DOING, THE "Z0" AND "CZIL" ASSIGNED HERE MUST CORRESPOND TO THEIR
!   VALUES, CALLED BY SFLX JUST BEFORE CALL SFCDIF.
!   THUS THE VALUE OF "Z0" ASSIGNED HERE MUST CORRESPOND TO THAT ASSIGNED
!   FOR THE CHOSEN VEG CLASS THAT WAS ALREADY INPUT EARLIER IN THE DRIVER.
!
!   BECAUSE OF THE IMPLICIT ITERATIVE NATURE OF THE "PAULSON" SURFACE-LAYER
!   SCHEME USED IN ROUTINE SFCDIF, CH AND CM ARE CO-DEPENDENT.  SIMILARLY,
!   THE IMPLICIT NATURE OF THE SFCDIF SCHEME ALSO REQUIRES THAT FOR EITHER
!   OPTION A OR B, CH AND CM MUST BE INITIALIZED EARLIER IN THE DRIVER BEFORE
!   THE START OF THE TIME-STEP LOOP, AS WELL AS BE CARRIED FORWARD FROM
!   TIME STEP TO TIME STEP AS "STATE VARIABLES", BECAUSE THE VALUES OF
!   CH AND CM FROM A PREVIOUS TIME STEP REPRESENT THE FIRST-GUESS VALUES FOR
!   THE CALL TO SFCDIF IN THE PRESENT TIME STEP.
!
!   SOME USERS MAY CHOOSE TO EXECUTE AN ENTIRELY DIFFERENT SCHEME IN PLACE OF
!   ROUTINE SFCDIF HERE, E.G. AN EXPLICIT SCHEME SUCH AS LOUIS (1979) THAT
!   EMPLOYS NO ITERATION AND HAS NO REQUIREMENT TO CARRY CH AND CM FORWARD
!   AS STATE VARIABLES FROM TIME STEP TO TIME STEP.  IN THAT CASE, IN
!   OPTION A, THE ROUTINE SHOULD BE CALLED HERE IN THE DRIVER AFTER ALL
!   NECESSARY INPUT ARGUMENTS FOR IT ARE DEFINED AT THIS POINT, OR CALLED IN
!   ROUTINE SFLX, AT THE POINT SFCDIF IS CALLED.
!
!      CALL SFCDIF ( Z, Z0, T1V, TH2V, SFCSPD,CZIL, CM, CH )
!
! ---------------------------------------------------------------------|
! INITIALIZE CH, CM (NOTE: initial these before time loop)
!      CH=1.E-4
!      CM=1.E-4
! 1998 May 22 0030 (Julian= 142) typical values initialization
!      CH=  0.0150022404
!      CM=  0.0205970779
! **NOTE: TRYING THESE VALUES AS TEST! 
! ---------------------------------------------------------------------|
#if (!defined COUPLED)     

! If the offline forcing is not providing the exchange coefficients,
! then let Noah2.7.1 calculate its own surface aerodynamic conductance.
        if ( noah271_struc(n)%forcing_ch == 0 ) then
! CH and CM are initialized in either coldstart or restart.
! Use current values to speed up the convergence of the iterative
! method in sfcdif.
!      noah271_struc(n)%noah(t)%CH=  1.0E-4
!      noah271_struc(n)%noah(t)%CM=  1.0E-4
           CALL SFCDIF(ZLVL,Z0,T1V,TH2V,                                    &
                SFCSPD,noah271_struc(n)%noah(t)%czil,                 &
                noah271_struc(n)%noah(t)%CM,noah271_struc(n)%noah(t)%CH,&
                USTAR)
        else
           ! Still call sfcdif because we need ustar, but do not overwrite
           ! the provided ch and cm values.
           cm = noah271_struc(n)%noah(t)%cm
           ch = noah271_struc(n)%noah(t)%ch/noah271_struc(n)%forc_count
           CALL SFCDIF(ZLVL,noah271_struc(n)%noah(t)%z0,T1V,TH2V,           &
                SFCSPD,noah271_struc(n)%noah(t)%czil,cm,ch,USTAR)
        endif
#else
        ! Still call sfcdif because we need ustar, but do not overwrite
        ! the provided ch and cm values.
        cm = noah271_struc(n)%noah(t)%cm
        ch = noah271_struc(n)%noah(t)%ch
        CALL SFCDIF(ZLVL,Z0,T1V,TH2V,                                       &
             SFCSPD,noah271_struc(n)%noah(t)%czil,cm,ch,USTAR)
#endif

!=== MAIN CALL TO LAND-SURFACE PHYSICS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! ---------------------------------------------------------------------|
! NEW FOR NOAH 2.7.1
! ---------------------------------------------------------------------|
!--  Determine Precipitation Type Before Call to SFLX (noah271_physics.F90)
        FFROZP = 0.0
        IF (PRCP .GT. 0.0) THEN
           IF (SFCTMP .LE. T0) THEN
              FFROZP = 1.0
           ENDIF
        ENDIF

        ALBEDO = ALB
        SNCOVR = 0.
        SOLNET = 0.
        ! ---------------------------------------------------------------------|
        ! QC FOR SNOW PROPERTIES. JESSE 20050106
        ! ---------------------------------------------------------------------|
        IF ( noah271_struc(n)%noah(t)%SNEQV .EQ. 0.0 ) THEN
           noah271_struc(n)%noah(t)%SNOWH = 0.0
        ELSE
           IF ( noah271_struc(n)%noah(t)%SNOWH .EQ. 0.0 ) THEN
              noah271_struc(n)%noah(t)%SNOWH  = noah271_struc(n)%noah(t)%SNEQV * 10.
           ENDIF
        ENDIF

        ! ---------------------------------------------------------------------|
        ! RESET FLUXES BEFORE SFLX
        ! ---------------------------------------------------------------------|
        if(noah271_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass) then
           ice = 1
        else
           ice = 0
        endif

#if (defined COUPLED)
        ice = 0 
        if(noah271_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass) then 
           ! open water case
           if(noah271_struc(n)%noah(t)%xice.eq.1) then 
              ice = 1
              noah271_struc(n)%noah(t)%smc = 1.0
              noah271_struc(n)%noah(t)%stc = 273.16
           endif
        else
           ! land or sea-ice case
           if(noah271_struc(n)%noah(t)%xice.eq.1) then
              ice = 1
              noah271_struc(n)%noah(t)%smc = 1.0
           endif
           DQSDT2 = Q2SAT*A23M4/(SFCTMP-A4)**2

           if(noah271_struc(n)%noah(t)%sneqv.gt.0.0) then
              ! snow on surface (use ice saturation properties)
              call SNFRAC(noah271_struc(n)%noah(t)%sneqv,                 &
                   SNUP,SALP,noah271_struc(n)%noah(t)%snowh,SNCOVR)
              noah271_struc(n)%noah(t)%sca = sncovr ! EMK NUWRF
              SFCTSNO = SFCTMP
              E2SAT = 611.2*EXP(6174.*(1./273.15 - 1./SFCTSNO))
              Q2SATI = 0.622*E2SAT/(SFCPRS-E2SAT)
              Q2SATI = Q2SATI/(1.0+Q2SATI)    ! spec. hum.
              if(noah271_struc(n)%noah(t)%T1.gt.273.15) then
                 ! warm ground temps, weight the saturation between ice and water according to SNOWC
                 Q2SAT = Q2SAT*(1.-SNCOVR) + Q2SATI*SNCOVR
                 DQSDT2 = DQSDT2*(1.-SNCOVR) + Q2SATI*6174./(SFCTSNO**2)*SNCOVR
              else
                 ! cold ground temps, use ice saturation only
                 Q2SAT = Q2SATI
                 DQSDT2 = Q2SATI*6174./(SFCTSNO**2)
              endif
              ! for snow cover fraction at 0 C, ground temp will not change, so DQSDT2 effectively zero
              if ((noah271_struc(n)%noah(t)%T1.gt.273.).and.              &
                   (SNCOVR.gt.0.)) DQSDT2 = DQSDT2*(1.-SNCOVR)
           endif

           if(ICE.eq.0) then
              !           TBOT= noah271_struc(n)%noah(t)%tempbot
           else
              TBOT=271.16
           endif
        endif
#endif

        evp = 0
        eta = 0 
        shtflx = 0 
        ec = 0 
        edir = 0
        et = 0 
        ett = 0 
        esnow = 0 
        drip = 0
        dew = 0
        beta = 0 
        etp = 0 
        gflx = 0 
        flx1 = 0 
        flx2 = 0 
        flx3 = 0
        snomlt = 0 
        sncovr = 0 
        noah271_struc(n)%noah(t)%runoff1 = 0
        noah271_struc(n)%noah(t)%runoff2 = 0
        runoff3 = 0 
        rc = 0 
        pc = 0 
        rcs = 0 
        rct = 0 
        rcq = 0 
        rcsoil = 0 
        mstavrz = 0 
        mstavtot = 0 
        soilm = 0 
        soilmtc = 0
        tsoil = 0
        q1 = 0
        do i=1,noah271_struc(n)%nslay
           soilmtc = soilmtc +                                            &
                noah271_struc(n)%noah(t)%smc(i)*LIS_CONST_RHOFW*          &
                noah271_struc(n)%lyrthk(i)
        enddo
        ! First attempt at calculating DelSurfHeat - D. Mocko
        !     soilhtc = (((NOAH271_STRUC(N)%NOAH(T)%SH2O(1)*CH2O)               &
        !               + ((1.0-SMCMAX)*CSOIL)                                  &
        !               + ((SMCMAX-NOAH271_STRUC(N)%NOAH(T)%SMC(1))*CAIR)       &
        !               + ((NOAH271_STRUC(N)%NOAH(T)%SMC(1)-NOAH271_STRUC(N)%NOAH(T)%SH2O(1))*CICE)) &
        !               * noah271_struc(n)%noah(t)%stc(1) * SLDPTH(1))
        soilhtc = 0.0
        startht = soilhtc
        startsm = soilmtc
        startswe = noah271_struc(n)%noah(t)%sneqv*LIS_CONST_RHOFW
        startint = noah271_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW
        IF (noah271_struc(n)%noah(t)%vegt .EQ. LIS_rc%urbanclass) then
           SHDFAC=0.05
           RSMIN=400.0
           SMCMAX = 0.45
           SMCREF = 0.42
           SMCWLT = 0.40
        endif

        LIS_gfrac(n)%greenness(t) = shdfac
        if ( shdfac .lt. 0 ) then
           write(LIS_logunit,*) 'SHDFAC < 0 ',t, noah271_struc(n)%noah(t)%vegt, &
                noah271_struc(n)%noah(t)%soiltype,shdfac
        endif
        !     print*, t, noah271_struc(n)%noah(t)%smcmax, & 
        !          noah271_struc(n)%noah(t)%psisat, &
        !          noah271_struc(n)%noah(t)%dksat, &
        !          noah271_struc(n)%noah(t)%bexp, &
        !          noah271_struc(n)%noah(t)%quartz, &
        !          smcwlt, &
        !          noah271_struc(n)%noah(t)%smc(1)
        !     if(t.eq.157293) then 
        !        write(231,*) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, noah271_struc(n)%noah(t)%snowh, &
        !             noah271_struc(n)%noah(t)%sneqv!, LIS_domain(n)%grid(LIS_domain(n)%gindex(&
        !             LIS_domain(n)%tile(t)%col, LIS_domain(n)%tile(t)%row))%lat, &
        !             LIS_domain(n)%grid(LIS_domain(n)%gindex(&
        !             LIS_domain(n)%tile(t)%col, LIS_domain(n)%tile(t)%row))%lon
        !     endif

        EMISSI_OUT = noah271_struc(n)%noah(t)%emiss
        CALL SFLX (T,                                                 &
             noah271_struc(n)%noah(t)%sbeta, &
             noah271_struc(n)%noah(t)%fxexp, &
             noah271_struc(n)%noah(t)%czil, &
             noah271_struc(n)%noah(t)%CFACTR,noah271_struc(n)%noah(t)%CMCMAX,&
             noah271_struc(n)%noah(t)%RSMAX,noah271_struc(n)%noah(t)%TOPT,&
             noah271_struc(n)%noah(t)%REFDK,noah271_struc(n)%noah(t)%REFKDT,&
             noah271_struc(n)%noah(t)%CSOIL,noah271_struc(n)%noah(t)%FRZK,&
             FFROZP,ICE,DT,Z,NSOIL,SLDPTH,                                &
             LWDN,SOLDN,SOLNET,SFCPRS,PRCP,SFCTMP,Q2,SFCSPD,              &
             TH2,Q2SAT,DQSDT2,                                            &
             SLOPE,SHDFAC,SHDMIN,PTU,ALB,SNOALB,                          &
             RSMIN,RGL,HS,SNUP,Z0,noah271_struc(n)%noah(t)%emiss,            &
             XLAI,NROOT,                                                  &
             PSISAT,BEXP,DKSAT,SMCMAX,QUARTZ,DWSAT,                       &
             SMCWLT,SMCREF,SMCDRY,F1,KDT,FRZX,FRZFACT,TBOT,               &
             noah271_struc(n)%noah(t)%CMC,noah271_struc(n)%noah(t)%T1,          &
             noah271_struc(n)%noah(t)%STC,noah271_struc(n)%noah(t)%SMC,         &
             noah271_struc(n)%noah(t)%SH2O,noah271_struc(n)%noah(t)%SNOWH,      &
             noah271_struc(n)%noah(t)%SNEQV,ALBEDO,noah271_struc(n)%noah(t)%CH, &
             noah271_struc(n)%noah(t)%CM,                                    &
             EVP,ETA,SHTFLX,ETA_KINEMATIC,                                &
             EC,EDIR,ET,ETT,ESNOW,DRIP,DEW,                               &
             BETA,ETP,GFLX,                                               &
             FLX1,FLX2,FLX3,                                              &
             SNOMLT,SNCOVR,                                               &
             NOAH271_STRUC(N)%NOAH(T)%RUNOFF1,NOAH271_STRUC(N)%NOAH(T)%RUNOFF2,RUNOFF3,                                     &
             RC,PC,RCS,RCT,RCQ,RCSOIL,Q1,                                 &
             MSTAVRZ,MSTAVTOT,SOILM,ETANRG,TSOIL,EMISSI_OUT)

       noah271_struc(n)%noah(t)%sca = sncovr ! EMK NUWRF

#if 0
        CALL NOAH271_SFLX ( &
             ICE,DT,Z,NSOIL,SLDPTH,&
             LWDN,SOLDN,SFCPRS,PRCP,SFCTMP,Q2,SFCSPD, &
             TH2,Q2SAT,DQSDT2,&
             SLOPE,SHDFAC,ALB,SNOALB, &
             RSMIN,RGL,HS,SNUP,Z0,XLAI,NROOT,&
             PSISAT,BEXP,DKSAT,SMCMAX,QUARTZ,DWSAT, &
             SMCWLT,SMCREF,SMCDRY,KDT,FRZX,TBOT, &
             NOAH271_STRUC(N)%NOAH(T)%CMC,NOAH271_STRUC(N)%NOAH(T)%T1,&
             NOAH271_STRUC(N)%NOAH(T)%STC,NOAH271_STRUC(N)%NOAH(T)%SMC,&
             NOAH271_STRUC(N)%NOAH(T)%SH2O, &
             NOAH271_STRUC(N)%NOAH(T)%SNOWH,NOAH271_STRUC(N)%NOAH(T)%SNEQV,&
             ALBEDO,NOAH271_STRUC(N)%NOAH(T)%CH,NOAH271_STRUC(N)%NOAH(T)%CM,&
             EVP,ETA,SHTFLX,&
             EC,EDIR,ET,ETT,ESNOW,DRIP,DEW, &
             BETA,ETP,GFLX, &
             FLX1,FLX2,FLX3,&
             SNOMLT,SNCOVR,&
             NOAH271_STRUC(N)%NOAH(T)%RUNOFF1,NOAH271_STRUC(N)%NOAH(T)%RUNOFF2,RUNOFF3, &
             RC,PC,RCS,RCT,RCQ,RCSOIL,&
             MSTAVRZ,MSTAVTOT,SOILM,ETANRG)
#endif
        rho = sfcprs/(287.04*t2v)
        noah271_struc(n)%noah(t)%tau =  -rho*ustar*ustar/sfcspd     
        noah271_struc(n)%noah(t)%tauu = noah271_struc(n)%noah(t)%tau*uwind
        noah271_struc(n)%noah(t)%tauv = noah271_struc(n)%noah(t)%tau*vwind

        !******************************************************************
        !   CALCULATE UPWARD LONGWAVE RAD USING UPDATED SKIN TEMPERATURE
        !   MIKE EK'S SNOW EMISSIVITY. JESSE 20041228
        T14 = noah271_struc(n)%noah(t)%T1 * noah271_struc(n)%noah(t)%T1 * noah271_struc(n)%noah(t)%T1 * noah271_struc(n)%noah(t)%T1 

        FUP = 5.67E-8 * T14 * (0.95*SNCOVR + 1.0*(1.0-SNCOVR))

        !-  CALCULATE SNOW DENSITY FROM OUTPUT:
        IF (noah271_struc(n)%NOAH(T)%SNEQV .GT. 0) THEN
           ! Commented out Jesse's snow treatment, as it
           ! violates mass balance calculations - D. Mocko
           !JESSE 20050901 VERY LOW SNEQV TREATMNENT
           !        IF(noah271_struc(n)%NOAH(T)%SNEQV.LE.0.001 .AND. noah271_struc(n)%NOAH(T)%T1.GT.273.15) THEN
           !           noah271_struc(n)%NOAH(T)%SNEQV = 0.
           !           noah271_struc(n)%NOAH(T)%SNOWH = 0.
           !        ENDIF
           !SNDENS=noah271_struc(n)%NOAH(T)%SNEQV/noah271_struc(n)%NOAH(T)%SNOWH
           continue
        ELSE
           noah271_struc(n)%NOAH(T)%SNOWH=0.0
           !SNDENS=0.0
        ENDIF

        !   CALCULATE RESIDUAL OF ALL SURFACE ENERGY BALANCE EQN TERMS.

        !      GFLX = -GFLX
        !      F = SOLDN*(1.0-ALBEDO) + LWDN
        !      RES    = F - SHTFLX - GFLX - ETA - FUP - FLX1 - FLX2 - FLX3
        !         ENDIF

        !      WRITE(LIS_logunit,*)'  --------------------------------------'
        !      WRITE(LIS_logunit,*)'  State Variables '
        !      WRITE(LIS_logunit,*)'  --------------------------------------'
        !      WRITE(*,*) NOAH271_STRUC(N)%NOAH(T)%T1,' T1...Skin temperature (K)'
        !      WRITE(*,*)(NOAH271_STRUC(N)%NOAH(T)%STC(IJ), IJ=1,NSOIL),' STC'
        !      WRITE(*,*)(NOAH271_STRUC(N)%NOAH(T)%SMC(IJ), IJ=1,NSOIL),' SMC'
        !      WRITE(*,*)(NOAH271_STRUC(N)%NOAH(T)%SH2O(IJ), IJ=1,NSOIL),' SH2O'
        !      WRITE(*,*) NOAH271_STRUC(N)%NOAH(T)%CMC,' CMC...Canopy water content (m)'
        !      WRITE(*,*) NOAH271_STRUC(N)%NOAH(T)%SNOWH,' SNOWH...Actual snow depth (m)'
        !      WRITE(*,*) NOAH271_STRUC(N)%NOAH(T)%SNEQV,' SNEQV...Water equiv snow depth (m)'
        !      WRITE(*,*) 'CH= ',NOAH271_STRUC(N)%NOAH(T)%CH,'   CM= ',NOAH271_STRUC(N)%NOAH(T)%CM
        !      WRITE(LIS_logunit,*)'  --------------------------------------'

        !=== Collect the output variables into NOAH271_STRUC(N)%NOAH(T)%RETURN
        noah271_struc(n)%noah(t)%swnet = soldn*(1.0-albedo)
        noah271_struc(n)%noah(t)%lwnet = EMISSI_OUT*(5.67E-8)*(NOAH271_STRUC(N)%NOAH(T)%T1**4.0)-LWDN
        noah271_struc(n)%noah(t)%qle = eta
        !     print*, 'qle ',eta, shtflx, gflx, noah271_struc(n)%noah(t)%swnet-&
        !          noah271_struc(n)%noah(t)%lwnet-eta-shtflx+gflx
        !     stop
        noah271_struc(n)%noah(t)%beta = beta

        !     noah271_struc(n)%noah(t)%etanrg = etanrg
        noah271_struc(n)%noah(t)%qh = shtflx
        noah271_struc(n)%noah(t)%qg = gflx
        noah271_struc(n)%noah(t)%qf = flx3
        noah271_struc(n)%noah(t)%qv = esnow
        noah271_struc(n)%noah(t)%qtau = SFCPRS / (R * T2V) * noah271_struc(n)%noah(t)%CM * USTAR
        noah271_struc(n)%noah(t)%qa = -flx1
        ! First attempt at calculating DelSurfHeat - D. Mocko
        !     soilhtc = (((NOAH271_STRUC(N)%NOAH(T)%SH2O(1)*CH2O)               &
        !               + ((1.0-SMCMAX)*CSOIL)                                  &
        !               + ((SMCMAX-NOAH271_STRUC(N)%NOAH(T)%SMC(1))*CAIR)       &
        !               + ((NOAH271_STRUC(N)%NOAH(T)%SMC(1)-NOAH271_STRUC(N)%NOAH(T)%SH2O(1))*CICE)) &
        !               * noah271_struc(n)%noah(t)%stc(1) * SLDPTH(1))
        soilhtc = 0.0
        noah271_struc(n)%noah(t)%delsurfheat = startht - soilhtc
        noah271_struc(n)%noah(t)%delcoldcont = (FLX2 - esnow) * DT
        noah271_struc(n)%noah(t)%eta_kinematic = eta_kinematic

        if (sfctmp .le. LIS_CONST_TKFRZ) then
           noah271_struc(n)%noah(t)%snowfall = prcp
           noah271_struc(n)%noah(t)%pcp = 0.0
        else
           noah271_struc(n)%noah(t)%snowfall = 0.0
           noah271_struc(n)%noah(t)%pcp = prcp
        endif
        noah271_struc(n)%noah(t)%evap = evp
        noah271_struc(n)%noah(t)%qs = noah271_struc(n)%noah(t)%runoff1*LIS_CONST_RHOFW
        noah271_struc(n)%noah(t)%qsb = noah271_struc(n)%noah(t)%runoff2*LIS_CONST_RHOFW
        noah271_struc(n)%noah(t)%qsm = snomlt*LIS_CONST_RHOFW/dt
        do i=1,noah271_struc(n)%nslay
           noah271_struc(n)%noah(t)%soilmoist(i) = &
                noah271_struc(n)%noah(t)%smc(i)*LIS_CONST_RHOFW*            &
                noah271_struc(n)%lyrthk(i)
        enddo
        soilmtc = 0
        do i=1,noah271_struc(n)%nslay
           soilmtc = soilmtc + noah271_struc(n)%noah(t)%soilmoist(i)
        enddo
        noah271_struc(n)%noah(t)%delsoilmoist = soilmtc - startsm
        noah271_struc(n)%noah(t)%delswe = (noah271_struc(n)%noah(t)%sneqv*LIS_CONST_RHOFW) - startswe
        noah271_struc(n)%noah(t)%delintercept = (noah271_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW) - startint
        NOAH271_STRUC(N)%NOAH(T)%albedo=ALBEDO
        !     if(noah271_struc(n)%noah(t)%sneqv.gt.1000) then 
        noah271_struc(n)%noah(t)%swe = noah271_struc(n)%noah(t)%sneqv*LIS_CONST_RHOFW
        noah271_struc(n)%noah(t)%sca = SNCOVR

        ! NOTE:  Soil temperature for each layer is passed directly to output routines.
        NOAH271_STRUC(N)%NOAH(T)%soilwet = MSTAVTOT
        NOAH271_STRUC(N)%NOAH(T)%ecanop =  ec
        NOAH271_STRUC(N)%NOAH(T)%tveg = ETT
        NOAH271_STRUC(N)%NOAH(T)%esoil =EDIR
        NOAH271_STRUC(N)%NOAH(T)%subsnow = esnow
        noah271_struc(n)%noah(t)%canopint = noah271_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW

        ! ROOT ZONE COLUMN SOIL MOISTURE IN METERS (SOILRZ)
        ! define root zone soil moisture as a 1km average. 
        !     do k = 1,nroot
        !        soilrz = soilrz+(noah271_struc(n)%noah(t)%smc(k)*sldpth(k)*LIS_CONST_RHOFW)
        !     end do     
        !     noah271_struc(n)%noah(t)%rootmoist = soilrz
        soilrz = 0.0
        do k=1,nroot
           soilrz = soilrz + (noah271_struc(n)%noah(t)%smc(k)*noah271_struc(n)%lyrthk(k)*LIS_CONST_RHOFW)
           soilrzmax = soilrzmax + (noah271_struc(n)%noah(t)%smcmax*sldpth(k)*   &
                              LIS_CONST_RHOFW) ! SY
        enddo
        noah271_struc(n)%noah(t)%rootmoist =  soilrz   !in volumetric
        roottemp = 0 
        do k=1,nroot
           roottemp = roottemp + noah271_struc(n)%noah(t)%stc(k)*noah271_struc(n)%lyrthk(k)
        enddo

        noah271_struc(n)%noah(t)%q1        = q1
        noah271_struc(n)%noah(t)%soilm     = soilm*1000.0

        noah271_struc(n)%noah(t)%chs2 = noah271_struc(n)%noah(t)%cqs2
        if(q1.gt.noah271_struc(n)%noah(t)%qsfc) then 
           noah271_struc(n)%noah(t)%cqs2 = noah271_struc(n)%noah(t)%chs2
        endif
        noah271_struc(n)%noah(t)%qsfc = q1/(1-q1)

        ! ADDITIONAL AFWA OUTPUTS
        if (SFCTMP .lt. noah271_struc(n)%noah(t)%tair_agl_min) then
           noah271_struc(n)%noah(t)%tair_agl_min = SFCTMP
           ! Also save RH corresponding to the minimum surface temperature.
           !         print*, 'noah main',t, q2, q2sat
           noah271_struc(n)%noah(t)%rhmin = Q2/Q2SAT
        endif
        if(noah271_struc(n)%noah(t)%albedo.lt.0) then 
           print*, 'neg. albedo ',t, noah271_struc(n)%noah(t)%albedo        
           stop
        endif

        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_SWNET,value=soldn *       &
             (1.0-noah271_struc(n)%noah(t)%albedo),vlevel=1,unit="W m-2",&
             direction="DN",valid_min=0.0, valid_max=1200.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_LWNET,vlevel=1,           &
             value=(-1.0*((5.67E-8)*(NOAH271_STRUC(N)%NOAH(T)%T1**4.0)*   &
             EMISSI_OUT - LWDN)),unit="W m-2",direction="DN",&
             valid_min=-500.0, valid_max=510.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QLE,value=eta,            &
             vlevel=1,unit="W m-2",direction="UP",&
             valid_min=-700.0, valid_max=700.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QH,value=shtflx,          &
             vlevel=1,unit="W m-2",direction="UP",&
             valid_min=-600.0,valid_max=600.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QG,value=-gflx,           &
             vlevel=1,unit="W m-2",direction="DN",&
             valid_min=-500.0,valid_max=500.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QF,value=flx3,            &
             vlevel=1,unit="W m-2",direction="S2L",&
             valid_min=-1200.0,valid_max=1200.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QV,value=esnow,           &
             vlevel=1,unit="W m-2",direction="S2V",&
             valid_min=-600.0,valid_max=600.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QTAU,                     &
             value=(SFCPRS/(R*T2V)*noah271_struc(n)%noah(t)%CM*USTAR),    &
             vlevel=1,unit="N m-2",direction="DN",&
             valid_min=-100.0,valid_max=100.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QA,value=-flx1,           &
             vlevel=1,unit="W m-2",direction="DN",&
             valid_min=-50.0,valid_max=50.0,&
             surface_type = LIS_rc%lsm_index)
        ! First attempt at calculating DelSurfHeat - D. Mocko
        !     soilhtc = (((NOAH271_STRUC(N)%NOAH(T)%SH2O(1)*CH2O)               &
        !               + ((1.0-SMCMAX)*CSOIL)                                  &
        !               + ((SMCMAX-NOAH271_STRUC(N)%NOAH(T)%SMC(1))*CAIR)       &
        !               + ((NOAH271_STRUC(N)%NOAH(T)%SMC(1)-NOAH271_STRUC(N)%NOAH(T)%SH2O(1))*CICE)) &
        !               * noah271_struc(n)%noah(t)%stc(1) * SLDPTH(1))
        soilhtc = 0.0
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_DELSURFHEAT,              &
             value=(startht-soilhtc),vlevel=1,unit="J m-2",direction="INC",&
             valid_min=-500.0,valid_max=500.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_DELCOLDCONT,              &
             value=((flx2 - esnow) * DT),vlevel=1,unit="J m-2",direction="INC",&
             valid_min=-600.0,valid_max=1000.0,&
             surface_type = LIS_rc%lsm_index)
        !Bowen Ratio - sensible/latent
        if (eta.gt.0) then 
           call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_BR,value=(shtflx/eta), &
                vlevel=1,unit="-",direction="-",&
                valid_min=0.0,valid_max=10.0,&
                surface_type = LIS_rc%lsm_index)
        else
           call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_BR,value=0.0,          &
                vlevel=1,unit="-",direction="-",&
                valid_min=0.0,valid_max=10.0,&
                surface_type = LIS_rc%lsm_index)
        endif
        !Evaporative Fraction
        if ((eta+shtflx).ne.0) then 
           call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_EF,                    &
                value=abs(eta/(eta+shtflx)),vlevel=1,unit="-",direction="-",&
                valid_min=0.0,valid_max=10.0,&
                surface_type = LIS_rc%lsm_index)
        else
           !double check
           call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_EF,value=1.0,          &
                vlevel=1,unit="-",direction="-",&
                valid_min=0.0,valid_max=10.0,&
                surface_type = LIS_rc%lsm_index)
        endif

        if (sfctmp .le. LIS_CONST_TKFRZ) then
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=prcp,       &
                vlevel=1,unit="kg m-2 s-1",direction="DN",&
                valid_min=0.0,valid_max=0.0085,&
                surface_type = LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=0.0,        &
                vlevel=1,unit="kg m-2 s-1",direction="DN",&
                valid_min=0.0,valid_max=0.02,&
                surface_type = LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=prcp*dt,    &
                vlevel=1,unit="kg m-2",direction="DN",&
                surface_type = LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=0.0,        &
                vlevel=1,unit="kg m-2",direction="DN",&
                surface_type = LIS_rc%lsm_index)
        else
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=0.0,        &
                vlevel=1,unit="kg m-2 s-1",direction="DN",&
                valid_min=0.0,valid_max=0.0085,&
                surface_type = LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=prcp,       &
                vlevel=1,unit="kg m-2 s-1",direction="DN",&
                valid_min=0.0,valid_max=0.02,&
                surface_type = LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=0.0,        &
                vlevel=1,unit="kg m-2",direction="DN",&
                surface_type = LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=prcp*dt,    &
                vlevel=1,unit="kg m-2",direction="DN",&
                surface_type = LIS_rc%lsm_index)
        endif
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TOTALPRECIP,               &
             value=prcp,vlevel=1,unit="kg m-2 s-1",direction="DN",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TOTALPRECIP,               &
             value=prcp*dt,vlevel=1,unit="kg m-2",direction="DN",&
             surface_type = LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EVAP,value=evp,            &
             vlevel=1,unit="kg m-2 s-1",direction="UP",&
             valid_min=-0.0003,valid_max=0.0003,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EVAP, value=(evp*3600.),   &
             vlevel=1,unit="mm hr-1",direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QS,                        &
             value= noah271_struc(n)%noah(t)%runoff1*LIS_CONST_RHOFW,vlevel=1,unit="kg m-2 s-1",&
             direction="OUT",valid_min=0.0,valid_max=5.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QS,                        &
             value= noah271_struc(n)%noah(t)%runoff1*LIS_CONST_RHOFW*dt,vlevel=1,unit="kg m-2",&
             direction="OUT",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSB,                       &
             value= noah271_struc(n)%noah(t)%runoff2*LIS_CONST_RHOFW,vlevel=1,unit="kg m-2 s-1",&
             direction="OUT",valid_min=0.0,valid_max=5.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSB,                       &
             value= noah271_struc(n)%noah(t)%runoff2*LIS_CONST_RHOFW*dt,vlevel=1,&
             unit="kg m-2",direction="OUT",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSM,                       &
             value= snomlt*LIS_CONST_RHOFW/dt,vlevel=1,unit="kg m-2 s-1",&
             direction="S2L",valid_min=0.0,valid_max=0.005,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSM,                       &
             value= snomlt*LIS_CONST_RHOFW,vlevel=1,unit="kg m-2",direction="S2L",&
             surface_type = LIS_rc%lsm_index)
        soilmtc = 0
        do i=1,noah271_struc(n)%nslay
           soilmtc = soilmtc + noah271_struc(n)%noah(t)%soilmoist(i)
        enddo
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSOILMOIST,              &
             value=(soilmtc-startsm),vlevel=1,unit="kg m-2",&
             direction="INC",valid_min=-2000.0,valid_max=2000.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSWE,                    &
             value=(noah271_struc(n)%noah(t)%sneqv*LIS_CONST_RHOFW-       &
             startswe),vlevel=1,unit="kg m-2",direction="INC",&
             valid_min=-2000.0,valid_max=2000.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELINTERCEPT,              &
             value=(noah271_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW-startint),&
             vlevel=1,unit="kg m-2",direction="INC",&
             valid_min=-100.0,valid_max=100.0,&
             surface_type = LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_BARESOILT,value=TSOIL,     &
             vlevel=1,unit="K",direction="-",&
             valid_min=213.0, valid_max=280.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_AVGSURFT,value=(TSOIL *    &
             (1.0 - SNCOVR)) + (NOAH271_STRUC(N)%NOAH(T)%T1 * SNCOVR),    &
             vlevel=1,unit="K",direction="-",&
             valid_min=213.0, valid_max=333.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RADT,                      &
             value=NOAH271_STRUC(N)%NOAH(T)%T1,vlevel=1,unit="K",&
             direction="-",valid_min=213.0, valid_max=353.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ALBEDO,                    &
             value=ALBEDO,vlevel=1,unit="-",direction="-",&
             valid_min=0.0, valid_max=1.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ALBEDO,                    &
             value=(ALBEDO*100.),vlevel=1,unit="%",direction="-",&
             valid_min=0.0, valid_max=100.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWE,vlevel=1,unit="kg m-2", &
             value=noah271_struc(n)%noah(t)%sneqv*LIS_CONST_RHOFW,&
             direction="-",valid_min=0.0,valid_max=2000.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWE,vlevel=1,unit="m",     &
             value=noah271_struc(n)%noah(t)%sneqv,direction="-",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWDEPTH,                 &
             vlevel=1,value = noah271_struc(n)%noah(t)%snowh,unit="m",&
             direction="-",&
             surface_type = LIS_rc%lsm_index)

        if(noah271_struc(n)%noah(t)%snowh.gt.0) then 
           bdsno = (noah271_struc(n)%noah(t)%sneqv*1000)/(noah271_struc(n)%noah(t)%snowh)
        else
           bdsno = LIS_rc%udef
        endif

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWDENSITY,                 &
             vlevel=1,value = bdsno,unit="kg m-3",&
             direction="-",&
             surface_type = LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWCOVER,value=SNCOVR,    &
             vlevel=1,unit="-",direction="-",valid_min=0.0,&
             valid_max=1.0,&
             surface_type = LIS_rc%lsm_index)

        ! NOTE:  Soil temperature for each layer is passed directly to output routines.
        do i=1,noah271_struc(n)%nslay
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILMOIST,vlevel=i,     &
                value=noah271_struc(n)%noah(t)%smc(i)*sldpth(i)*          &
                LIS_CONST_RHOFW,unit='kg m-2',direction="-",&
                valid_min=0.0,valid_max=2000.0,&
                surface_type = LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILMOIST,vlevel=i,     &
                value=noah271_struc(n)%noah(t)%smc(i), &
                unit='m^3 m-3',direction="-",valid_min=0.0,valid_max=0.55,&
                surface_type = LIS_rc%lsm_index)
        enddo
        do i=1,noah271_struc(n)%nslay
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILTEMP,vlevel=i,      &
                value=noah271_struc(n)%noah(t)%stc(i),unit="K",&
                direction="-",valid_min=213.0,valid_max=333.0,&
                surface_type = LIS_rc%lsm_index)
        enddo
        do i=1,noah271_struc(n)%nslay
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SMLIQFRAC,vlevel=i,     &
                value=(noah271_struc(n)%noah(t)%sh2o(i)/                  &
                noah271_struc(n)%noah(t)%smc(i)),unit='-',&
                direction="-",valid_min=0.0,valid_max=1.0,&
                surface_type = LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SMLIQFRAC,vlevel=i,     &
                value=noah271_struc(n)%noah(t)%sh2o(i),unit='m^3 m-3',direction="-",&
                surface_type = LIS_rc%lsm_index)
        enddo
        do i=1,noah271_struc(n)%nslay
           temp = &
                (noah271_struc(n)%noah(t)%smc(i) - noah271_struc(n)%noah(t)%smcwlt) / &
                (noah271_struc(n)%noah(t)%smcmax - noah271_struc(n)%noah(t)%smcwlt)

           if ( temp > 1.0 ) then
              temp  = 1.0
           endif
           if ( temp  < 0.01 ) then
              temp  = 0.01
           endif

           noah271_struc(n)%noah(t)%relsmc(i) = temp

           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RELSMC,vlevel=i,        &
                value=temp, unit='-',direction="-",&
                surface_type = LIS_rc%lsm_index)
        enddo
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILWET,value=MSTAVTOT,    &
             vlevel=1,unit="-",direction="-",&
             surface_type = LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP,                   &
             value=(etp/(((1.-SNCOVR)*LIS_CONST_LATVAP) +                 &
             SNCOVR*CONST_LATSUB)),vlevel=1,unit="kg m-2 s-1",                &
             direction="UP",valid_min=-0.0006,valid_max=0.0006,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP,                   &
             value=((etp/(((1.-SNCOVR)*LIS_CONST_LATVAP) +                &
             SNCOVR*CONST_LATSUB))*3600.),vlevel=1,unit="mm hr-1 ",         &
             direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP,                   &
             value=etp,vlevel=1,unit="W m-2",direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ECANOP,                    &
             value=(ec/LIS_CONST_LATVAP),vlevel=1,unit="kg m-2 s-1",&
             direction="UP",valid_min=-0.0003,valid_max=0.0003,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ECANOP,                    &
             value=((ec/LIS_CONST_LATVAP)*3600.),vlevel=1,&
             unit="mm hr-1",direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ECANOP,                    &
             value=ec,vlevel=1,unit="W m-2",direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TVEG,                      &
             value=(ett/LIS_CONST_LATVAP),vlevel=1,unit="kg m-2 s-1",&
             direction="UP",valid_min=-0.0003,valid_max=0.0003,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TVEG,                      &
             value=((ett/LIS_CONST_LATVAP)*3600.),vlevel=1,unit="mm hr-1",direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TVEG,                      &
             value=ett,vlevel=1,unit="W m-2",direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ESOIL,                     &
             value=(edir/LIS_CONST_LATVAP),vlevel=1,unit="kg m-2 s-1",&
             direction="UP",valid_min=-0.0003,valid_max=0.0003,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ESOIL,                     &
             value=((edir/LIS_CONST_LATVAP)*3600.),vlevel=1,&
             unit="mm hr-1",direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ESOIL,                     &
             value=edir,vlevel=1,unit="W m-2",direction="UP",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SUBSNOW,                   &
             value=(esnow/2.83E+6),vlevel=1,unit="kg m-2 s-1",direction="-",&
             valid_min=-0.0003,valid_max=0.0003,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SUBSNOW,                   &
             value=((esnow/2.83E+6)*3600.0),vlevel=1,unit="mm hr-1",direction="-",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SUBSNOW,                   &
             value=esnow,vlevel=1,unit="W m-2",direction="-",&
             surface_type = LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CANOPINT,vlevel=1,         &
             value=noah271_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW,          &
             unit="kg m-2",direction="-",valid_min=0.0,valid_max=100.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROOTMOIST,vlevel=1,        &
             value=soilrz,unit="kg m-2",direction="-",valid_min=0.0,&
             valid_max=2000.0,&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROOTMOIST,vlevel=1,        &
             value=soilrz/1000.0,unit="m^3 m-3",direction="-",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROOTTEMP,vlevel=1,         &
             value=roottemp,unit="K",direction="-",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RHMIN,vlevel=1,            &
             value=noah271_struc(n)%noah(t)%rhmin,unit="-",direction="-",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RHMIN,vlevel=1,            &
             value=(noah271_struc(n)%noah(t)%rhmin*100.),unit="%",direction="-",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EMISSFORC,                 &
             value=EMISSI_OUT,vlevel=1,unit="-",direction="-",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GREENNESS,                 &
             value=shdfac,vlevel=1,unit="-",direction="-",&
             surface_type = LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GREENNESS,                 &
             value=shdfac*100.0,vlevel=1,unit="%",direction="-",&
             surface_type = LIS_rc%lsm_index)

! EMK Bug fix...Output deep soil temperature                                
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TEMPBOT,                   &
          value=noah271_struc(n)%noah(t)%tempbot,vlevel=1,unit="K",direction="-",&
          surface_type = LIS_rc%lsm_index)

     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROUGHNESS,                 &
          value=z0,vlevel=1,unit="m",direction="-",&
          surface_type = LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LAI,                       &
          value=xlai,vlevel=1,unit="-",direction="-",&
          surface_type = LIS_rc%lsm_index)
     
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CH,                        &
          value=noah271_struc(n)%noah(t)%CH,vlevel=1,unit="m s-1",direction="-",&
          surface_type = LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CM,                        &
          value=noah271_struc(n)%noah(t)%CM,vlevel=1,unit="m s-1",direction="-",&
          surface_type = LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_MIXRATIO,                  &
          value=Q1,vlevel=1,unit="kg kg-1",direction="-",&
          surface_type = LIS_rc%lsm_index)

     ! SY: Begin FLDAS
     WR_TimeStep = (etp/(((1.-SNCOVR)*LIS_CONST_LATVAP) +              &
          SNCOVR*CONST_LATSUB))*LIS_rc%ts
     AET_TimeStep  = evp*LIS_rc%ts
     WRSI_TimeStep = LIS_rc%udef
     if( (WR_TimeStep  > SUMWR_EXITWRSIBELOWTHISQTY) .and. &
         (WR_TimeStep  < SUMWR_EXITWRSIABOVETHISQTY) .and. &
         (AET_TimeStep < SUMET_EXITWRSIABOVETHISQTY) ) then
       if (AET_TimeStep > WR_TimeStep) then
         WRSI_TimeStep = 100
       else
         WRSI_TimeStep = 100*AET_TimeStep/WR_TimeStep
       endif
     endif
     IF( WR_TimeStep <= 0.0 ) THEN
       WRSI_TimeStep = 0.0
       AET_TimeStep = WR_TimeStep
       IF ( WR_TimeStep < 0.0 ) THEN
          WRSI_TimeStep = 100
       END IF
     END IF
     IF ( 100*soilrz/soilrzmax < noah271_struc(n)%noah(t)%smcwlt ) THEN
        WRSI_TimeStep = 1
     END IF
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WR_TimeStep,               &
          value=WR_TimeStep,vlevel=1,unit="kg m-2",direction="-",&
          surface_type = LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_AET_TimeStep,              &
          value=AET_TimeStep,vlevel=1,unit="kg m-2",direction="-",&
          surface_type = LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WRSI_TimeStep,              &
          value=WRSI_TimeStep,vlevel=1,unit="-",direction="-",&
          surface_type = LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SurplusWater_TimeStep,     &
          value= noah271_struc(n)%noah(t)%runoff1*LIS_CONST_RHOFW*dt +  &
          noah271_struc(n)%noah(t)%runoff2*LIS_CONST_RHOFW*dt,   &
          vlevel=1,unit="kg m-2",direction="-",&
          surface_type = LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWI,vlevel=1,        &
          value=100*soilrz/soilrzmax,unit="%",direction="-",&
          surface_type = LIS_rc%lsm_index)
     ! SY: End FLDAS

        deallocate(sldpth)
!Adjust deep soil temperature as lagged function of skin temperature
    if (LIS_rc%tbot_update_lag .eq. 1) then
!       write(LIS_logunit,*)'EMK: domain,tile,sncovr,tsoil,t1 = ', &
!            n,t,sncovr,tsoil,NOAH32_STRUC(N)%NOAH(T)%T1
       call LIS_updateTbot(n,t,julian_in,yr,dt,&
            (TSOIL * (1.0 - SNCOVR)) + (NOAH271_STRUC(N)%NOAH(T)%T1 * SNCOVR), &
            tbot)
       noah271_struc(n)%noah(t)%tempbot = tbot
    end if
!reset 
        noah271_struc(n)%noah(t)%tair = 0 
        noah271_struc(n)%noah(t)%qair = 0 
        noah271_struc(n)%noah(t)%swdown = 0 
        noah271_struc(n)%noah(t)%lwdown = 0
        noah271_struc(n)%noah(t)%uwind = 0
        noah271_struc(n)%noah(t)%vwind = 0 
        noah271_struc(n)%noah(t)%psurf = 0 
        noah271_struc(n)%noah(t)%rainf = 0 
        noah271_struc(n)%noah(t)%rainf_c = 0 
        
        if(noah271_struc(n)%forcing_ch .ne.0) then 
           noah271_struc(n)%noah(t)%ch = 0 
        endif


     enddo
     noah271_struc(n)%forc_count = 0 
     
  endif
end subroutine noah271_main


!*** NOAH FUNCTIONS ****************************************************

FUNCTION noah271_DQS (T)
  
  use LIS_constantsMod, only : LIS_CONST_LATVAP, LIS_CONST_TKFRZ
  IMPLICIT NONE
  
!
!  PURPOSE:  TO CALCULATE VALUES OF VAPOR PRESSURE (E)
!            AND P * DQS/DT (P TIMES CHG IN SAT MXG RATIO WITH RESPECT
!            TO THE CHG IN TEMP) IN SUBSTITUTION TO THE LOOK-UP TABLES.
!
!            FORMULAS AND LIS_CONSTANTS FROM ROGERS AND YAU, 1989.
!                         ADDED BY PABLO J. GRUNMANN, 6/30/97.
!

  REAL DESDT
  REAL noah271_DQS
!    REAL ESD
  REAL LW
  REAL T
  REAL ES
  
  !      REAL, PARAMETER:: CP = 1005.
  !      REAL, PARAMETER:: CV = 718.
  !      REAL, PARAMETER:: CVV = 1410.
  REAL, PARAMETER:: CPV = 1870.
  REAL, PARAMETER:: RV = 461.5
  REAL, PARAMETER:: CW = 4187.
  REAL, PARAMETER:: EPS = 0.622
  REAL, PARAMETER:: ESO = 611.2
  REAL, PARAMETER:: TO = 273.15
  REAL, PARAMETER:: LVH2O = 2.501000E+6
  
  
!     ABOUT THE PARAMETERS:
!
!     EPS ---------- WATER - DRY AIR MOLECULAR MASS RATIO, EPSILON
!
!   VALUES FOR SPECIFIC HEAT CAPACITY AND INDIVIDUAL GAS LIS_CONSTANTS
!   IN [JOULES/(KG*KELVIN)] UNITS.
!
!     DRY AIR:
!             CP, CV
!     WATER VAPOR:
!                 CVV = 1410.
!                 CPV = 1870.
!                 RV  =  461.5
!     LIQUID WATER:
!                  CW = 4187.
!
!     ESO = ES(T=273.15 K) = SAT. VAPOR PRESSURE (IN PASCAL) AT T=TO
!      TO = 273.15
!
!     SAT. MIXING  RATIO: QS ~= EPS*ES/P
!     CLAUSIUS-CLAPEYRON: DES/DT = L*ES/(RV*T^2)
!     @QS/@T =  (EPS/P)*DES/DT

  LW = LIS_CONST_LATVAP - ( CW - CPV ) * ( T - LIS_CONST_TKFRZ )
  ES = ESO*EXP (LW*(1/LIS_CONST_TKFRZ - 1/T)/RV)
  DESDT = LW*ES/(RV*T*T)
  
  !    FOR INSERTION IN DQSDT FUNCTION:
  !    DQSDT = DQS/P , WHERE DQS = EPS*DESDT
  
  noah271_DQS = EPS*DESDT
  
  RETURN
END FUNCTION NOAH271_DQS

!----------------------------------------------------------------------

FUNCTION DQSDT ( SFCTMP, SFCPRS )

  IMPLICIT NONE

!
!    PURPOSE:  TO RETRIEVE THE APPROPRIATE VALUE OF DQSDT (THE CHANGE
!    =======   OF THE SATURATION MIXING RATIO WITH RESPECT TO THE
!              CHANGE IN TEMPERATURE) FROM:
!
!              FORMULAS INTRODUCED IN FUNCTION DQS
!                                  (MODIFIED BY PABLO GRUNMANN, 7/9/97).
!

  REAL SFCTMP
  REAL SFCPRS
  REAL NOAH271_DQS
  REAL DQSDT
  
  IF ((SFCTMP .GE. 173.0) .AND. (SFCTMP  .LE.  373.0)) THEN

!  IF THE INPUT SFC AIR TEMP IS BTWN 173 K AND 373 K, USE
!   FUNCTION DQS TO DETERMINE THE SLOPE OF SAT.MIX RATIO FUNCTION

     DQSDT = NOAH271_DQS (SFCTMP) / SFCPRS

  ELSE

!  OTHERWISE, SET DQSDT EQUAL TO ZERO

     DQSDT = 0.0
     
  END IF
  
  RETURN
END FUNCTION DQSDT

!---------------------------------------------------------------

FUNCTION E(T)
  
  use LIS_constantsMod, only : LIS_CONST_LATVAP, LIS_CONST_TKFRZ
  IMPLICIT NONE

!
!  PURPOSE:  TO CALCULATE VALUES OF SAT. VAPOR PRESSURE (E)
!            FORMULAS AND LIS_CONSTANTS FROM ROGERS AND YAU, 1989.
!                         ADDED BY PABLO J. GRUNMANN, 7/9/97.
!

  REAL LW
  REAL T
  REAL E

!      REAL, PARAMETER:: EPS = 0.622 
!      REAL, PARAMETER:: CP = 1005.
!      REAL, PARAMETER:: CV = 718.
!      REAL, PARAMETER:: CVV = 1410.
  REAL, PARAMETER:: CPV = 1870.
  REAL, PARAMETER:: RV = 461.5
  REAL, PARAMETER:: CW = 4187.
  REAL, PARAMETER:: ESO = 611.2
  REAL, PARAMETER:: TO = 273.15
  REAL, PARAMETER:: LVH2O = 2.501000E+6
  
!   ABOUT THE PARAMETERS:
!
!    EPS --- WATER - DRY AIR MOLECULAR MASS RATIO, EPSILON
!
!    VALUES FOR SPECIFIC HEAT CAPACITY AND INDIVIDUAL GAS LIS_CONSTANTS
!    IN [JOULES/(KG*KELVIN)] UNITS.
!
!     DRY AIR:
!             CP, CV
!     WATER VAPOR:
!             CVV = 1410.
!             CPV = 1870.
!             RV  =  461.5
!     LIQUID WATER:
!             CW = 4187.
!
!     ESO = ES(TO) = SAT. VAPOR PRESSURE (IN PASCAL) AT T=TO
!      TO = 273.15
!
!     CLAUSIUS-CLAPEYRON: DES/DT = L*ES/(RV*T^2)

  LW = LIS_CONST_LATVAP - ( CW - CPV ) * ( T - LIS_CONST_TKFRZ )
  E = ESO*EXP (LW*(1/LIS_CONST_TKFRZ - 1/T)/RV)

  RETURN
END FUNCTION E


!CC 1. DRIVER SUBROUTINE ==> SUBROUTINE OBTLWDN CCCCCCCCCCCCCCCCC

!      SUBROUTINE OBTLWDN(SFCTMP,LWDN)

!                      RADIATION
!
! The following step (OBTENTION OF LWDN) is used if
! user wants to calculate longwave downward.
!
! OBTENTION OF LWDN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! COMPUTATION OF LWDN (INCOMING LW RADIATION) FROM TAIR AND Q:
!
!...............LWDN = EMISS*SIGMA*(TAK)^4.
!
! WHERE:   TAK = AIR TEMP IN KELVIN
!        EMISS = 0.7  OR  (IDSO AND JACKSON, 1969):
!
!        EMISS = (1 - 0.261 EXP(-7.77*10^(-4)X(273-TAK)^2)
!
!      NEED STEFAN-BOLTZMANN LIS_CONSTANT, SIGMA
!         SIGMA = 5.672 * 10^-8  W M^-2 T^-4
!
!           SIGMA = 5.672E-8
!           TAK = SFCTMP
!           EMISS = 1 - 0.261*EXP((-7.77E-4)*(273-TAK)^2.)
!
!           LWDN = EMISS*SIGMA*TAK^4.
!
!        RETURN
!        END

!CCCC  END OF DRIVER SUBROUTINES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!BOP
!
! !ROUTINE: sh2o_init.F90
!
! !DESCRIPTION:
!  To do a 'cold start' initialization of soil liquid water SH2O for 
!   NOAH LSM, though also adaptable for other land-surface models,
!   using either GDAS or Eta forcing data. 
!
! !REVISION HISTORY:
!  NCEP;  Original code developed by NCEP for Eta model
!           subroutines to initialize soil liquid water, SH2O   
!  04 Nov 2002: Kristi Arsenault; Modified code to be used with noahrst.f
!                                 to initialize NOAH with NCEP forcing data  
!
!
! !INTERFACE:
subroutine sh2oinit(smc,stc,smcmax,psis,beta,sh2o)

  implicit none      
! !ARGUMENTS:
  REAL STC      ! NOAH Soil Layer Temperature (K)
  REAL SMC      ! NOAH Soil Layer Total Moisture (liq+frzn) 
  REAL SH2O     ! NOAH Soil Layer Liquid Moisture
  
  REAL  PSIS                    ! Saturated soil potential
  REAL  BETA                    ! B-parameter 
  REAL  SMCMAX                  ! Max soil moisture content (porosity)
  REAL  BX
!EOP
  REAL  FK
  REAL  FRH2O 

  REAL :: HLICE=3.335E5         ! Ice parameter
  REAL :: GRAV=9.81             ! Gravity (m s-1) 
  REAL :: T0=273.15             ! Freezing point of water
  REAL :: BLIM=5.5              ! B-parameter upper limit

!=== End Variable Definition =============================================
! ----------------------------------------------------------------------
! COLD START:  determine liquid soil water content (SH2O)
! NSOIL number of soil layers
! ----------------------------------------------------------------------
!  SH2O <= SMC for T < 273.149K (-0.001C)
  IF (STC .LT. 273.149) THEN
! ----------------------------------------------------------------------
! first guess following explicit solution for Flerchinger Eqn from Koren 
! et al, JGR, 1999, Eqn 17 (KCOUNT=0 in FUNCTION FRH2O). 
! ----------------------------------------------------------------------
     BX = BETA
     IF ( BETA .GT. BLIM )  BX = BLIM
     FK=(((HLICE/(GRAV*(-PSIS)))* & 
          ((STC-T0)/STC))**(-1/BX))*SMCMAX
     IF (FK .LT. 0.02) FK = 0.02
     SH2O = MIN ( FK, SMC )
! ----------------------------------------------------------------------
! now use iterative solution for liquid soil water content using 
! FUNCTION FRH2O with the initial guess for SH2O from above explicit
! first guess.
! ----------------------------------------------------------------------
     SH2O = FRH2O(STC,SMC,SH2O,SMCMAX,BETA,PSIS)

  ELSE
! ----------------------------------------------------------------------
!  SH2O = SMC for T => 273.149K (-0.001C)
     SH2O=SMC
! ----------------------------------------------------------------------
  ENDIF
  
  RETURN 
END subroutine sh2oinit

!===========================================================================

FUNCTION FRH2O(TKELV,SMC,SH2O,SMCMAX,BEXP,PSIS)

  IMPLICIT NONE
  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  PURPOSE:  CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL WATER CONTENT 
!  IF TEMPERATURE IS BELOW 273.15K (T0).  REQUIRES NEWTON-TYPE ITERATION 
!  TO SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF 
!  KOREN ET AL. (1999, JGR, VOL 104(D16), 19569-19585).
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! New version (JUNE 2001): much faster and more accurate newton iteration 
! achieved by first taking log of eqn cited above -- less than 4
! (typically 1 or 2) iterations achieves convergence.  Also, explicit  
! 1-step solution option for special case of parameter Ck=0, which reduces 
! the original implicit equation to a simpler explicit form, known as the 
! ""Flerchinger Eqn". Improved handling of solution in the limit of  
! freezing point temperature T0. 
! 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! INPUT:
!   TKELV.........Temperature (Kelvin)
!   SMC...........Total soil moisture content (volumetric)
!   SH2O..........Liquid soil moisture content (volumetric)
!   SMCMAX........Saturation soil moisture content (from REDPRM)
!   B.............Soil type "B" parameter (from REDPRM)
!   PSIS..........Saturated soil matric potential (from REDPRM)
!
! OUTPUT:
!   FRH2O.........supercooled liquid water content.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       REAL BEXP
       REAL BLIM
       REAL BX
       REAL CK
       REAL DENOM
       REAL DF
!       REAL DH2O
!       REAL DICE
       REAL DSWL
       REAL ERROR
       REAL FK
       REAL FRH2O
       REAL GS
       REAL HLICE
       REAL PSIS
       REAL SH2O
       REAL SMC
       REAL SMCMAX
       REAL SWL
       REAL SWLK
       REAL TKELV
       REAL T0

       INTEGER NLOG
       INTEGER KCOUNT

       PARAMETER (CK=8.0)
!      PARAMETER (CK=0.0)
       PARAMETER (BLIM=5.5)
!      PARAMETER (BLIM=7.0)
       PARAMETER (ERROR=0.005)

       PARAMETER (HLICE=3.335E5)
       PARAMETER (GS = 9.81)
!       PARAMETER (DICE=920.0)
!       PARAMETER (DH2O=1000.0)
       PARAMETER (T0=273.15)

!  ###   LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)  ####
!  ###   SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT  ####
!  ###   IS NON-REALISTICALLY HIGH AT VERY LOW TEMPERATURES    ####
!##################################################################

      BX = BEXP
      IF ( BEXP .GT. BLIM ) BX = BLIM
!------------------------------------------------------------------

! INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
      NLOG=0
      KCOUNT=0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (T0), SH2O = SMC 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      IF (TKELV .GT. (T0 - 1.E-3)) THEN
         FRH2O=SMC

      ELSE

        IF (CK .NE. 0.0) THEN

! -------------------------------------------------------------
! OPTION 1: ITERATED SOLUTION FOR NONZERO CK
! IN KOREN ET AL, JGR, 1999, EQN 17
! -------------------------------------------------------------
! INITIAL GUESS FOR SWL (frozen content)
         SWL = SMC-SH2O

! KEEP WITHIN BOUNDS.
          IF (SWL .GT. (SMC-0.02)) SWL=SMC-0.02
          IF(SWL .LT. 0.) SWL=0. 
!--------------------------------------------------------------
!  START OF ITERATIONS 
!--------------------------------------------------------------
         DO WHILE (NLOG .LT. 10 .AND. KCOUNT .EQ. 0)
          NLOG = NLOG+1
          DF = ALOG(( PSIS*GS/HLICE ) * ( ( 1.+CK*SWL )**2. ) * & 
               ( SMCMAX/(SMC-SWL) )**BX) - ALOG(-(TKELV-T0)/TKELV)
          DENOM = 2. * CK / ( 1.+CK*SWL ) + BX / ( SMC - SWL )
          SWLK = SWL - DF/DENOM
! BOUNDS USEFUL FOR MATHEMATICAL SOLUTION.
          IF (SWLK .GT. (SMC-0.02)) SWLK = SMC - 0.02
          IF(SWLK .LT. 0.) SWLK = 0.
! MATHEMATICAL SOLUTION BOUNDS APPLIED.
          DSWL=ABS(SWLK-SWL)
          SWL=SWLK 

!---------------------------------------------------------------
! IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)  
! WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED. 
!---------------------------------------------------------------
          IF ( DSWL .LE. ERROR )  THEN
            KCOUNT=KCOUNT+1
          END IF
         END DO 
!---------------------------------------------------------------
!  END OF ITERATIONS 
!---------------------------------------------------------------
! BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION.
         FRH2O = SMC - SWL

!CCCCCCCCCCCCCCCCCCCCCCC END OPTION 1 CCCCCCCCCCCCCCCCCCCCCCCCCCC

        ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0
! IN KOREN ET AL., JGR, 1999, EQN 17
!----------------------------------------------------------------
        IF (KCOUNT .EQ. 0) THEN
!      Print*,'Flerchinger used in NEW version. Iterations=',NLOG
          FK=(((HLICE/(GS*(-PSIS)))*((TKELV-T0)/TKELV))** & 
           (-1/BX))*SMCMAX
!  APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION
          IF (FK .LT. 0.02) FK = 0.02
          FRH2O = MIN ( FK, SMC )

!CCCCCCCCCCCCCCCCCCCCCCCCC END OPTION 2 CCCCCCCCCCCCCCCCCCCCCCCCCC

        ENDIF

      ENDIF

      RETURN
      END
