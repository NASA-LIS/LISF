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
! !ROUTINE: noah36_main
! \label{noah36_main}
! 
! !REVISION HISTORY:
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!  20 Jan 2011: David Mocko, fixes to output diagnostics for energy
!                            & water balance and fix to conversion of
!                            U,V forcing height to same level as T,q
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  14 Jan 2014: David Mocko, reconfirmed Noah3.3 in LIS7.0
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!  05 Jan 2021: Augusto Getirana: 2-way coupling
!
! !INTERFACE:
subroutine noah36_main(n)
! !USES:
!$ use omp_lib
  use LIS_coreMod
  use LIS_timeMgrMod,    only : LIS_isAlarmRinging
  use LIS_albedoMod,     only : LIS_alb
  use LIS_constantsMod,  only : LIS_CONST_RHOFW, LIS_CONST_TKFRZ, &
                                LIS_CONST_LATVAP
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use LIS_histDataMod
  use LIS_FORC_AttributesMod 
  use module_sf_noah36lsm, only : SFLX, SFCDIF_OFF, CALHUM, SNFRAC
  use module_sfcdif_wrf_36, only: SFCDIF_MYJ
  use module_sf_noah36lsm_glacial, only : SFLX_glacial
  use noah36_lsmMod
  use LIS_tbotAdjustMod, only: LIS_tbotTimeUtil,LIS_updateTbot

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n

!
! !DESCRIPTION:
!  This is the entry point for calling the Noah-3.6 LSM physics.
!  This routine calls the {\tt SFLX} routine that performs the
!  land surface computations, to solve for water and energy equations.
!  For documentation of the {\tt SFLX} and routines from some
!  previous versions of Noah, please see:
!      http://www.ral.ucar.edu/research/land/technology/lsm.php
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  integer, parameter :: nsold = 20 !maximum number of soil layers
  real        :: ffrozp     ! flag used for precipitation type
  real        :: dt         ! time step length (sec) of physical integration
  real        :: zlvl       ! observation height (m) for T,q
  real        :: zlvl_wind  ! observation height (m) for wind (u,v)
!  REAL        :: ZFAC       ! Factor to convert U,V to same level as T,q  
  real        :: sfctmp     ! air temperature at zlvl above ground (K)
  real        :: q2         ! mixing ration at height zlvl above ground (kg kg-1)
  real        :: sfcprs     ! PRESSURE AT HEIGHT ZLVL ABOVE GROUND (PASCALS)
  real        :: prcp       ! PRECIP RATE (KG M-2 S-1) (NOTE, THIS IS A RATE)
  real        :: uwind      ! U-Wind component (m s-1)
  real        :: vwind      ! V-Wind component (m s-1)
  real        :: cpcp       ! Convective Precipitation (kg m-2 s-1)
  real        :: soldn      ! SOLAR DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET SOLAR)
  real        :: q2sat      ! SAT SPECIFIC HUMIDITY AT HEIGHT ZLVL ABOVE GROUND (KG KG-1)
  real        :: lwdn       ! LW DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET LONGWAVE)
  real        :: sfcspd     ! WIND SPEED (M S-1) AT HEIGHT ZLVL_WIND ABOVE GROUND
  real        :: esat       ! Saturation vapor pressure for water (Pa)
  integer     :: nsoil      ! number of soil layers
  real, allocatable :: sldpth(:)  ! thickness values for soil layers(meters)
  integer     :: ice        ! flag for ice physics (0=land, 1=sea-ice, -1=glacier-ice)
  integer     :: isurban    !
  real        :: solnet     ! NET DOWNWARD SOLAR RADIATION ((W M-2; POSITIVE)
!  real        :: alb        ! Background Snow-free albedo (fraction)
!  real        :: albedo     ! Total surface albedo fraction (incl. any snow cover effect)
  logical     :: local 
  real        :: solardirect ! Direct component of downward solar radiation (W M-2) (not used)
  real        :: prcprain  ! Liquid-precipitation rate (KG M-2 S-1) (not used)
  real        :: cosz      ! Solar zenith angle (not used for now)
  real        :: th2       ! Potential temperature (K) at hgt z above grnd
  real        :: th2v      ! Virtual potential temperature (K) at z
  real        :: t1v       ! Virtual temperature at ground (sub 1)
  real        :: t2v       ! Virtual temp. at 1st mid. lev. above grnd (2)
  real        :: dqsdt2    ! Slope of sat. specific hum. curve for Penman
  real        :: slope     ! slope estimate of linear reservoir coefficient  
  real        :: shdfac    ! 12-month green vegetation fraction
  real        :: shdmin    ! minimum areal fractional coverage of annual green vegetation
  real        :: shdmax    ! maximum areal fractional coverage of annual green vegetation
  real        :: tbot      ! annually-fixed, soil-temp condition at zbot
  real        :: evp
  real        :: eta       ! actual latent heat flux (w m-2) (neg: up from sfc)
  real        :: eta_kinematic
  real        :: shtflx    ! Sensible heat flux (W m-2) 
  real        :: fdown     ! Radiation forcing at the surface
  real        :: ec        ! canopy evaporation (m s-1)
  real        :: edir      ! direct soil evaporation (m s-1)
  real        :: et(nsold) ! plant transp. from each root/soil layer(m s-1)
  real        :: SMAV(nsold) ! Soil Moisture Availability at each level, fraction between
                             ! SMCWLT (SMAV=0.0) and SMCMAX (SMAV=1.0)
  real        :: ett       ! accum plant transpiration (m s-1)
  real        :: esnow     ! sublimation from snowpack (kg m-2 s-1) (or dep.)
  real        :: drip      ! excess canopy moisture (m)
  real        :: dew       ! dewfall amount  (m s-1)
  real        :: beta      ! ratio of actual/potential evap (dimensionless) 
  real        :: etp       ! final potential evapotransp. (w m-2)
  real        :: gflx      ! soil heat flux (w m-2)
  real        :: flx1      ! precip-snow surface  (w m-2)
  real        :: flx2      ! freezing rain latent heat flux (w m-2)
  real        :: flx3      ! phase-change heat flux from snowmelt (w m-2)
  real        :: flx4      ! UA snow-physics energy added to sensible heat flux from canopy (w m-2)
  real        :: fvb       ! UA snow-physics fraction of vegetation with snow below (unitless fraction, 0-1)
  real        :: gama      ! UA snow-physics exp(-XLAI)
  real        :: fbur      ! UA snow-physics fraction of vegetation covered by snow (unitless fraction, 0-1)
  real        :: fgsn      ! UA snow-physics fraction of ground covered by snow (unitless fraction, 0-1)
! The following three variables are for WRF-HYDRO,
!   which is not currently implemented within LIS:
  real        :: sfhead1rt, infxs1rt, etpnd1
  real        :: snomlt    ! snow melt (m) (water equivalent)
  real        :: sncovr    ! fractional snow cover (unitless fraction, 0-1)
!  real        :: runoff1   ! ground surface runoff (m s-1)
!  real        :: runoff2   ! underground runoff (m s-1)
  real        :: runoff3   ! runoff within soil layers (m s-1)
  real        :: rc        ! canopy resistance (s m-1)
  real        :: pc        ! plant coefficient (unitless fraction, 0-1) where
  real        :: rcs       ! incoming solar rc factor (dimensionless)
  real        :: rct       ! air temperature rc factor (dimensionless)
  real        :: rcq       ! atmos vapor pressure deficit rc factor (dimensionless)
  real        :: rcsoil    ! soil moisture rc factor (dimensionless)
  real        :: soilw     ! avail root zone soil moisture (unitless fraction)
  real        :: soilm     ! total soil column moisture content (frozen+unfrozen)(m)
  REAL        :: TSOIL     ! SOIL SURFACE TEMPERATURE (K)
!  real        :: q1        ! Effective mixing ratio at surface (kg kg-1), used for
!                           ! diagnosing the mixing ratio at 2 m for coupled model
  real        :: snoalb         ! maximum albedo expected over deep snow
  real        :: soilrz         ! root zone soil column water content (kg m-2)
  real        :: soilrzmax      ! SY
  logical     :: rdlai2d
  logical     :: usemonalb
  real        :: ribb 
  real        :: ptu
  real        :: ustar
  integer       ::  soiltyp
  character*256 ::  llanduse, lsoil

  real        :: frzk, frzfact
  real        :: t2diag, q2diag, rho
  real        :: relsmc   ! Relative "wetness" soil moisture (fraction) (smc - smcwlt) / (porosity - smcwlt)
  integer     :: t,i,k
  real        :: wchange_prev

  REAL, PARAMETER:: R = 287.04
! Noah-3.6's value for CP is 1004.6 - D. Mocko
  REAL, PARAMETER:: CP = 1004.6
!  REAL, PARAMETER:: CH2O = 4.2E6
!  REAL, PARAMETER:: CSOIL = 2.00E+6
!  REAL, PARAMETER:: CAIR = 1004.0
!  REAL, PARAMETER:: CICE = 2.106E6
! Noah-3.6 uses this value instead of 273.16 - D. Mocko
  REAL, PARAMETER:: T0 = 273.15       ! Freezing point in Kelvin (273.15 K)
  REAL, PARAMETER:: LVH2O = 2.501000E+6 ! Latent heat for evapo for water  
  REAL, PARAMETER:: EPS = 0.622 ! Water/(dry air) molec mass ratio (epsilon)
  !--  PARAMETER USED TO CALCULATE ROUGHNESS LENGTH OF HEAT.
  !      PARAMETER(CZIL = 0.2)
  ! - Changed for version 2.6 June 2nd 2003 *
  !      PARAMETER(CZIL = 0.075)
  ! - Changed for version 2.7.1 Jesse 20041225
  !      PARAMETER(CZIL = 0.1)
  ! - CZIL now set in general parameter table
  REAL, PARAMETER :: A2=17.67,A3=273.15,A4=29.65, A23M4=A2*(A3-A4)
  REAL     :: soilhtc, soilmtc
  REAL     :: startht, startsm, startswe, startint
  REAL     :: SFCTSNO
  REAL     :: E2SAT
  REAL     :: Q2SATI
  real     :: cm, ch
  logical  :: alarmCheck, Bondvillecheck
  integer  :: iret

  real     :: julian_in
  integer  :: yr

  real, parameter :: CONST_LATSUB=2.83E+6 ! to be consistent with Noah SFLX
  integer         :: row, col

  REAL     :: WRSI_TimeStep
  REAL     :: WR_TimeStep
  REAL     :: AET_TimeStep
  REAL, parameter    :: SUMWR_EXITWRSIBELOWTHISQTY = 0
  REAL, parameter    :: SUMWR_EXITWRSIABOVETHISQTY = 10000
  REAL, parameter    :: SUMET_EXITWRSIABOVETHISQTY = 10000
  character*3        :: fnest

#if 0 
  integer :: svk_col,svk_row,ii,jj
  real    :: svk_statebf(LIS_rc%lnc(n),LIS_rc%lnr(n))
#endif
! ______________________________________________________________________

! J.Case (9/11/2014)  -- new variable declarations for computing RELSMC and 
! integrated relative soil moisture (i.e. the old "MSTAVTOT" from v2.7.1).
  real :: tempval, soiltm, soiltw, soilt
! J.Case -- end mod (9/11/2014)

#if 0 
  svk_statebf = 0.0
  
  do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     svk_col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
     svk_row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
     
     svk_statebf(svk_col,svk_row) =  svk_statebf(svk_col,svk_row) + &
          noah36_struc(n)%noah(t)%snowh
  enddo

  do jj=1,LIS_rc%lnr(n)
     do ii=1,LIS_rc%lnc(n)
        if(svk_statebf(ii,jj).gt.0) then 
           svk_statebf(ii,jj) = svk_statebf(ii,jj)/LIS_rc%nensem(n)
        endif
     enddo
  enddo

  open(100,file='state.bin',form='unformatted')
  write(100) svk_statebf
  close(100)
#endif

  write(fnest,'(i3.3)') n
  alarmCheck = LIS_isAlarmRinging(LIS_rc,"Noah36 model alarm "//trim(fnest))
  if(alarmCheck) then
     ! Get Julian day of year
     call LIS_tbotTimeUtil(julian_in,yr)
     Bondvillecheck = .false.
     do i=1,LIS_rc%nmetforc
        if (trim(LIS_rc%metforc(i)).eq."Bondville") Bondvillecheck = .true.
     enddo
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ffrozp,dt,zlvl,zlvl_wind,sfctmp,q2,sfcprs,prcp,uwind,vwind,cpcp,soldn,q2sat,lwdn,sfcspd,esat,nsoil,sldpth,ice,isurban,solnet,local,solardirect,prcprain,cosz,th2,th2v,t1v,t2v,dqsdt2,slope,shdfac,shdmin,shdmax,tbot,evp,eta,eta_kinematic,shtflx,fdown,ec,edir,et,smav,ett,esnow,drip,dew,beta,etp,gflx,flx1,flx2,flx3,snomlt,sncovr,runoff3,rc,pc,rcs,rct,rcq,rcsoil,soilw,soilm,tsoil,snoalb,soilrz,soilrzmax,rdlai2d,usemonalb,ribb,ptu,ustar,soiltyp,llanduse,lsoil,frzk,frzfact,t2diag,q2diag,rho,i,k,soilhtc,soilmtc,startht,startsm,startswe,startint,sfctsno,e2sat,q2sati,ch,cm,row,col,WRSI_TimeStep,WR_TimeStep,AET_TimeStep)

!$OMP DO
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)

        dt = LIS_rc%nts(n) ! EMK
        sfctmp = noah36_struc(n)%noah(t)%tair/noah36_struc(n)%forc_count
        q2     = noah36_struc(n)%noah(t)%qair/noah36_struc(n)%forc_count
        soldn  = noah36_struc(n)%noah(t)%swdown/noah36_struc(n)%forc_count
        lwdn   = noah36_struc(n)%noah(t)%lwdown*&
                 noah36_struc(n)%noah(t)%emiss/noah36_struc(n)%forc_count
        uwind  = (noah36_struc(n)%noah(t)%uwind)*(noah36_struc(n)%noah(t)%uwind)
        vwind  = (noah36_struc(n)%noah(t)%vwind)*(noah36_struc(n)%noah(t)%vwind)
        sfcspd = sqrt( uwind + vwind )/noah36_struc(n)%forc_count
        sfcprs = noah36_struc(n)%noah(t)%psurf/noah36_struc(n)%forc_count
        prcp   = noah36_struc(n)%noah(t)%rainf/noah36_struc(n)%forc_count
        cpcp   = noah36_struc(n)%noah(t)%rainf_c/noah36_struc(n)%forc_count

        if (noah36_struc(n)%noah(t)%q1.lt.0.0) then
           noah36_struc(n)%noah(t)%q1 = q2
        endif
        if((SFCTMP.eq.LIS_rc%udef).or.&
             (SFCTMP.eq.0).or.&
             (Q2.eq.LIS_rc%udef).or.&
             (SOLDN.eq.LIS_rc%udef).or.&
             (LWDN.eq.LIS_rc%udef).or.&
             (UWIND.eq.LIS_rc%udef).or.&
             (VWIND.eq.LIS_rc%udef).or.&
             (SFCPRS.eq.LIS_rc%udef).or.&
             (PRCP.eq.LIS_rc%udef)) then 
           write(LIS_logunit,*) 'Undefined forcing values found in Noah-3.6'
           write(LIS_logunit,*) 'for tile: ',t, &
                LIS_domain(n)%grid(LIS_domain(n)%gindex( & 
                                       LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
                                       LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row))%lat, &
                LIS_domain(n)%grid(LIS_domain(n)%gindex( & 
                                       LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
                                       LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row))%lon

           write(LIS_logunit,*) " ... Forcing values at this point: ",&
                 sfctmp,q2,soldn,lwdn,uwind,vwind,sfcspd,sfcprs,prcp,cpcp
           write(LIS_logunit,*) " Stopping program ..."
           call LIS_endrun()
        endif

        ffrozp = 0.0
        if (prcp .gt. 0.0) then
           if (sfctmp .lt. t0) then
              ffrozp = 1.0
           endif
        endif

        isurban = LIS_rc%urbanclass

        if(LIS_rc%plevel.eq.3) then 
           write(LIS_logunit,*) 'forcing ',sfctmp,q2,soldn,lwdn,uwind,vwind,&
                sfcspd,sfcprs,prcp,cpcp
        endif
        !-- Prevent Numerical Instability for Wind Speed
! For a true benchmark against the Noah-3.6 testcase from NCAR,
! comment out the line below "if(SFCSPD.le.0.01) SFCSPD=0.01".
! This line prevents numerical instability in all conditions.
! The Bondville testcase does have several hours with zero or
! near-zero wind speeds, and commenting this line out will
! reproduce the output of the testcase. - dmm
        if (.not.(Bondvillecheck)) then
        if(SFCSPD.le.0.01) SFCSPD=0.01
        endif

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
        !     ESAT = E(SFCTMP)

        !-- CALCULATE SATURATION MIXING RATIO (PABLO GRUNMANN, 05/28/98)

#if (!defined COUPLED)
        !     Q2SAT = 0.622 * ESAT /(SFCPRS - (1.-0.622)*ESAT)
        !     IF (Q2 .GE.  Q2SAT)  Q2 = Q2SAT*0.99

        ! Commented out above Q2SAT calculation and below DQSTD2 calculation
        ! using the routine and code from the Noah-3.6 released code - D. Mocko
        CALL CALHUM(SFCTMP, SFCPRS, Q2SAT, DQSDT2) ! Returns Q2SAT, DQSDT2

        ! Set the heights of the forcing variables (J. Santanello; D. Mocko)
        if ((noah36_struc(n)%zh.lt.0.0).or.(noah36_struc(n)%zm.lt.0.0)) then
           noah36_struc(n)%zh = noah36_struc(n)%noah(t)%fheight
           noah36_struc(n)%zm = noah36_struc(n)%noah(t)%fheight
           if (noah36_struc(n)%noah(t)%fheight.le.0.0) then
              write(LIS_logunit,*) 'Forcing height less than or'
              write(LIS_logunit,*) 'equal to zero!  Stopping run'
              call LIS_endrun()
           endif
        endif
        ZLVL = noah36_struc(n)%zh
        ZLVL_WIND = noah36_struc(n)%zm
#else
        Q2SAT = noah36_struc(n)%noah(t)%q2sat
        ZLVL = 0.5*noah36_struc(n)%noah(t)%z
        ZLVL_WIND = 0.5*noah36_struc(n)%noah(t)%z
#endif

        !-------------------------------------------------------------------------
        !   calculate slope of sat. specific humidity curve for penman: dqsdt2
        !-------------------------------------------------------------------------
        !     dqsdt2 = dqsdt (sfctmp, sfcprs)

        !-------------------------------------------------------------------------
        !   calc virtual temps and potential temps at grnd (sub 1) and at
        !    the 1st mdl lvl abv the grnd (sub 2). expon is cp divd by r.
        !-------------------------------------------------------------------------
        th2 = sfctmp + (0.0098 * zlvl)
        t2v = sfctmp * (1.0 + 0.61 * q2)

        t1v  =  noah36_struc(n)%noah(t)%t1 * (1.0 + 0.61 * q2)
        th2v = th2 * (1.0 + 0.61 * q2)

        local = .false. 

        !-------------------------------------------------------------------------
        ! Soil layer thicknesses (m)
        !-------------------------------------------------------------------------
        nsoil = noah36_struc(n)%nslay
        allocate(sldpth(nsoil))
        if(noah36_struc(n)%usedsoilmap.ne.0) then 
           do i=1,nsoil
              sldpth(i) = noah36_struc(n)%noah(t)%lyrthk(i)
           enddo
        else
           do i = 1,nsoil
              sldpth(i) = noah36_struc(n)%lyrthk(i)
           enddo
        endif

        !-------------------------------------------------------------------------
        !   These variables are not being used right now. 
        !-------------------------------------------------------------------------
        cosz = 0.0
        prcprain = 0.0
        solardirect = 0.0
        ! not used (soil parameters set in noah36_setsoils)
        slope = noah36_struc(n)%noah(t)%slope
        soiltyp = noah36_struc(n)%noah(t)%soiltype

        !-------------------------------------------------------------------------
        ! Maximum Albedo over very Deep Snow
        !-------------------------------------------------------------------------
        snoalb = noah36_struc(n)%noah(t)%mxsnalb
        solnet = soldn*(1.0-noah36_struc(n)%noah(t)%albedo)

        tbot = noah36_struc(n)%noah(t)%tempbot

        rdlai2d   = .true. !unused 
        usemonalb = .true. !unused
        ribb      = 0 !?? 
        !--  Phota Thermal Unit (PTU)
        PTU    = 0.0

        !==  OPTIONAL SUBROUTINE:  Calculate LW Radiation (Down)  =================
        !      CALL OBTLWDN(SFCTMP,LWDN)


        !--  Initialize SOILM for 1st timestep water balance
        !       SOILM = 0.0
        !--  Initialize ROOT ZONE COLUMN SOIL MOISTURE IN METERS (SOILRZ)
        SOILRZ = 0.0
        SOILRZMAX = 0.0 ! SY

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
        ! then let Noah-3.6 calculate its own surface aerodynamic conductance.
        if ((noah36_struc(n)%forcing_ch.eq.0).or.(noah36_struc(n)%forcing_cm.eq.0)) then
           if (noah36_struc(n)%sfcdifoption.eq.0) then
              CALL SFCDIF_OFF(ZLVL,ZLVL_WIND,noah36_struc(n)%noah(t)%z0_old,&
                   T1V,TH2V,SFCSPD,noah36_struc(n)%noah(t)%czil,         &
                   noah36_struc(n)%noah(t)%cm,noah36_struc(n)%noah(t)%ch,&
                   noah36_struc(n)%noah(t)%vegt,ISURBAN,                 &
                   noah36_struc(n)%iz0tlnd,USTAR)
           else if (noah36_struc(n)%sfcdifoption.eq.1) then
              CALL SFCDIF_MYJ(ZLVL,ZLVL_WIND,noah36_struc(n)%noah(t)%z0_old,&
                   noah36_struc(n)%noah(t)%z0brd_old,       &
                   SFCPRS,noah36_struc(n)%noah(t)%t1,              &
                   SFCTMP,noah36_struc(n)%noah(t)%Q1,              &
                   Q2,SFCSPD,noah36_struc(n)%noah(t)%czil,            &
                   RIBB,noah36_struc(n)%noah(t)%cm,              &
                   noah36_struc(n)%noah(t)%ch,              &
                   noah36_struc(n)%noah(t)%vegt,            &
                   ISURBAN,noah36_struc(n)%iz0tlnd,                 &
                   noah36_struc(n)%ztmax2,                  &
                   noah36_struc(n)%dzeta2,                  &
                   noah36_struc(n)%psih2,                   &
                   noah36_struc(n)%psim2,USTAR)
           endif
        else
           ! Still call sfcdif because we need ustar, but do not overwrite
           ! the provided ch and cm values.
           cm = noah36_struc(n)%noah(t)%cm/noah36_struc(n)%forc_count
           ch = noah36_struc(n)%noah(t)%ch/noah36_struc(n)%forc_count
           if (noah36_struc(n)%sfcdifoption.eq.0) then
              CALL SFCDIF_OFF(ZLVL,ZLVL_WIND,noah36_struc(n)%noah(t)%z0_old,&
                   T1V,TH2V,SFCSPD,noah36_struc(n)%noah(t)%czil,cm,ch,   &
                   noah36_struc(n)%noah(t)%vegt,ISURBAN,                 &
                   noah36_struc(n)%iz0tlnd,USTAR)
           else if (noah36_struc(n)%sfcdifoption.eq.1) then
              CALL SFCDIF_MYJ(ZLVL,ZLVL_WIND,noah36_struc(n)%noah(t)%z0_old,&
                   noah36_struc(n)%noah(t)%z0brd_old,       &
                   SFCPRS,noah36_struc(n)%noah(t)%t1,              &
                   SFCTMP,noah36_struc(n)%noah(t)%Q1,              &
                   Q2,SFCSPD,noah36_struc(n)%noah(t)%czil,            &
                   RIBB,cm,                                      &
                   ch,                                      &
                   noah36_struc(n)%noah(t)%vegt,            &
                   ISURBAN,noah36_struc(n)%iz0tlnd,                 &
                   noah36_struc(n)%ztmax2,                  &
                   noah36_struc(n)%dzeta2,                  &
                   noah36_struc(n)%psih2,                   &
                   noah36_struc(n)%psim2,USTAR)
           endif
        endif
#else
        ! Still call sfcdif because we need ustar, but do not overwrite
        ! the provided ch and cm values.
        cm = noah36_struc(n)%noah(t)%cm
        ch = noah36_struc(n)%noah(t)%ch
        if (noah36_struc(n)%sfcdifoption.eq.0) then
!           CALL SFCDIF_OFF(ZLVL,noah36_struc(n)%noah(t)%z0_old,&
!                T1V,TH2V,SFCSPD,noah36_struc(n)%noah(t)%czil,cm,ch,     &
!                noah36_struc(n)%noah(t)%vegt,ISURBAN,     &
!                noah36_struc(n)%iz0tlnd,USTAR)
           CALL SFCDIF_OFF(ZLVL,noah36_struc(n)%noah(t)%z0_old,&
                T1V,TH2V,SFCSPD,noah36_struc(n)%noah(t)%czil,cm,ch,     &
                USTAR)

        else if (noah36_struc(n)%sfcdifoption.eq.1) then
           CALL SFCDIF_MYJ(ZLVL,noah36_struc(n)%noah(t)%z0_old, &
                noah36_struc(n)%noah(t)%z0brd_old,          &
                SFCPRS,noah36_struc(n)%noah(t)%t1,                 &
                SFCTMP,noah36_struc(n)%noah(t)%Q1,                 &
                Q2,SFCSPD,noah36_struc(n)%noah(t)%czil,               &
                RIBB,cm,                                         &
                ch,                                         &
                noah36_struc(n)%noah(t)%vegt,               &
                ISURBAN,noah36_struc(n)%iz0tlnd,                    &
                noah36_struc(n)%ztmax2,                     &
                noah36_struc(n)%dzeta2,                     &
                noah36_struc(n)%psih2,                      &
                noah36_struc(n)%psim2,USTAR)
        endif
#endif

        noah36_struc(n)%noah(t)%kdt = noah36_struc(n)%noah(t)%refkdt* &
             noah36_struc(n)%noah(t)%dksat/&
             noah36_struc(n)%noah(t)%refdk
        frzk = noah36_struc(n)%noah(t)%frzk
        if(noah36_struc(n)%noah(t)%smcref.eq.0) then 
           noah36_struc(n)%noah(t)%frzx = 0.0
        else
           FRZFACT = (noah36_struc(n)%noah(t)%SMCMAX / noah36_struc(n)%noah(t)%SMCREF) &
                * (0.412 / 0.468)
           noah36_struc(n)%noah(t)%FRZX = FRZK * FRZFACT
        endif

! "xice" should be set in NU-WRF for a COUPLED simulation.
! For offline, "xice" is set in noah36_setvegparms.F90 - dmm
        ice = noah36_struc(n)%noah(t)%xice

#if (defined COUPLED)
        DQSDT2 = Q2SAT*A23M4/(SFCTMP-A4)**2
        
        if(noah36_struc(n)%noah(t)%sneqv.gt.0.0) then
           ! snow on surface (use ice saturation properties)
           call SNFRAC(noah36_struc(n)%noah(t)%sneqv,                  &
                noah36_struc(n)%noah(t)%snup,                          &
                noah36_struc(n)%noah(t)%salp,                          &
                noah36_struc(n)%noah(t)%snowh,SNCOVR,                  &
                noah36_struc(n)%noah(t)%lai,                           &
                noah36_struc(n)%noah(t)%shdfac,                        &
                fvb,gama,fbur,fgsn,                                    &
                noah36_struc(n)%noah(t)%ztopv,                         &
                noah36_struc(n)%noah(t)%zbotv,                         &
                noah36_struc(n)%ua_phys)
           noah36_struc(n)%noah(t)%sca = sncovr ! EMK NUWRF
           SFCTSNO = SFCTMP
           E2SAT = 611.2*EXP(6174.*(1./273.15 - 1./SFCTSNO))
           Q2SATI = 0.622*E2SAT/(SFCPRS-E2SAT)
           Q2SATI = Q2SATI/(1.0+Q2SATI)    ! spec. hum.
           if(noah36_struc(n)%noah(t)%T1.gt.273.15) then
                 ! warm ground temps, weight the saturation between ice and water according to SNOWC
              Q2SAT = Q2SAT*(1.-SNCOVR) + Q2SATI*SNCOVR
              DQSDT2 = DQSDT2*(1.-SNCOVR) + Q2SATI*6174./(SFCTSNO**2)*SNCOVR
           else
              ! cold ground temps, use ice saturation only
              Q2SAT = Q2SATI
              DQSDT2 = Q2SATI*6174./(SFCTSNO**2)
           endif
              ! for snow cover fraction at 0 C, ground temp will not change, so DQSDT2 effectively zero
           if ((noah36_struc(n)%noah(t)%T1.gt.273.).and.               &
                (SNCOVR.gt.0.)) DQSDT2 = DQSDT2*(1.-SNCOVR)
        endif
#endif

        soilmtc = 0
        do i=1,noah36_struc(n)%nslay
           soilmtc = soilmtc +                                              &
                noah36_struc(n)%noah(t)%smc(i)*LIS_CONST_RHOFW*sldpth(i)
        enddo
        ! First attempt at calculating DelSurfHeat - D. Mocko
        !   soilhtc = (((NOAH36_STRUC(N)%NOAH(T)%SH2O(1)*CH2O)                  &
        !             + ((1.0-noah36_struc(n)%noah(t)%SMCMAX)*CSOIL)            &
        !             + ((noah36_struc(n)%noah(t)%SMCMAX-NOAH36_STRUC(N)%NOAH(T)%SMC(1))*CAIR) &
        !             + ((NOAH36_STRUC(N)%NOAH(T)%SMC(1)-NOAH36_STRUC(N)%NOAH(T)%SH2O(1))*CICE)) &
        !             * noah36_struc(n)%noah(t)%stc(1) * SLDPTH(1))
        soilhtc = 0.0
        startht = soilhtc
        startsm = soilmtc
        startswe = noah36_struc(n)%noah(t)%sneqv*LIS_CONST_RHOFW
        startint = noah36_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW

#define _DEBUG_PRINT_ 0
#if _DEBUG_PRINT_
        print*, "Before SFLX."
        print*, 'b  FFROZP = ', FFROZP
        print*, 'b  ICE = ', ICE
        print*, 'b  ISURBAN = ', ISURBAN
        print*, 'b  DT = ', DT
        print*, 'b  ZLVL = ', ZLVL
        print*, 'b  ZLVL_WIND = ', ZLVL_WIND
        print*, 'b  NSOIL = ', NSOIL
        print*, 'b  SLDPTH = ', SLDPTH
#if 0
        print*, 'b  LLANDUSE = ', trim(LLANDUSE)
        print*, 'b  LSOIL = ', trim(LSOIL)
#endif
        print*, 'b  LWDN = ', LWDN
        print*, 'b  SOLDN = ', SOLDN
        print*, 'b  SOLNET = ', SOLNET
        print*, 'b  SFCPRS = ', SFCPRS
        print*, 'b  PRCP = ', PRCP
        print*, 'b  SFCTMP = ', SFCTMP
        print*, 'b  Q2 = ', Q2
        print*, 'b  SFCSPD = ', SFCSPD
#if 0
        print*, 'b  COSZ = ', COSZ
        print*, 'b  PRCPRAIN = ', PRCPRAIN
        print*, 'b  SOLARDIRECT = ', SOLARDIRECT
#endif
        print*, 'b  TH2 = ', TH2
        print*, 'b  Q2SAT = ', Q2SAT
        print*, 'b  DQSDT2 = ', DQSDT2
        print*, 'b  VEGTYP = ', noah36_struc(n)%noah(t)%vegt
        print*, 'b  SOILTYP = ', SOILTYP
        print*, 'b  SLOPETYP = ',1
        print*, 'b  SHDFAC = ', noah36_struc(n)%noah(t)%shdfac
        print*, 'b  SHDMIN = ', noah36_struc(n)%noah(t)%shdmin
        print*, 'b  SHDMAX = ', noah36_struc(n)%noah(t)%shdmax
        print*, 'b  ALB = ', noah36_struc(n)%noah(t)%alb
        print*, 'b  SNOALB = ', SNOALB
        print*, 'b  TBOT = ', TBOT
        print*, 'b  Z0BRD = ', noah36_struc(n)%noah(t)%z0brd
        print*, 'b  Z0 = ', noah36_struc(n)%noah(t)%z0
        print*, 'b  EMISSI = ', noah36_struc(n)%noah(t)%emiss
        print*, 'b  EMBRD = ', noah36_struc(n)%noah(t)%embrd
        print*, 'b  CMC = ', noah36_struc(n)%noah(t)%cmc
        print*, 'b  T1 = ', noah36_struc(n)%noah(t)%t1
        print*, 'b  STC = ', noah36_struc(n)%noah(t)%stc
        print*, 'b  SMC = ', noah36_struc(n)%noah(t)%smc
        print*, 'b  SH2O = ', noah36_struc(n)%noah(t)%sh2o
        print*, 'b  SNOWH = ', noah36_struc(n)%noah(t)%snowh
        print*, 'b  SNEQV = ', noah36_struc(n)%noah(t)%sneqv
        print*, 'b  ALBEDO = ', noah36_struc(n)%noah(t)%albedo
        print*, 'b  CH = ', noah36_struc(n)%noah(t)%ch
        print*, 'b  CM = ', noah36_struc(n)%noah(t)%cm
        print*, 'b  ETA = ', ETA
        print*, 'b  SHEAT = ', shtflx
        print*, 'b  ETAKIN = ', eta_kinematic
        print*, 'b  FDOWN = ', FDOWN
        print*, 'b  EC = ', EC
        print*, 'b  EDIR = ', EDIR
        print*, 'b  ET = ', ET(1), ET(2), ET(3), ET(4)
        print*, 'b  ETT = ', ETT
        print*, 'b  ESNOW = ', ESNOW
        print*, 'b  DRIP = ', DRIP
        print*, 'b  DEW = ', DEW
        print*, 'b  BETA = ', BETA
        print*, 'b  ETP = ', ETP
        print*, 'b  SSOIL = ', gflx
        print*, 'b  FLX1 = ', FLX1
        print*, 'b  FLX2 = ', FLX2
        print*, 'b  FLX3 = ', FLX3
        print*, 'b  SNOMLT = ', SNOMLT
        print*, 'b  SNCOVR = ', SNCOVR
        print*, 'b  RUNOFF1 = ', noah36_struc(n)%noah(t)%RUNOFF1
        print*, 'b  RUNOFF2 = ', noah36_struc(n)%noah(t)%RUNOFF2
        print*, 'b  RUNOFF3 = ', RUNOFF3
        print*, 'b  RC = ', RC
        print*, 'b  PC = ', PC
        print*, 'b  RSMIN = ', noah36_struc(n)%noah(t)%rsmin
        print*, 'b  XLAI = ', noah36_struc(n)%noah(t)%lai
        print*, 'b  RCS = ', RCS
        print*, 'b  RCT = ', RCT
        print*, 'b  RCQ = ', RCQ
        print*, 'b  RCSOIL = ', RCSOIL
        print*, 'b  SOILW = ', SOILW
        print*, 'b  SOILM = ', SOILM
        print*, 'b  Q1 = ', noah36_struc(n)%noah(t)%Q1
#if 0
        print*, 'b  RDLAI2D = ', RDLAI2D
        print*, 'b  USEMONALB = ', USEMONALB
#endif
        print*, 'b  SNOTIME1 = ', noah36_struc(n)%noah(t)%snotime1
        print*, 'b  RIBB = ', RIBB
        print*, 'b  SMCWLT = ', noah36_struc(n)%noah(t)%smcwlt
        print*, 'b  SMCDRY = ', noah36_struc(n)%noah(t)%smcdry
        print*, 'b  SMCREF = ', noah36_struc(n)%noah(t)%smcref
        print*, 'b  SMCMAX = ', noah36_struc(n)%noah(t)%smcmax
        print*, 'b  NROOT = ', noah36_struc(n)%noah(t)%nroot
        print *,' '
        print *,'CALLING SFLX'
        print *,' '
#endif

        if (ice.eq.0) then 
           col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
           row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
           call sflx(                                                                 &
                ffrozp, isurban, dt, zlvl, nsoil, sldpth,                             &
                local,                                                                &
                llanduse, lsoil,                                                      & 
                lwdn, soldn, solnet, sfcprs, prcp, sfctmp, q2, sfcspd,                &
                cosz, prcprain, solardirect,                                          & 
                th2, q2sat, dqsdt2,                                                   & 
                noah36_struc(n)%noah(t)%vegt, soiltyp, slope, & 
                noah36_struc(n)%noah(t)%shdfac, noah36_struc(n)%noah(t)%shdmin,       & 
                noah36_struc(n)%noah(t)%shdmax, & 
                noah36_struc(n)%noah(t)%alb, snoalb, tbot,noah36_struc(n)%noah(t)%z0brd, & 
                noah36_struc(n)%noah(t)%z0, noah36_struc(n)%noah(t)%emiss,            & 
                noah36_struc(n)%noah(t)%embrd,    &
                noah36_struc(n)%noah(t)%cmc, noah36_struc(n)%noah(t)%t1,              & 
                noah36_struc(n)%noah(t)%stc, noah36_struc(n)%noah(t)%smc,             & 
                noah36_struc(n)%noah(t)%sh2o, noah36_struc(n)%noah(t)%snowh,          & 
                noah36_struc(n)%noah(t)%sneqv, noah36_struc(n)%noah(t)%albedo,        & 
                noah36_struc(n)%noah(t)%ch, noah36_struc(n)%noah(t)%cm,               & 
                evp,eta, shtflx, eta_kinematic, fdown,                                & 
                ec, edir, et, ett, esnow, drip, dew,                                  & 
                beta, etp, gflx, flx1, flx2, flx3,                                    & 
                flx4, fvb, fbur, fgsn, noah36_struc(n)%ua_phys,                       &
                noah36_struc(n)%noah(t)%ztopv, noah36_struc(n)%noah(t)%zbotv,         &
                snomlt, sncovr, noah36_struc(n)%noah(t)%runoff1,                      &
                noah36_struc(n)%noah(t)%runoff2, runoff3,                             & 
                rc, pc, noah36_struc(n)%noah(t)%rsmin,                                & 
                noah36_struc(n)%noah(t)%lai, rcs, rct, rcq, rcsoil,                   &      
                soilw, soilm, noah36_struc(n)%noah(t)%q1, smav, rdlai2d,              &
                usemonalb, noah36_struc(n)%noah(t)%snotime1, ribb,                    & 
                noah36_struc(n)%noah(t)%smcwlt, noah36_struc(n)%noah(t)%smcdry,       & 
                noah36_struc(n)%noah(t)%smcref, noah36_struc(n)%noah(t)%smcmax,       &
                noah36_struc(n)%noah(t)%nroot,                                   & 
                noah36_struc(n)%noah(t)%hs, noah36_struc(n)%noah(t)%rgl,         & 
                noah36_struc(n)%noah(t)%snup,                                    & 
                noah36_struc(n)%noah(t)%cmcmax, noah36_struc(n)%noah(t)%rsmax,   & 
                noah36_struc(n)%noah(t)%topt, noah36_struc(n)%noah(t)%bexp,      & 
                noah36_struc(n)%noah(t)%dksat, noah36_struc(n)%noah(t)%dwsat,    & 
                noah36_struc(n)%noah(t)%f1, noah36_struc(n)%noah(t)%quartz,      & 
                noah36_struc(n)%noah(t)%psisat, noah36_struc(n)%noah(t)%czil,    & 
                noah36_struc(n)%noah(t)%sbeta, noah36_struc(n)%noah(t)%fxexp,    & 
                noah36_struc(n)%noah(t)%csoil, noah36_struc(n)%noah(t)%salp,     & 
                noah36_struc(n)%noah(t)%kdt, noah36_struc(n)%noah(t)%cfactr, & 
                noah36_struc(n)%noah(t)%zbot, noah36_struc(n)%noah(t)%refkdt, ptu, & 
                noah36_struc(n)%noah(t)%frzx,                                      &
                noah36_struc(n)%noah(t)%sndens, & !added for use in SCF DA, yliu 
                noah36_struc(n)%noah(t)%lvcoef, tsoil,                           &
                sfhead1rt, infxs1rt, etpnd1, noah36_struc(n)%snowfix,&
                !ag (05Jan2021)
                noah36_struc(n)%noah(t)%rivsto,& ! in   - river storage [m/s] 
                noah36_struc(n)%noah(t)%fldsto,& ! in   - flood storage [m/s]
                noah36_struc(n)%noah(t)%fldfrc)  ! in   - flooded fraction [-]

           noah36_struc(n)%noah(t)%sca = sncovr ! EMK NUWRF
           noah36_struc(n)%noah(t)%z0_old = noah36_struc(n)%noah(t)%z0

#if _DEBUG_PRINT_
        print *,' '
        print*, "*****--------------After SFLX.------------****"
        print *,' '
        print*, 'a  FFROZP = ', FFROZP
        print*, 'a  ICE = ', ICE
        print*, 'a  ISURBAN = ', ISURBAN
        print*, 'a  DT = ', DT
        print*, 'a  ZLVL = ', ZLVL
        print*, 'a  ZLVL_WIND = ', ZLVL_WIND
        print*, 'a  NSOIL = ', NSOIL
        print*, 'a  SLDPTH = ', SLDPTH
#if 0
        print*, 'a  LLANDUSE = ', trim(LLANDUSE)
        print*, 'a  LSOIL = ', trim(LSOIL)
#endif
        print*, 'a  LWDN = ', LWDN
        print*, 'a  SOLDN = ', SOLDN
        print*, 'a  SOLNET = ', SOLNET
        print*, 'a  SFCPRS = ', SFCPRS
        print*, 'a  PRCP = ', PRCP
        print*, 'a  SFCTMP = ', SFCTMP
        print*, 'a  Q2 = ', Q2
        print*, 'a  SFCSPD = ', SFCSPD
#if 0
        print*, 'a  COSZ = ', COSZ
        print*, 'a  PRCPRAIN = ', PRCPRAIN
        print*, 'a  SOLARDIRECT = ', SOLARDIRECT
#endif
        print*, 'a  TH2 = ', TH2
        print*, 'a  Q2SAT = ', Q2SAT
        print*, 'a  DQSDT2 = ', DQSDT2
        print*, 'a  VEGTYP = ', noah36_struc(n)%noah(t)%vegt
        print*, 'a  SOILTYP = ', SOILTYP
        print*, 'a  SLOPETYP = ',1
        print*, 'a  SHDFAC = ', noah36_struc(n)%noah(t)%shdfac
        print*, 'a  SHDMIN = ', noah36_struc(n)%noah(t)%shdmin
        print*, 'a  SHDMAX = ', noah36_struc(n)%noah(t)%shdmax
        print*, 'a  ALB = ', noah36_struc(n)%noah(t)%alb
        print*, 'a  SNOALB = ', SNOALB
        print*, 'a  TBOT = ', TBOT
        print*, 'a  Z0BRD = ', noah36_struc(n)%noah(t)%z0brd
        print*, 'a  Z0 = ', noah36_struc(n)%noah(t)%z0
        print*, 'a  EMISSI = ', noah36_struc(n)%noah(t)%emiss
        print*, 'a  EMBRD = ', noah36_struc(n)%noah(t)%embrd
        print*, 'a  CMC = ', noah36_struc(n)%noah(t)%cmc
        print*, 'a  T1 = ', noah36_struc(n)%noah(t)%t1
        print*, 'a  STC = ', noah36_struc(n)%noah(t)%stc
        print*, 'a  SMC = ', noah36_struc(n)%noah(t)%smc
        print*, 'a  SH2O = ', noah36_struc(n)%noah(t)%sh2o
        print*, 'a  SNOWH = ', noah36_struc(n)%noah(t)%snowh
        print*, 'a  SNEQV = ', noah36_struc(n)%noah(t)%sneqv
        print*, 'a  ALBEDO = ', noah36_struc(n)%noah(t)%albedo
        print*, 'a  CH = ', noah36_struc(n)%noah(t)%ch
        print*, 'a  CM = ', noah36_struc(n)%noah(t)%cm
        print*, 'a  ETA = ', ETA
        print*, 'a  SHEAT = ', shtflx
        print*, 'a  ETAKIN = ', eta_kinematic
        print*, 'a  FDOWN = ', FDOWN
        print*, 'a  EC = ', EC
        print*, 'a  EDIR = ', EDIR
        print*, 'a  ET = ', ET(1), ET(2), ET(3), ET(4)
        print*, 'a  ETT = ', ETT
        print*, 'a  ESNOW = ', ESNOW
        print*, 'a  DRIP = ', DRIP
        print*, 'a  DEW = ', DEW
        print*, 'a  BETA = ', BETA
        print*, 'a  ETP = ', ETP
        print*, 'a  SSOIL = ', gflx
        print*, 'a  FLX1 = ', FLX1
        print*, 'a  FLX2 = ', FLX2
        print*, 'a  FLX3 = ', FLX3
        print*, 'a  SNOMLT = ', SNOMLT
        print*, 'a  SNCOVR = ', SNCOVR
        print*, 'a  RUNOFF1 = ', noah36_struc(n)%noah(t)%RUNOFF1
        print*, 'a  RUNOFF2 = ', noah36_struc(n)%noah(t)%RUNOFF2
        print*, 'a  RUNOFF3 = ', RUNOFF3
        print*, 'a  RC = ', RC
        print*, 'a  PC = ', PC
        print*, 'a  RSMIN = ', noah36_struc(n)%noah(t)%rsmin
        print*, 'a  XLAI = ', noah36_struc(n)%noah(t)%lai
        print*, 'a  RCS = ', RCS
        print*, 'a  RCT = ', RCT
        print*, 'a  RCQ = ', RCQ
        print*, 'a  RCSOIL = ', RCSOIL
        print*, 'a  SOILW = ', SOILW
        print*, 'a  SOILM = ', SOILM
        print*, 'a  Q1 = ', noah36_struc(n)%noah(t)%Q1
#if 0
        print*, 'a  RDLAI2D = ', RDLAI2D
        print*, 'a  USEMONALB = ', USEMONALB
#endif
        print*, 'a  SNOTIME1 = ', noah36_struc(n)%noah(t)%snotime1
        print*, 'a  RIBB = ', RIBB
        print*, 'a  SMCWLT = ', noah36_struc(n)%noah(t)%smcwlt
        print*, 'a  SMCDRY = ', noah36_struc(n)%noah(t)%smcdry
        print*, 'a  SMCREF = ', noah36_struc(n)%noah(t)%smcref
        print*, 'a  SMCMAX = ', noah36_struc(n)%noah(t)%smcmax
        print*, 'a  NROOT = ', noah36_struc(n)%noah(t)%nroot
#endif

        elseif (ice.eq.1) then 
           noah36_struc(n)%noah(t)%stc = 273.16
           noah36_struc(n)%noah(t)%smc = 1.0
           noah36_struc(n)%noah(t)%sh2o = 1.0
           noah36_struc(n)%noah(t)%lai = 0.01
           TBOT=271.16
        elseif (ice.eq.-1) then 
       !
       ! Set values that the LSM is expected to update,
       ! but don't get updated for glacial points.
       !
           SOILM = 0.0 !BSINGH(PNNL)- SOILM is undefined for this case, it is used for diagnostics so setting it to zero
           noah36_struc(n)%noah(t)%lai = 0.01 ! KWM Should this be Zero over land ice?  Does this value matter?
           noah36_struc(n)%noah(t)%RUNOFF2 = 0.0
           RUNOFF3 = 0.0
        !
        ! Soil moisture fields set to 1.0 for glacial points.
        !
           noah36_struc(n)%noah(t)%smc = 1.0
           noah36_struc(n)%noah(t)%sh2o = 1.0
           SMAV = 0.0
        !
        ! EDIR, ETT, EC, EVP, and RC need to be set to zero to be sane.
        !
           EDIR = 0.0
           ETT  = 0.0
           EC   = 0.0
           EVP  = 0.0 
           RC   = 0.0
        !
        ! SHDFAC, set elsewhere in the driver, is overwritten here.
        ! Not that SFLX_GLACIAL uses it, but for consistency, in that
        ! we assume there is no vegetation on a glacial point.
        !
           noah36_struc(n)%noah(t)%shdfac = 0.0
        !
        ! For glacial, set ALB to the ALBEDOMAX (from VEGPARM.TBL) for
        ! snow/ice points.  Similarly, set Z0BRT to the Z0MIN for
        ! snow/ice points.
        !
! Set albedo/z0brd for "ice=-1" in noah36_setvegparms.F90 - dmm
!        noah36_struc(n)%noah(t)%alb = ALBEDOMAXTBL(VEGTYP)
!        noah36_struc(n)%noah(t)%z0brd = Z0MINTBL(VEGTYP)

        !
        ! Call the Noah LSM routines for Glacial Ice points.
        !
           ! SVK/EMK...Make sure sncovr is initialized!
           sncovr = noah36_struc(n)%noah(t)%sca 

           CALL SFLX_GLACIAL(ICE,FFROZP,DT,ZLVL,NSOIL,SLDPTH,          &
                LWDN,SOLNET,SFCPRS,PRCP,SFCTMP,Q2,                     &
                TH2,Q2SAT,DQSDT2,                                      &
                noah36_struc(n)%noah(t)%alb,SNOALB,TBOT,               &
                noah36_struc(n)%noah(t)%Z0BRD,                         &
                noah36_struc(n)%noah(t)%Z0,                            &
                noah36_struc(n)%noah(t)%EMISS,                         &
                noah36_struc(n)%noah(t)%EMBRD,                         &
                noah36_struc(n)%noah(t)%T1,                            &
                noah36_struc(n)%noah(t)%STC,                           &
                noah36_struc(n)%noah(t)%SNOWH,                         &
                noah36_struc(n)%noah(t)%SNEQV,                         &
                noah36_struc(n)%noah(t)%ALBEDO,                        &
                noah36_struc(n)%noah(t)%CH,                            &
                ETA,SHTFLX,ETA_KINEMATIC,FDOWN,                        &
                ESNOW,DEW,                                             &
                ETP,GFLX,                                              &
                FLX1,FLX2,FLX3,                                        &
                SNOMLT,SNCOVR,                                         &
                noah36_struc(n)%noah(t)%RUNOFF1,                       &
                noah36_struc(n)%noah(t)%Q1,                            &
                noah36_struc(n)%noah(t)%SNOTIME1,                      &
                RIBB,noah36_struc(n)%noah(t)%lvcoef,tsoil)
        endif

! Save variables for passing to WRF or for output
        noah36_struc(n)%noah(t)%qh = shtflx
        noah36_struc(n)%noah(t)%qg = gflx
        noah36_struc(n)%noah(t)%qle = eta
        noah36_struc(n)%noah(t)%eta_kinematic = eta_kinematic

! J.Case (9/11/2014) -- This block of code should be done only in coupled runs.
#if (defined COUPLED)
        noah36_struc(n)%noah(t)%chs2 = noah36_struc(n)%noah(t)%cqs2
! J.Case (9/11/2014) -- Moved 2 lines below to ensure qsfc is defined before the if-block after the line
        noah36_struc(n)%noah(t)%qsfc = noah36_struc(n)%noah(t)%q1/&
             (1-noah36_struc(n)%noah(t)%q1)
        if(noah36_struc(n)%noah(t)%q1.gt.noah36_struc(n)%noah(t)%qsfc) then 
           noah36_struc(n)%noah(t)%cqs2 = noah36_struc(n)%noah(t)%chs2
        endif
#endif
! J.Case (9/11/2014, end mods)

        if (SFCTMP .lt. noah36_struc(n)%noah(t)%tair_agl_min) then
           noah36_struc(n)%noah(t)%tair_agl_min = SFCTMP
           noah36_struc(n)%noah(t)%rhmin = Q2/Q2SAT
        endif

        noah36_struc(n)%noah(t)%soilm = soilm*1000.0
        noah36_struc(n)%noah(t)%qs    = noah36_struc(n)%noah(t)%runoff1 *  &
                                        LIS_CONST_RHOFW
        noah36_struc(n)%noah(t)%qsb   = noah36_struc(n)%noah(t)%runoff2 *  &
                                        LIS_CONST_RHOFW
        noah36_struc(n)%noah(t)%qsm   = snomlt*LIS_CONST_RHOFW/dt

! Save variables to output
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_SWNET,value=soldn * &
             (1.0-noah36_struc(n)%noah(t)%albedo),vlevel=1,unit="W m-2",&
             direction="DN",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_LWNET,vlevel=1,  &
             value=(-1.0*(noah36_struc(n)%noah(t)%emiss*                &
             ((5.67E-8)*(NOAH36_STRUC(N)%NOAH(T)%T1**4.0))-             &
             (noah36_struc(n)%noah(t)%lwdown*                           &
             noah36_struc(n)%noah(t)%emiss))),unit="W m-2",&
             direction="DN",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QLE,value=eta,            &
             vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QH,value=shtflx,          &
             vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QG,value=-gflx,           &
             vlevel=1,unit="W m-2",direction="DN",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QF,value=flx3,            &
             vlevel=1,unit="W m-2",direction="S2L",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QV,value=esnow,           &
             vlevel=1,unit="W m-2",direction="S2V",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QTAU,                     &
             value=(SFCPRS/(R*T2V)*noah36_struc(n)%noah(t)%CM*USTAR),     &
             vlevel=1,unit="N m-2",direction="DN",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QA,value=-flx1,           &
             vlevel=1,unit="W m-2",direction="DN",surface_type=LIS_rc%lsm_index)
        ! First attempt at calculating DelSurfHeat - D. Mocko
        !     soilhtc = (((NOAH36_STRUC(N)%NOAH(T)%SH2O(1)*CH2O)                &
        !               + ((1.0-noah36_struc(n)%noah(t)%SMCMAX)*CSOIL)          &
        !               + ((noah36_struc(n)%noah(t)%SMCMAX-NOAH36_STRUC(N)%NOAH(T)%SMC(1))*CAIR) &
        !               + ((NOAH36_STRUC(N)%NOAH(T)%SMC(1)-NOAH36_STRUC(N)%NOAH(T)%SH2O(1))*CICE)) &
        !               * noah36_struc(n)%noah(t)%stc(1) * SLDPTH(1))
        soilhtc = 0.0
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_DELSURFHEAT,              &
             value=(startht-soilhtc),vlevel=1,unit="J m-2",&
             direction="INC",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_DELCOLDCONT,              &
             value=((flx2 - esnow) * DT),vlevel=1,unit="J m-2",&
             direction="INC",surface_type=LIS_rc%lsm_index)
        !Bowen Ratio - sensible/latent
        if (eta.gt.0) then 
           call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_BR,value=(shtflx/eta), &
                vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        else
           call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_BR,value=0.0,          &
                vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        endif
        !Evaporative Fraction
        if ((eta+shtflx).ne.0) then 
           call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_EF,                    &
                value=abs(eta/(eta+shtflx)),vlevel=1,unit="-",&
                direction="-",surface_type=LIS_rc%lsm_index)
        else
           !double check
           call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_EF,value=1.0,          &
                vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        endif

        ! Noah-3.6 uses this value instead of 273.16 - D. Mocko
        if (sfctmp .lt. T0) then
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=prcp,       &
                vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=0.0,        &
                vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=prcp*dt,    &
                vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=0.0,        &
                vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
        else
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=0.0,        &
                vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=prcp,       &
                vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWF,value=0.0,        &
                vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINF,value=prcp*dt,    &
                vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)
        endif

!        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TOTALPRECIP,               &
!             value=prcp,vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
!        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TOTALPRECIP,               &
!             value=prcp*dt,vlevel=1,unit="kg m-2",direction="DN",surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EVAP,value=evp,            &
             vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EVAP, value=(evp*3600.),   &
             vlevel=1,unit="mm hr-1",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QS,                        &
             value= noah36_struc(n)%noah(t)%runoff1*LIS_CONST_RHOFW,      &
             vlevel=1,unit="kg m-2 s-1",direction="OUT",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QS,                        &
             value= noah36_struc(n)%noah(t)%runoff1*LIS_CONST_RHOFW*dt,   &
             vlevel=1,unit="kg m-2",direction="OUT",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSB,                       &
             value= noah36_struc(n)%noah(t)%runoff2*LIS_CONST_RHOFW,      &
             vlevel=1,unit="kg m-2 s-1",direction="OUT",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSB,                       &
             value= noah36_struc(n)%noah(t)%runoff2*LIS_CONST_RHOFW*dt,   &
             vlevel=1,unit="kg m-2",direction="OUT",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSM,                       &
             value= snomlt*LIS_CONST_RHOFW/dt,vlevel=1,unit="kg m-2 s-1",&
             direction="S2L",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_QSM,                       &
             value= snomlt*LIS_CONST_RHOFW,vlevel=1,unit="kg m-2",&
             direction="S2L",surface_type=LIS_rc%lsm_index)
        soilmtc = 0
        do i=1,noah36_struc(n)%nslay
           soilmtc = soilmtc + noah36_struc(n)%noah(t)%smc(i) *           &
                LIS_CONST_RHOFW * sldpth(i)
        enddo
        wchange_prev = prcp - evp - (noah36_struc(n)%noah(t)%runoff1+&
             noah36_struc(n)%noah(t)%runoff2)*LIS_CONST_RHOFW

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSURFSTOR,              &
             value=(noah36_struc(n)%noah(t)%wchange_prev - wchange_prev),vlevel=1,unit="kg m-2 s-1",&
             direction="INC",surface_type=LIS_rc%lsm_index)
        noah36_struc(n)%noah(t)%wchange_prev = wchange_prev

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSOILMOIST,              &
             value=(soilmtc-startsm),vlevel=1,unit="kg m-2",&
             direction="INC",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSWE,                    &
             value=(noah36_struc(n)%noah(t)%sneqv*LIS_CONST_RHOFW-        &
             startswe),vlevel=1,unit="kg m-2",direction="INC",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELINTERCEPT,              &
             value=(noah36_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW-startint),&
             vlevel=1,unit="kg m-2",direction="INC",surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_BARESOILT,value=TSOIL,     &
             vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_AVGSURFT,value=(TSOIL *    &
             (1.0 - SNCOVR)) + (NOAH36_STRUC(N)%NOAH(T)%T1 * SNCOVR),     &
             vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RADT,                      &
             value=NOAH36_STRUC(N)%NOAH(T)%T1,vlevel=1,unit="K",&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ALBEDO,                    &
             value=noah36_struc(n)%noah(t)%albedo,vlevel=1,unit="-",&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ALBEDO,                    &
             value=(noah36_struc(n)%noah(t)%albedo*100.),vlevel=1,&
             unit="%",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWE,vlevel=1,unit="kg m-2", &
             value=noah36_struc(n)%noah(t)%sneqv*LIS_CONST_RHOFW,&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWE,vlevel=1,unit="m",     &
             value=noah36_struc(n)%noah(t)%sneqv,&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWDEPTH,                 &
             vlevel=1,value=noah36_struc(n)%noah(t)%snowh,unit="m",&
             direction="-",surface_type=LIS_rc%lsm_index)
        noah36_struc(n)%noah(t)%sca = sncovr
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWCOVER,value=SNCOVR,    &
             vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        ! EMK...Support percentage
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNOWCOVER,value=SNCOVR*100.0,    &
             vlevel=1,unit="%",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RHMIN,                     &
             value=noah36_struc(n)%noah(t)%rhmin,vlevel=1,unit="-",direction="-",&
             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RHMIN,                     &
             value=(noah36_struc(n)%noah(t)%rhmin*100.),vlevel=1,unit="%",&
             direction="-",surface_type=LIS_rc%lsm_index)

        do i=1,noah36_struc(n)%nslay
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILMOIST,vlevel=i,     &
                value=noah36_struc(n)%noah(t)%smc(i)*sldpth(i)*           &
                LIS_CONST_RHOFW,unit='kg m-2',direction="-",surface_type=LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILMOIST,vlevel=i,     &
                value=noah36_struc(n)%noah(t)%smc(i),unit='m^3 m-3',&
                direction="-",surface_type=LIS_rc%lsm_index)
        enddo
        do i=1,noah36_struc(n)%nslay
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILTEMP,vlevel=i,      &
                value=noah36_struc(n)%noah(t)%stc(i),unit="K",&
                direction="-",surface_type=LIS_rc%lsm_index)
        enddo
        do i=1,noah36_struc(n)%nslay
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SMLIQFRAC,vlevel=i,     &
                value=(noah36_struc(n)%noah(t)%sh2o(i)/                   &
                (noah36_struc(n)%noah(t)%sh2o(i) +                 &
                noah36_struc(n)%noah(t)%smc(i))),unit='-',&
                direction="-",surface_type=LIS_rc%lsm_index)
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SMLIQFRAC,vlevel=i,     &
                value=noah36_struc(n)%noah(t)%sh2o(i),unit='m^3 m-3',&
                direction="-",surface_type=LIS_rc%lsm_index)
        enddo

        do i=1,noah36_struc(n)%nslay
           relsmc = &
                (noah36_struc(n)%noah(t)%smc(i) - noah36_struc(n)%noah(t)%smcwlt) / &
                (noah36_struc(n)%noah(t)%smcmax - noah36_struc(n)%noah(t)%smcwlt)

           if ( relsmc > 1.0 ) then
              relsmc  = 1.0
           endif
           if ( relsmc  < 0.01 ) then
              relsmc  = 0.01
           endif

           noah36_struc(n)%noah(t)%relsmc(i) = relsmc

! J.Case (9/11/2014) -- Set relative soil moisture to missing (LIS_rc%udef)
! if the vegetation type is urban class.
           if (noah36_struc(n)%noah(t)%vegt .eq. isurban) then
             relsmc = LIS_rc%udef
           endif
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RELSMC,vlevel=i, &
                value=relsmc, unit='-',direction="-",surface_type=&
                LIS_rc%lsm_index)
           if ( relsmc == LIS_rc%udef ) then
              tempval = relsmc
           else
              tempval = relsmc*100.0
           endif
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RELSMC,vlevel=i, &
                value=tempval, unit='%',direction="-",surface_type=&
                LIS_rc%lsm_index)
        enddo

        do i=1,noah36_struc(n)%nslay
           call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILET,vlevel=i,        &
                value=et(i),unit="W m-2",direction="-",surface_type=LIS_rc%lsm_index)
        enddo


! J.Case -- borrowed from LIS6/module_sf_noah271lsm.F90.
! ----------------------------------------------------------------------
! TOTAL COL SOIL MOISTURE AVAIL RELATIVE TO POROSITY/SATURATION (SOILT) 
!  (aka, MSTAVTOT)
! ----------------------------------------------------------------------
     soiltm = 0.0
     soiltw = 0.0
     soilt  = 0.0
     soiltm = &
       (noah36_struc(n)%noah(t)%smcmax - noah36_struc(n)%noah(t)%smcwlt) * &
        noah36_struc(n)%lyrthk(1)
     soiltw = &
       (noah36_struc(n)%noah(t)%smc(1) - noah36_struc(n)%noah(t)%smcwlt) * &
        noah36_struc(n)%lyrthk(1)
     do i=2,noah36_struc(n)%nslay
       soiltm = soiltm + &
         (noah36_struc(n)%noah(t)%smcmax - noah36_struc(n)%noah(t)%smcwlt) * &
         (noah36_struc(n)%lyrthk(i))
       soiltw = soiltw + &
         (noah36_struc(n)%noah(t)%smc(i) - noah36_struc(n)%noah(t)%smcwlt) * &
         (noah36_struc(n)%lyrthk(i))
     enddo
! J.Case (4/25/2013) -- Modified to correct values out of bounds (mainly for urban points).
     if (soiltm .gt. 0) then
       soilt = soiltw / soiltm
     else
       soilt = 0.01
     endif
     if (soilt > 1.0)  soilt = 1.0
     if (soilt < 0.01) soilt = 0.01

! J.Case (5/1/2013) -- Set column-integrated relative soil moisture to missing
! if the vegetation type is urban class.
        if (noah36_struc(n)%noah(t)%vegt .eq. isurban) then
          soilt = LIS_rc%udef
        endif

! J.Case NOTE: The variable below is not the same as integrated relative soil moisture 
! (i.e. MSTAVTOT) from noah2.7.1, even though it is being stored in the same variable 
! "LIS_MOC_SOILWET".  The unit label of "-" is not correct, because "soilm" is actually 
! the depth of the integrated soil moisture in millimeters.
! So, I added the integrated relative soil moisture (soilt) in the first call with 
! units of "%", while the 2nd call is the same as the original, except changed units 
! to "mm" (and multiplying by 1000), as it is defined in line #986 in LIS_histDataMod.F90.
!
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILWET, &
             value=noah36_struc(n)%noah(t)%soilm, &
             vlevel=1,unit="mm",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILWET,value=soilt,       &
             vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        if ( soilt == LIS_rc%udef ) then
           tempval = soilt
        else
           tempval = soilt*100.0
        endif
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SOILWET,value=tempval,  &
             vlevel=1,unit="%",direction="-",surface_type=LIS_rc%lsm_index)
! J.Case (9/11/2014) -- end mods from LIS6/Noah2.7.1

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP,                   &
             value=(etp/(((1.-SNCOVR)*LIS_CONST_LATVAP) +                 &
             SNCOVR*CONST_LATSUB)),vlevel=1,unit="kg m-2 s-1",                &
             direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP,                   &
             value=((etp/(((1.-SNCOVR)*LIS_CONST_LATVAP) +                &
             SNCOVR*CONST_LATSUB))*3600.),vlevel=1,unit="mm hr-1 ",         &
             direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_POTEVAP,                   &
             value=etp,vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ECANOP,                    &
             value=(ec/LIS_CONST_LATVAP),vlevel=1,unit="kg m-2 s-1",&
             direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ECANOP,                    &
             value=((ec/LIS_CONST_LATVAP)*3600.),vlevel=1,unit="mm hr-1",&
             direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ECANOP,                    &
             value=ec,vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TVEG,                      &
             value=(ett/LIS_CONST_LATVAP),vlevel=1,unit="kg m-2 s-1",&
             direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TVEG,                      &
             value=((ett/LIS_CONST_LATVAP)*3600.),vlevel=1,unit="mm hr-1",&
             direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TVEG,                      &
             value=ett,vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ESOIL,                     &
             value=(edir/LIS_CONST_LATVAP),vlevel=1,unit="kg m-2 s-1",&
             direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ESOIL,                     &
             value=((edir/LIS_CONST_LATVAP)*3600.),vlevel=1,unit="mm hr-1",&
             direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ESOIL,                     &
             value=edir,vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SUBSNOW,                   &
             value=(esnow/2.83E+6),vlevel=1,unit="kg m-2 s-1",direction="-",&
             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SUBSNOW,                   &
             value=((esnow/2.83E+6)*3600.0),vlevel=1,unit="mm hr-1",&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SUBSNOW,                   &
             value=esnow,vlevel=1,unit="W m-2",direction="-",surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CANOPINT,vlevel=1,         &
             value=noah36_struc(n)%noah(t)%cmc*LIS_CONST_RHOFW,           &
             unit="kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
        ! ROOT ZONE COLUMN SOIL MOISTURE IN METERS (SOILRZ)
        do k = 1,noah36_struc(n)%noah(t)%nroot
           soilrz    = soilrz + (noah36_struc(n)%noah(t)%smc(k)*sldpth(k) *    &
                                 LIS_CONST_RHOFW)
           soilrzmax = soilrzmax + (noah36_struc(n)%noah(t)%smcmax*sldpth(k) * &
                                    LIS_CONST_RHOFW) ! SY
        end do
        noah36_struc(n)%noah(t)%rootmoist = soilrz   !in volumetric
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROOTMOIST,vlevel=1,        &
             value=soilrz,unit="kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROOTMOIST,vlevel=1,        &
             value=soilrz/1000.0,unit="m^3 m-3",direction="-",surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CCOND,vlevel=1,            &
             value=RC,unit="m s-1",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CHFORC,                    &
             value=noah36_struc(n)%noah(t)%ch,vlevel=1,unit="m s-1",&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EMISSFORC,                 &
             value=noah36_struc(n)%noah(t)%emiss,vlevel=1,unit="-",&
             direction="-",surface_type=LIS_rc%lsm_index)
! J.Case (9/11/2014) -- Added support for outputting GVF in % for AWIPS II display.
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GREENNESS,                 &
             value=noah36_struc(n)%noah(t)%shdfac,vlevel=1,unit="-",             &
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GREENNESS,                 &
             value=noah36_struc(n)%noah(t)%shdfac*100.,vlevel=1,unit="%",        &
             direction="-",surface_type=LIS_rc%lsm_index)

! EMK Bug fix...Output deep soil temperature
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TEMPBOT,                   &
             value=noah36_struc(n)%noah(t)%tempbot,vlevel=1,unit="K",            &
             direction="-",surface_type=LIS_rc%lsm_index)
          
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_Z0BRD,                     &
             value=noah36_struc(n)%noah(t)%z0brd,vlevel=1,unit="m",&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROUGHNESS,                 &
             value=noah36_struc(n)%noah(t)%z0,vlevel=1,unit="m",&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LAI,                       &
             value=noah36_struc(n)%noah(t)%lai,vlevel=1,unit="-",&
             direction="-",surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CH,                        &
             value=noah36_struc(n)%noah(t)%CH,vlevel=1,unit="m s-1",direction="-",&
             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CM,                        &
             value=noah36_struc(n)%noah(t)%CM,vlevel=1,unit="m s-1",direction="-",&
             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_MIXRATIO,                  &
             value=noah36_struc(n)%noah(t)%Q1,vlevel=1,unit="kg kg-1",direction="-",&
             surface_type=LIS_rc%lsm_index)

     ! SY: Begin FLDAS
     WR_TimeStep = (etp/(((1.-SNCOVR)*LIS_CONST_LATVAP) +              &
          SNCOVR*CONST_LATSUB))*LIS_rc%ts
     AET_TimeStep = evp*LIS_rc%ts
     WRSI_TimeStep = LIS_rc%udef
     if( (WR_TimeStep > SUMWR_EXITWRSIBELOWTHISQTY) .and. &
         (WR_TimeStep < SUMWR_EXITWRSIABOVETHISQTY) .and. &
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
       IF( WR_TimeStep < 0.0 ) THEN
         WRSI_TimeStep = 100
       END IF
     END IF
     if(soilrzmax.gt.0) then 
        IF( 100*soilrz/soilrzmax < noah36_struc(n)%noah(t)%smcwlt ) THEN
           WRSI_TimeStep = 1
        END IF
     endif

     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WR_TimeStep,               &
          value=WR_TimeStep,vlevel=1,unit="kg m-2",direction="-",&
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_AET_TimeStep,              &
          value=AET_TimeStep,vlevel=1,unit="kg m-2",direction="-",&
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_WRSI_TimeStep,              &
          value=WRSI_TimeStep,vlevel=1,unit="-",direction="-",&
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SurplusWater_TimeStep,     &
          value= noah36_struc(n)%noah(t)%runoff1*LIS_CONST_RHOFW*dt +  &
          noah36_struc(n)%noah(t)%runoff2*LIS_CONST_RHOFW*dt,   &
          vlevel=1,unit="kg m-2",direction="-",&
          surface_type=LIS_rc%lsm_index)
     if(soilrzmax.gt.0) then 
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWI,vlevel=1,        &
             value=100*soilrz/soilrzmax,unit="%",direction="-",&
             surface_type=LIS_rc%lsm_index)
     else
        call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWI,vlevel=1,        &
             value=LIS_rc%udef,unit="%",direction="-",&
             surface_type=LIS_rc%lsm_index)
     endif
     
        ! SY: End FLDAS

        deallocate(sldpth)
!T2 and Q2 diagnostics 
#if (defined COUPLED)
     rho = noah36_struc(n)%noah(t)%psurf/(287.04*noah36_struc(n)%noah(t)%t1)
     q2diag =  noah36_struc(n)%noah(t)%q1 - noah36_struc(n)%noah(t)%qle/&
          (noah36_struc(n)%noah(t)%cqs2* rho)
     t2diag = noah36_struc(n)%noah(t)%t1 - noah36_struc(n)%noah(t)%qh/&
          (rho*noah36_struc(n)%noah(t)%chs2*CP)

     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_T2DIAG,                       &
          value=t2diag,vlevel=1,unit="K",direction="-",&
          surface_type=LIS_rc%lsm_index)
     call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_Q2DIAG,                       &
          value=q2diag,vlevel=1,unit="kg kg-1",direction="-",&
          surface_type=LIS_rc%lsm_index)
#endif

!Adjust deep soil temperature as lagged function of skin temperature
    if (LIS_rc%tbot_update_lag .eq. 1) then
!       write(LIS_logunit,*)'EMK: domain,tile,sncovr,tsoil,t1 = ', &
!            n,t,sncovr,tsoil,NOAH36_STRUC(N)%NOAH(T)%T1
       call LIS_updateTbot(n,t,julian_in,yr,dt,&
            (TSOIL * (1.0 - SNCOVR)) + (NOAH36_STRUC(N)%NOAH(T)%T1 * SNCOVR), &
            tbot)
       noah36_struc(n)%noah(t)%tempbot = tbot
    end if

!reset 
        noah36_struc(n)%noah(t)%tair = 0 
        noah36_struc(n)%noah(t)%qair = 0 
        noah36_struc(n)%noah(t)%swdown = 0 
        noah36_struc(n)%noah(t)%lwdown = 0
        noah36_struc(n)%noah(t)%uwind = 0
        noah36_struc(n)%noah(t)%vwind = 0 
        noah36_struc(n)%noah(t)%psurf = 0 
        noah36_struc(n)%noah(t)%rainf = 0 
        noah36_struc(n)%noah(t)%rainf_c = 0 
        
        if(noah36_struc(n)%forcing_ch .ne.0) then 
           noah36_struc(n)%noah(t)%ch = 0 
        endif
        
        if(LIS_FORC_GVF%selectOpt.eq.1) then 
           noah36_struc(n)%noah(t)%shdfac=0.0
        endif
        
        if(LIS_FORC_ALB%selectOpt.eq.1) then 
           noah36_struc(n)%noah(t)%alb=0.0
        endif
        
        if(LIS_FORC_Z0%selectOpt.eq.1) then 
           noah36_struc(n)%noah(t)%z0=0.0
        endif
        
     enddo
!$OMP END DO
!$OMP END PARALLEL 
     noah36_struc(n)%forc_count = 0 
  endif
end subroutine noah36_main

!
!!*** NOAH FUNCTIONS ****************************************************
!
!  FUNCTION DQS (T)
!    
!    use LIS_constantsMod, only : LIS_CONST_LATVAP, LIS_CONST_TKFRZ
!    IMPLICIT NONE
!    
!!
!!  PURPOSE:  TO CALCULATE VALUES OF VAPOR PRESSURE (E)
!!            AND P * DQS/DT (P TIMES CHG IN SAT MXG RATIO WITH RESPECT
!!            TO THE CHG IN TEMP) IN SUBSTITUTION TO THE LOOK-UP TABLES.
!!
!!            FORMULAS AND LIS_CONSTANTS FROM ROGERS AND YAU, 1989.
!!                         ADDED BY PABLO J. GRUNMANN, 6/30/97.
!!
!
!    REAL DESDT
!    REAL DQS
!!    REAL ESD
!    REAL LW
!    REAL T
!    REAL ES
!  
!  !      REAL, PARAMETER:: CP = 1005.
!  !      REAL, PARAMETER:: CV = 718.
!  !      REAL, PARAMETER:: CVV = 1410.
!    REAL, PARAMETER:: CPV = 1870.
!    REAL, PARAMETER:: RV = 461.5
!    REAL, PARAMETER:: CW = 4187.
!    REAL, PARAMETER:: EPS = 0.622
!    REAL, PARAMETER:: ESO = 611.2
!    REAL, PARAMETER:: TO = 273.15
!    REAL, PARAMETER:: LVH2O = 2.501000E+6
!  
!  
!!     ABOUT THE PARAMETERS:
!!
!!     EPS ---------- WATER - DRY AIR MOLECULAR MASS RATIO, EPSILON
!!
!!   VALUES FOR SPECIFIC HEAT CAPACITY AND INDIVIDUAL GAS LIS_CONSTANTS
!!   IN [JOULES/(KG*KELVIN)] UNITS.
!!
!!     DRY AIR:
!!             CP, CV
!!     WATER VAPOR:
!!                 CVV = 1410.
!!                 CPV = 1870.
!!                 RV  =  461.5
!!     LIQUID WATER:
!!                  CW = 4187.
!!
!!     ESO = ES(T=273.15 K) = SAT. VAPOR PRESSURE (IN PASCAL) AT T=TO
!!      TO = 273.15
!!
!!     SAT. MIXING  RATIO: QS ~= EPS*ES/P
!!     CLAUSIUS-CLAPEYRON: DES/DT = L*ES/(RV*T^2)
!!     @QS/@T =  (EPS/P)*DES/DT
!
!    LW = LIS_CONST_LATVAP - ( CW - CPV ) * ( T - LIS_CONST_TKFRZ )
!    ES = ESO*EXP (LW*(1/LIS_CONST_TKFRZ - 1/T)/RV)
!    DESDT = LW*ES/(RV*T*T)
!  
!  !    FOR INSERTION IN DQSDT FUNCTION:
!  !    DQSDT = DQS/P , WHERE DQS = EPS*DESDT
!  
!    DQS = EPS*DESDT
!   
!    RETURN
!  END FUNCTION DQS
