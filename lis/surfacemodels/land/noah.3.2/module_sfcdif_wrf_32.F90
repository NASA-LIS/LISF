!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! #define _BARLAGE_VERSION_
#define _KWM_VERSION_

#ifdef _KWM_VERSION_
module module_sfcdif_wrf_32
  use module_model_constants_32

  integer, parameter, private :: ITRMX = 5

  real,    parameter, private :: EXCML   = 0.0001
  real,    parameter, private :: EXCMS   = 0.0001
  real,    parameter, private :: VKARMAN = 0.4
  real,    parameter, private :: ZTFC    = 1.0
  real,    parameter, private :: ELOCP   = 2.72E6 / CP
  real,    parameter, private :: EPSU2   = 1.E-6
  real,    parameter, private :: EPSUST  = 1.e-9
  real,    parameter, private :: SQVISC  = 258.2
  real,    parameter, private :: RIC     = 0.505
  real,    parameter, private :: EPSZT   = 1.E-28

  integer, parameter, private :: KZTM  = 10001
  integer, parameter, private :: KZTM2 = KZTM-2

  real,    parameter, private :: WWST  = 1.2
  real,    parameter, private :: WWST2 = WWST*WWST
  real,    parameter, private :: ztmin2 = -5.0

contains

  SUBROUTINE SFCDIF_MYJ ( ZSL, Z0, Z0BASE, SFCPRS, TZ0, TLOW, QZ0,   &
       QLOW, SFCSPD, CZIL, RIB, AKMS, AKHS, VEGTYP, ISURBAN, IZ0TLND,&
       ZTMAX2, DZETA2, PSIH2, PSIM2, USTAR)
!     ****************************************************************
!     *                                                              *
!     *                       SURFACE LAYER                          *
!     *                                                              *
!     ****************************************************************
!----------------------------------------------------------------------
!
    IMPLICIT NONE
!
!----------------------------------------------------------------------
!
    REAL, INTENT(IN) :: ZSL  ! Height above ground (m) of atmospheric variables.
    REAL, INTENT(IN) :: Z0   ! Roughness length (m); may be adjusted for snow
    REAL, INTENT(IN) :: Z0BASE ! Background roughness length (m)
    REAL, INTENT(IN) :: SFCPRS
    REAL, INTENT(IN) :: TZ0  ! Temperature (K) at level Z0 above ground -- assumed to be the Skin Temperature.
    REAL, INTENT(IN) :: TLOW ! Temperature (K) at level <Za> above ground.
    REAL, INTENT(IN) :: QZ0  ! Specific humidity (kg kg-1) at level Z0 above ground
    REAL, INTENT(IN) :: QLOW ! Specific humidity (kg kg-1) at level ZSL above ground
    REAL, INTENT(IN) :: SFCSPD
    REAL, INTENT(IN) :: CZIL
    INTEGER, INTENT(IN) :: VEGTYP
    INTEGER, INTENT(IN) :: ISURBAN
    INTEGER, INTENT(IN) :: IZ0TLND
!
    REAL,INTENT(INOUT) :: AKHS
    REAL,INTENT(INOUT) :: AKMS
    REAL,INTENT(OUT)   :: RIB

    REAL, INTENT(IN) :: ZTMAX2
    REAL, INTENT(IN) :: DZETA2
! These need to be initialized somewhere.
  real, dimension(KZTM), INTENT(IN) :: PSIH2
  real, dimension(KZTM), INTENT(IN) :: PSIM2

!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
    REAL :: THLOW ! Potential Temperature (K) at level ZSL above ground.
    REAL :: THZ0  ! Potential Temperature (K) at level Z0 above ground.
    REAL :: THELOW !
    REAL :: CWMLOW ! Cloud water (assumed to be zero for HRLDAS)
    REAL :: USTAR
    REAL :: RLMO
    INTEGER :: ITR,K
    REAL :: CZIL_LOCAL
    REAL :: ZILFC
!
    REAL :: A,B,BTGH,BTGX,CXCHL,CXCHS,CZETMAX,DTHV,DU2,ELFC      &
         & ,PSHZ,PSHZL,PSMZ,PSMZL   &
         & ,RDZT                                  &
         & ,RLOGT,RLOGU,RZ,RZST,RZSU,SIMH,SIMM,TEM,THM          &
         & ,USTARK,WSTAR2               &
         & ,ZETALT,ZETALU          &
         & ,ZETAT,ZETAU,ZQ,ZSLT,ZSLU,ZT,ZU,ZZIL

!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Compute potential temperatures from input <Za>-level temperature TLOW
! and input Skin-level (or Z0-level) temperature.

    THLOW = TLOW * (P0/SFCPRS) ** RCP
    THZ0  = TZ0  * (P0/SFCPRS) ** RCP  ! Should actually use a pressure at Z0?

    ! For THELOW, WRF uses (1-QC*ELOCP/T)*TH.  If QC is assumed to be zero,
    ! Then THELOW is just TH ?
    THELOW = THLOW
!----------------------------------------------------------------------

    CXCHL=EXCML/ZSL
 !
    BTGX=G/THLOW
    ELFC=VKARMAN*BTGX

!MBKWM    IF(PBLH>1000.)THEN
!MBKWM       BTGH=BTGX*PBLH
!MBKWM    ELSE
    BTGH=BTGX*1000.
!MBKWM    ENDIF

!----------------------------------------------------------------------

    THM=(THELOW+THZ0)*0.5
    TEM=(TLOW+TZ0)*0.5
!
    A=THM*P608
    B=(ELOCP/TEM-1.-P608)*THM
    CWMLOW = 0.0 ! MB/KWM Assume no cloud water
!
    DTHV=((THELOW-THZ0)*((QLOW+QZ0+CWMLOW)*(0.5*P608)+1.)          &
         &        +(QLOW-QZ0+CWMLOW)*A+CWMLOW*B)
!
!    DU2=MAX((ULOW-UZ0)**2+(VLOW-VZ0)**2,EPSU2)
    DU2 = MAX(SFCSPD*SFCSPD,EPSU2)
    RIB=BTGX*DTHV*ZSL/DU2
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    ZU=Z0
    ZT=ZU*ZTFC
    ZSLU=ZSL+ZU
    RZSU=ZSLU/ZU
    RLOGU=LOG(RZSU)
    ZSLT=ZSL+ZU ! u,v and t are at the same level

    IF ( (IZ0TLND==0) .or. (VEGTYP == ISURBAN) ) THEN
       ! Just use the incoming CZIL value.
       ZILFC=-CZIL*VKARMAN*SQVISC
    ELSE
       ! Modify CZIL according to Chen & Zhang, 2009
       ! CZIL = 10 ** -0.40 H, ( where H = 10*Zo )
       CZIL_LOCAL = 10.0 ** ( -0.40 * ( Z0 / 0.07 ) )
       ZILFC=-CZIL_LOCAL*VKARMAN*SQVISC
    ENDIF

!----------------------------------------------------------------------
!
!  RIB modification to ZILFC term
!
!      CZETMAX = 20.
    CZETMAX = 10.
! stable
    IF(DTHV>0.)THEN
       IF (RIB<RIC) THEN
          ZZIL=ZILFC*(1.0+(RIB/RIC)*(RIB/RIC)*CZETMAX)
       ELSE
          ZZIL=ZILFC*(1.0+CZETMAX)
       ENDIF
! unstable
    ELSE
       ZZIL=ZILFC
    ENDIF

!----------------------------------------------------------------------
    ! First guess of USTAR based on input (previous timestep) AKHS and AKMS?
    IF (BTGH * AKHS * DTHV .ne. 0.0) THEN
       WSTAR2 = WWST2* ABS (BTGH * AKHS * DTHV)** (2./3.)
    ELSE
       WSTAR2 = 0.0
    END IF
    USTAR = MAX (SQRT (AKMS * SQRT (DU2+ WSTAR2)),EPSUST)
!----------------------------------------------------------------------

!
!----------------------------------------------------------------------
!
    land_point_iteration: DO ITR=1,ITRMX
!
!----------------------------------------------------------------------
!***  ZILITINKEVITCH FIX FOR ZT
!----------------------------------------------------------------------
!
       ZT=MAX(EXP(ZZIL*SQRT(USTAR*Z0BASE))*Z0BASE,EPSZT)
!
       RZST=ZSLT/ZT
       RLOGT=LOG(RZST)
!
!----------------------------------------------------------------------
!***  1./MONIN-OBUKHOV LENGTH-SCALE
!----------------------------------------------------------------------
!
       RLMO=ELFC*AKHS*DTHV/USTAR**3
       ZETALU=ZSLU*RLMO
       ZETALT=ZSLT*RLMO
       ZETAU=ZU*RLMO
       ZETAT=ZT*RLMO
!
       ZETALU=MIN(MAX(ZETALU,ZTMIN2),ZTMAX2)
       ZETALT=MIN(MAX(ZETALT,ZTMIN2),ZTMAX2)
       ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)
       ZETAT=MIN(MAX(ZETAT,ZTMIN2/RZST),ZTMAX2/RZST)
!
!----------------------------------------------------------------------
!***  LAND FUNCTIONS
!----------------------------------------------------------------------
!
       RZ=(ZETAU-ZTMIN2)/DZETA2
       K=INT(RZ)
       RDZT=RZ-REAL(K)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSMZ=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
       RZ=(ZETALU-ZTMIN2)/DZETA2
       K=INT(RZ)
       RDZT=RZ-REAL(K)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSMZL=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
       SIMM=PSMZL-PSMZ+RLOGU
!
       RZ=(ZETAT-ZTMIN2)/DZETA2
       K=INT(RZ)
       RDZT=RZ-REAL(K)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSHZ=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
       RZ=(ZETALT-ZTMIN2)/DZETA2
       K=INT(RZ)
       RDZT=RZ-REAL(K)
       K=MIN(K,KZTM2)
       K=MAX(K,0)
       PSHZL=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
       SIMH=(PSHZL-PSHZ+RLOGT)
!----------------------------------------------------------------------
       USTARK=USTAR*VKARMAN
       AKMS=MAX(USTARK/SIMM,CXCHL)
       AKHS=MAX(USTARK/SIMH,CXCHL)
!
!----------------------------------------------------------------------
!***  BELJAARS CORRECTION FOR USTAR
!----------------------------------------------------------------------
!
       IF(DTHV<=0.)THEN                                           !zj
          WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)               !zj
       ELSE                                                       !zj
          WSTAR2=0.                                               !zj
       ENDIF                                                      !zj
!
       USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)

!----------------------------------------------------------------------
    ENDDO land_point_iteration
!----------------------------------------------------------------------

  END SUBROUTINE SFCDIF_MYJ

  SUBROUTINE MYJSFCINIT(ZTMAX2,DZETA2,PSIH2,PSIM2)
!----------------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------------
    INTEGER :: K
!
    REAL, INTENT(OUT) :: ZTMAX2,DZETA2
  real, dimension(KZTM), INTENT(OUT) :: PSIH2
  real, dimension(KZTM), INTENT(OUT) :: PSIM2
    REAL :: X,ZETA1,ZETA2,ZRNG1,ZRNG2,ZTMAX1,DZETA1
!
    REAL, parameter :: PIHF=3.1415926/2.
    real, parameter :: EPS=1.E-6
    real, parameter :: ztmin1 = -5.0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!***  COMPUTE SURFACE LAYER INTEGRAL FUNCTIONS
!
!----------------------------------------------------------------------

      ZTMAX1=1.0
      ZTMAX2=1.0

      ZRNG1=ZTMAX1-ZTMIN1
      ZRNG2=ZTMAX2-ZTMIN2

      DZETA1=ZRNG1/(KZTM-1)
      DZETA2=ZRNG2/(KZTM-1)
!
!----------------------------------------------------------------------
!***  FUNCTION DEFINITION LOOP
!----------------------------------------------------------------------
!
      ZETA1=ZTMIN1
      ZETA2=ZTMIN2
!
      DO K=1,KZTM
!
!----------------------------------------------------------------------
!***  UNSTABLE RANGE
!----------------------------------------------------------------------
!
         IF(ZETA2<0.)THEN
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
!
            X=SQRT(SQRT(1.-16.*ZETA2))
!
            PSIM2(K)=-2.*LOG((X+1.)/2.)-LOG((X*X+1.)/2.)+2.*ATAN(X)-PIHF
            PSIH2(K)=-2.*LOG((X*X+1.)/2.)
!----------------------------------------------------------------------
!***  STABLE RANGE
!----------------------------------------------------------------------
!
         ELSE
!
!----------------------------------------------------------------------
!***  HOLTSLAG AND DE BRUIN 1988
!----------------------------------------------------------------------
!
            PSIM2(K)=0.7*ZETA2+0.75*ZETA2*(6.-0.35*ZETA2)*EXP(-0.35*ZETA2)
            PSIH2(K)=0.7*ZETA2+0.75*ZETA2*(6.-0.35*ZETA2)*EXP(-0.35*ZETA2)
!----------------------------------------------------------------------
!
         ENDIF
!
!----------------------------------------------------------------------
         IF(K==KZTM)THEN
            ZTMAX1=ZETA1
            ZTMAX2=ZETA2
         ENDIF
!
         ZETA1=ZETA1+DZETA1
         ZETA2=ZETA2+DZETA2
!----------------------------------------------------------------------
      ENDDO
!----------------------------------------------------------------------
      ZTMAX1=ZTMAX1-EPS
      ZTMAX2=ZTMAX2-EPS
!----------------------------------------------------------------------
!
    END SUBROUTINE MYJSFCINIT

END MODULE MODULE_SFCDIF_WRF_32
#endif

#ifdef _BARLAGE_VERSION_
module module_sfcdif_wrf_32
  use module_model_constants_32
contains
  subroutine myjsfcinit()
    implicit none
  end subroutine myjsfcinit

!----------------------------------------------------------------------
  SUBROUTINE SFCDIF_MYJ ( ZSL , Z0 , Z0BRD , SFCPRS , TZ0 , TLOW , QZ0 , &
       QLOW , SFCSPD , RIB , AKMS , AKHS , VEGTYP , ISURBAN , IZ0TLND )
!----------------------------------------------------------------------
!
    IMPLICIT NONE
!
!----------------------------------------------------------------------
!
    REAL,    INTENT(IN) :: ZSL    ! Atmo height of temperature, wind, etc
    REAL,    INTENT(IN) :: Z0     ! Roughness length, can contain snow effects
    REAL,    INTENT(IN) :: TZ0    ! Temperature at heigh z0, assumed skin T
    REAL,    INTENT(IN) :: TLOW   ! Atmospheric temperature
    REAL,    INTENT(IN) :: SFCSPD ! Atmospheric wind speed
    REAL,    INTENT(IN) :: Z0BRD  ! Background z0 of vegetation
    REAL,    INTENT(IN) :: QLOW   ! Atmospheric specific humidity
    REAL,    INTENT(IN) :: QZ0    ! Surface specific humidity
    REAL,    INTENT(IN) :: SFCPRS ! Atmospheric pressure (model level)
    INTEGER, INTENT(IN) :: VEGTYP
    INTEGER, INTENT(IN) :: ISURBAN
    INTEGER, INTENT(IN) :: IZ0TLND
!
    REAL,    INTENT(INOUT) :: AKHS ! Heat exchange coefficient
    REAL,    INTENT(INOUT) :: AKMS ! Momentum exchange coefficient
    REAL,    INTENT(OUT)   :: RIB
!----------------------------------------------------------------------
#ifdef _ITERATION50_
    INTEGER :: ITRMX=50 ! Iteration count for sfc layer computations
#else
    INTEGER :: ITRMX=5  ! Iteration count for sfc layer computations
#endif
!
    REAL, PARAMETER :: PI = 3.1415926    ! Pi, used in MO functions
    REAL, PARAMETER :: EXCML = 0.0001    ! Used to calculate min exchange coeff
    REAL, PARAMETER :: VKARMAN = 0.4     ! von Karman constant
    REAL, PARAMETER :: PBLH = 1000.0     ! Assume PBL height equals 1000m
    REAL, PARAMETER :: ELOCP = 2.72E6/CP ! Lv/cp; used for virtual temp calc
    REAL, PARAMETER :: EPSU2 = 1.E-6     ! Minimum wind speed
    REAL, PARAMETER :: EPSUST = 1.e-9    ! Minimum ustar
    REAL, PARAMETER :: ZTFC = 1.         ! Initial assumption of z0T/z0m
    REAL, PARAMETER :: TOPOFAC = 9.0e-6  ! Stable BL topo adjustment factor
    REAL, PARAMETER :: EPSZT = 1.E-28    ! minimum heat roughness length
    REAL, PARAMETER :: ZTMIN = -5.0      ! minimum z/L (MO similarity parameter)
    REAL, PARAMETER :: ZTMAX = 1.0       ! maximum z/L (MO similarity parameter)
    REAL, PARAMETER :: WWST = 1.2        ! Constant for ustar convective correction
    REAL, PARAMETER :: WWST2 = WWST*WWST ! Square WWST
!    REAL, PARAMETER :: CZIL = 0.1        ! Zilitinkevitch coefficient
    REAL, PARAMETER :: SQVISC = 258.2    ! 1/sqrt(visc) from Reynolds number
!    REAL, PARAMETER :: ZILFC = -CZIL*VKARMAN*SQVISC ! Zilitinkevitch factor
    REAL, PARAMETER :: P608 = 461.6/287.0-1.0 ! Rv/Rd - 1
    REAL, PARAMETER :: RIC = 0.505       ! Critical Richardson number
    REAL :: CZIL
    REAL :: ZILFC

    INTEGER :: ITR,K
!
    REAL :: RDZ,PSIMU,PSIHU,PSIMS,PSIHS,CXCHL,CXCHS,BTGH,BTGX,ELFC,DTHV, &
         DU2,WSTAR2,USTAR,ZT,ZU,ZSLU,RZSU,RLOGU,ZSLT,TOPOTERM,ZZIL,   &
         RZST,RLOGT,RLMO,ZETALT,ZETALU,ZETAT,ZETAU,XLU4,XLT4,XU4,XT4, &
         XLT,XLU,XT,XU,PSMZ,SIMM,SIMH,USTARK,X,PSHZ,THLOW,THZ0,  &
         THM,TEM,A,B,CWMLOW,CZETMAX

!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  UNSTABLE RANGE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
    PSIMU(X)=-2.*LOG((X+1.)/2.)-LOG((X*X+1.)/2.)+2.*ATAN(X)-PI/2.0
    PSIHU(X)=-2.*LOG((X*X+1.)/2.)
!----------------------------------------------------------------------
!***  STABLE RANGE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  HOLTSLAG AND DE BRUIN 1988
!----------------------------------------------------------------------
!
    PSIMS(X)=0.7*X+0.75*X*(6.-0.35*X)*EXP(-0.35*X)
    PSIHS(X)=0.7*X+0.75*X*(6.-0.35*X)*EXP(-0.35*X)
!----------------------------------------------------------------------

    THLOW = TLOW * (100000.0/SFCPRS) ** (2.0/7.0)
! approx THZ0 since no sfc pressure
    THZ0 = TZ0 * (100000.0/SFCPRS) ** (2.0/7.0) 

    RDZ=1./ZSL
    CXCHL=EXCML*RDZ
!
    BTGX=G/THLOW
    ELFC=VKARMAN*BTGX
!
!      IF(PBLH>1000.)THEN
!        BTGH=BTGX*PBLH
!      ELSE
    BTGH=BTGX*1000.
!      ENDIF
!
!----------------------------------------------------------------------
    THM=(THLOW+THZ0)*0.5
    TEM=(TLOW+TZ0)*0.5
!
    A=THM*P608
    B=(ELOCP/TEM-1.-P608)*THM
    CWMLOW = 0.0
!

! from WRF MYJSFC code
    DTHV=((THLOW-THZ0)*((QLOW+QZ0+CWMLOW)*(0.5*P608)+1.)          &
         &        +(QLOW-QZ0+CWMLOW)*A+CWMLOW*B) 

    DU2=MAX(SFCSPD*SFCSPD,EPSU2)
    IF (BTGH * AKHS * DTHV .ne. 0.0) THEN
       WSTAR2 = WWST2* ABS (BTGH * AKHS * DTHV)** (2./3.)
    ELSE
       WSTAR2 = 0.0
    END IF
    USTAR = MAX (SQRT (AKMS * SQRT (DU2+ WSTAR2)),EPSUST)
!----------------------------------------------------------------------
    ZU=Z0
    ZT=ZU*ZTFC
    ZSLU=ZSL+ZU
    RZSU=ZSLU/ZU
    RLOGU=LOG(RZSU)

    ZSLT=ZSL+ZU ! u,v and t are at the same level

    IF ( (IZ0TLND==0) .or. (VEGTYP == ISURBAN) ) THEN
       ! Just use the original CZIL value.
       CZIL = 0.1
    ELSE
       ! Modify CZIL according to Chen & Zhang, 2009
       ! CZIL = 10 ** -0.40 H, ( where H = 10*Zo )
       CZIL = 10.0 ** ( -0.40 * ( Z0 / 0.07 ) )
    ENDIF
    ZILFC = -CZIL*VKARMAN*SQVISC
!----------------------------------------------------------------------
!
!  RIB modification to ZILFC term
!
!       CZETMAX = 20.
    CZETMAX = 10.
#ifdef _CZETMAX20_
    CZETMAX = 20.
#endif

    RIB=BTGX*DTHV*ZSL/DU2
! stable
    IF(DTHV>0.)THEN
!          ZZIL=ZILFC*TOPOTERM
       IF (RIB<RIC) THEN
          ZZIL=ZILFC*(1.0+(RIB/RIC)*(RIB/RIC)*CZETMAX)
       ELSE
          ZZIL=ZILFC*(1.0+CZETMAX)
       ENDIF
! unstable
    ELSE
       ZZIL=ZILFC
    ENDIF
!      write(12,'(8f10.5)') dthv,thlow,btgx,zsl,du2,rib,zzil,zilfc
!----------------------------------------------------------------------
!
    land_point_iteration: DO ITR=1,ITRMX
!
!----------------------------------------------------------------------
!***  ZILITINKEVITCH FIX FOR ZT
!----------------------------------------------------------------------
!
! oldform   ZT=MAX(EXP(ZZIL*SQRT(USTAR*ZU))*ZU,EPSZT)
       ZT=MAX(EXP(ZZIL*SQRT(USTAR*Z0BRD))*Z0BRD,EPSZT)
!
       RZST=ZSLT/ZT
       RLOGT=LOG(RZST)
!
!----------------------------------------------------------------------
!***  1./MONIN-OBUKHOV LENGTH-SCALE
!----------------------------------------------------------------------
!
       RLMO=ELFC*AKHS*DTHV/USTAR**3
       ZETALU=ZSLU*RLMO
       ZETALT=ZSLT*RLMO
       ZETAU=ZU*RLMO
       ZETAT=ZT*RLMO
!
       ZETALU=MIN(MAX(ZETALU,ZTMIN),ZTMAX)
       ZETALT=MIN(MAX(ZETALT,ZTMIN),ZTMAX)
       ZETAU=MIN(MAX(ZETAU,ZTMIN/RZSU),ZTMAX/RZSU)
       ZETAT=MIN(MAX(ZETAT,ZTMIN/RZST),ZTMAX/RZST)
!
!----------------------------------------------------------------------
!***  LAND FUNCTIONS
!----------------------------------------------------------------------
!
       IF (RLMO .lt. 0.)THEN
          XLU4 = 1. -16.* ZETALU
          XLT4 = 1. -16.* ZETALT
          XU4 = 1. -16.* ZETAU

          XT4 = 1. -16.* ZETAT
          XLU = SQRT (SQRT (XLU4))
          XLT = SQRT (SQRT (XLT4))
          XU = SQRT (SQRT (XU4))

          XT = SQRT (SQRT (XT4))
          PSMZ = PSIMU (XU)
          SIMM = PSIMU (XLU) - PSMZ + RLOGU
          PSHZ = PSIHU (XT)
          SIMH = PSIHU (XLT) - PSHZ + RLOGT
       ELSE
          PSMZ = PSIMS (ZETAU)
          SIMM = PSIMS (ZETALU) - PSMZ + RLOGU
          PSHZ = PSIHS (ZETAT)
          SIMH = PSIHS (ZETALT) - PSHZ + RLOGT
       END IF
!----------------------------------------------------------------------
       USTARK=USTAR*VKARMAN
       AKMS=MAX(USTARK/SIMM,CXCHL)
       AKHS=MAX(USTARK/SIMH,CXCHL)
!
!----------------------------------------------------------------------
!***  BELJAARS CORRECTION FOR USTAR
!----------------------------------------------------------------------
!
       IF(DTHV<=0.)THEN                                           !zj
          WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)                !zj
       ELSE                                                       !zj
          WSTAR2=0.                                                !zj
       ENDIF                                                      !zj
!
       USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
!
!----------------------------------------------------------------------
    ENDDO land_point_iteration
!----------------------------------------------------------------------
!
  END SUBROUTINE SFCDIF_MYJ
! ----------------------------------------------------------------------
end module module_sfcdif_wrf_32
#endif
