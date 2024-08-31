!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! $Id: surfacelayer.F90,v 1.2.4.1.10.3.22.1.2.6.2.1.2.1.2.1.2.5.6.3 2011-06-10 19:26:22 amolod Exp $
! $Name: reichle-LDASsa_m2-10 $

module clsmf25_surfacelayer

! !USES:

use clsmf25_MAPL_constants
use clsmf25_DragCoefficientsMod

implicit none
private
PUBLIC louissurface,z0sea,helfsurface,psi,phi,linadj,zcsub

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: louissurface
! !INTERFACE:
      subroutine louissurface(ISTYPE,N,UU,WW,PS,TA,TS,QA,QS,PCU,LAI, &
                               Z0,DZ,CM,CN,RI,ZT,ZQ,CH,CQ,UUU,UCN,RE,DCH,DCQ)
      integer,           intent(IN ) :: N
      integer,           intent(IN ) :: ISTYPE
      real,    intent(IN)    :: UU (:)
      real,    intent(IN)    :: WW (:,:)
      real,    intent(IN)    :: PS (:)
      real,    intent(IN)    :: TA (:)
      real,    intent(IN)    :: TS (:,:)
      real,    intent(IN)    :: QA (:)
      real,    intent(IN)    :: QS (:,:)
      real,    intent(IN)    :: PCU(:)
      real,    intent(IN)    :: LAI(:)
      real,    intent(INOUT) :: Z0 (:,:)
      real,    intent(IN)    :: DZ (:)
      real,    intent(INOUT) :: CM (:,:)
      real,    intent(OUT)   :: CN (:)
      real,    intent(OUT)   :: RI (:)
      real,    intent(OUT)   :: ZT (:)
      real,    intent(OUT)   :: ZQ (:)
      real,    intent(OUT)   :: CH (:,:)
      real,    intent(OUT)   :: CQ (:,:)
      real,    intent(OUT)   :: UUU (:)
      real,    intent(OUT)   :: UCN (:)
      real,    intent(OUT)   :: RE  (:)
      real,    optional, intent(OUT)   :: DCH (:,:)
      real,    optional, intent(OUT)   :: DCQ (:,:)


      real, parameter :: OCEANICEZ0 = 1.0e-3
      real, parameter :: LAKEZ0     = 1.0E-5
      real, parameter :: LAKEICEZ0  = 0.01
      real, parameter :: LANDICEZ0  = 0.002
      integer, parameter :: ICE   = 1
      integer, parameter :: WATER = 2
      integer,parameter  :: FSNW=4  !  Snowcover subtile
      integer  :: irun

      real, allocatable  :: UST(:)
      real, allocatable  :: TVS(:)
      real, allocatable  :: URA(:)
      real, allocatable  :: TVA(:)

      real :: LAICRT

      irun = size(UU,1)
      LAICRT = 0.3

      allocate( UST(1:irun)  )
      allocate( TVS(1:irun)  )
      allocate( URA(1:irun)  )
      allocate( TVA(1:irun)  )
      
      if(ISTYPE==1 .or. ISTYPE==3) then
        UCN = PCU*8640. !  cm/day
        UCN = ((19.8*UCN*UCN)/(1.5 + UCN + UCN*UCN))**0.4 !  REDELSPERGER et al. (2000)
        UCN = UCN*UCN
!       UCN = 0.0  ! Remove Redelsperger
      endif
!
      TVA = TA*( 1.0 + MAPL_VIREPS*QA)
      if(ISTYPE==3)TVA = TA*(1+MAPL_VIREPS*QA)
!
      if(ISTYPE==1 .or. ISTYPE==3) then
        UUU = max( sqrt(UU*UU + WW(:,N) + UCN), MAPL_USMIN )
      elseif (ISTYPE==2 .or. ISTYPE==4) then
        UUU = max(UU,MAPL_USMIN)
      endif

      if(ISTYPE==1) then
      URA = UUU*( PS/(MAPL_RGAS*TVA) )
      elseif (ISTYPE==2 .or. ISTYPE==3 .or. ISTYPE==4) then
      URA = UUU*PS/(MAPL_RGAS*TVA)
      endif
      TVS = TS(:,N)*(1+MAPL_VIREPS*QS(:,N))

!  Estimate of friction velocity from previous step's drag coefficient
!---------------------------------------------------------------------

      UST = UUU*sqrt(CM(:,N)/URA)
      UST = min( max(UST,1.e-6) , 5.0 )

!  Iterate for aerodynamic roughness, based on previous step's stability
!-----------------------------------------------------------------------

      if(ISTYPE==1) then
        if (N==WATER) then
          call Z0SEA (Z0(:,N),UST,DZ)
        else
          Z0(:,N) = OCEANICEZ0
        end if
      elseif (ISTYPE==2) then
        if(N==ICE) then
          Z0(:,N) = LAKEICEZ0
        else
          Z0(:,N) = LAKEZ0
        end if
      elseif (ISTYPE==4) then
        Z0(:,N) = LANDICEZ0
      end if

! Get a new drag coefficient for the updated roughness
!-----------------------------------------------------

      call GETCDM(TVA,UUU,DZ,TVS,Z0(:,N), CM(:,N),CN,RI)
      if(ISTYPE==1) UST = UUU*sqrt(CM(:,N))

!  Reynolds number at the top of the viscous sublayer
!----------------------------------------------------

      if(ISTYPE==1) then
        RE  = max(min(Z0(:,N)*UST/MAPL_NUAIR, 50.),0.1)
      elseif (ISTYPE==3) then
        RE = max(min(Z0(:,N)*(sqrt(CM(:,N))*UUU) / MAPL_NUAIR, 1000.),0.4)
      elseif (ISTYPE==2 .or. ISTYPE==4) then
        RE = min(Z0(:,N) * (sqrt(CM(:,N))*UUU) / MAPL_NUAIR, 1000.)
      end if


!  Calculate scalar roughnesses
!------------------------------

      if(ISTYPE==1) then
        if(N==WATER) then
! This is based on Zilitinkevich et al.
!         ZT = Z0(:,N) * exp(MAPL_KARMAN*(3.2 - 4.0*sqrt(max(RE,0.1))))
! This is based on Garrett 1992
!        ZT = Z0(:,N) * exp(2.0 - 2.48*RE**.25)
! This is as used in ecmwf model, based on:
!     Beljaars, A. C. M., 1994: 
!    The parametrization of surface fluxes in large-scale models under
!    free convection. Q. J. R. Meteorol. Soc., 121, 255-270.
!
          ZT = Z0(:,N) * (0.4/RE)
        else
! This is used over snow and bare soil in CSM
          ZT = Z0(:,N) * exp(-0.13*sqrt(RE))
        end if
      elseif (ISTYPE==2) then
        if(N==ICE) then 
          ZT = Z0(:,N) * exp(-0.13*sqrt(RE))
        else
          ZT = Z0(:,N) * (0.4/RE)
        end if
      elseif (ISTYPE==3) then
! Viscous sub-layer effects on Z0 vanish as the surface
!  becomes very vegetated (LAI -> infinity), assuming they are handled
!  in the vegetation resistance to heat and moisture fluxes.

!!!  Note: Comment out special code for SNOW condition to alleviate excessive night-time cooling
!!!  -------------------------------------------------------------------------------------------
!!!     if (N==FSNW) then
!!!       ZT = Z0(:,N)*exp( -0.13 * sqrt(RE) )
!!!       ZQ = ZT
!!!     else
          ZT = 1.0 - (1.0-exp(-0.13*sqrt(RE)))*exp(-max(LAI,LAICRT))
          where(LAI<LAICRT)
!!  Brutsaert, W., 1982: Evaporation into the Atmosphere. D. Reidel, 299pp.
            ZT = ZT + (exp( -2.5*sqrt(sqrt(RE)) + log(7.4) ) - ZT)*(LAICRT-LAI)/LAICRT
          end where

          ZT = Z0(:,N)*ZT
          ZQ = ZT
!!!     end if
      elseif (ISTYPE==4) then
        ZT = Z0(:,N)*exp( -0.13 * sqrt(RE) )
      endif

      if (ISTYPE==1) then
        ZQ = ZT * 1.5
      elseif (ISTYPE==2 .or. ISTYPE==4) then
        ZQ = ZT
      endif

!  Compute the the heat and moisture surface coefficients
!--------------------------------------------------------

      if (present(DCH) .and. present(DCQ)) then

         call GETCDH(TVA,UUU,DZ,TVS,ZT,ZQ,CN,RI,  CH(:,N),CQ(:,N),DCH(:,N),DCQ(:,N))

      else

         call GETCDH(TVA,UUU,DZ,TVS,ZT,ZQ,CN,RI,  CH(:,N),CQ(:,N))

      end if

      CH(:,N) = CH(:,N)*URA
      CM(:,N) = CM(:,N)*URA
      CQ(:,N) = CQ(:,N)*URA

      deallocate( UST )
      deallocate( TVS )
      deallocate( URA )
      deallocate( TVA )

end subroutine louissurface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Z0SEA - Computes $z_o$ as a function of $u^*$ over water surfaces
! !INTERFACE:

  subroutine Z0SEA (Z0, USTAR, DZ)

! !ARGUMENTS:

    real,    intent(INOUT) :: Z0   (:)
    real,    intent(INOUT) :: USTAR(:)
    real,    intent(IN   ) :: DZ   (:)

! !DESCRIPTION:
!        Compute roughness length for ocean points
!          based on functions of Large and Pond
!          and of Kondo

    real               :: UF(size(Z0))
    integer, parameter :: NumIter = 3
    integer            :: K


! Begin

! UF is Karman*|U|*sqrt[f_m(Ri)] = |U|*sqrt[C_D]*log(DZ/Z0+1). 
! The stability factor is from the previous time step.
!------------------------------------------------------------------------------------

    UF    = USTAR * ALOG(DZ/Z0 + 1.0)

! Iterate roughness and U*
!-------------------------

    ITERATION: do K=1,NumIter

! Evaluate z0 using Beljaars' approx.
!------------------------------------

       Z0 = (0.11*MAPL_NUAIR)/USTAR + (0.018/MAPL_GRAV)*USTAR**2

! Update U* including stability effects
!--------------------------------------

       USTAR = min(max(UF/ALOG( DZ/Z0 + 1.0 ) , 1.e-6),5.0)

    end do ITERATION

  end subroutine Z0SEA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: helfsurface
! !INTERFACE:
       SUBROUTINE helfsurface(VUS,VVS,VT1,VT2,VSH1,VSH2,VP,VPE, &
        VZ0,LAI,IVWATER,VHS,N,IRUN, &
        VRHO,VKH,VKM,VUSTAR,VXX,VYY,VCU,VCT,VRIB,VZETA,VWS, &
        t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0)
!**********************************************************************
!  SUBROUTINE helfsurface - COMPUTES SURFACE TRANSFER COEFFICIENTS
!
!   ARGUMENTS ::
!
!     INPUT:
!     ------
!    US            -         U - COMPONENT OF SURFACE WIND
!    VS            -         V - COMPONENT OF SURFACE WIND
!    THV1          -         VIRTUAL POTENTIAL TEMPERATURE AT NLAY
!    THV2          -         VIRTUAL POTENTIAL TEMPERATURE AT GROUND
!    TH1           -         POTENTIAL TEMPERATURE AT NLAY
!    TH2           -         POTENTIAL TEMPERATURE AT GROUND
!    SH1           -         SPECIFIC HUMIDITY AT NLAY
!    SH2           -         SPECIFIC HUMIDITY AT GROUND
!    PK            -         EVEN LEVEL PRESSURE ** KAPPA AT LEVEL NLAY
!    PKE           -         EDGE LEVEL PRESSURE ** KAPPA AT GROUND
!    PE            -         SURFACE PRESSURE
!    Z0            -         SURFACE ROUGHNESS
!    WATER         -         ARRAY WITH '1' OVER OCEANS
!    HS            -         DEPTH OF SURFACE LAYER
!    N             -         NUMBER OF helfsurface ITERATIONS
!    CHOOSEZ0      -         INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                                          1 - Edson Z0 for mom. and heat, high wind limit
!                                          2 - L&P Z0, high wind limit
!                                          3 - Edson Z0 for mom. only, high wind limit
!     OUTPUT:
!     -------
!    RHO           -         DENSITY AT SURFACE
!    KH            -         HEAT TRANSFER COEFFICIENT (CT*USTAR)
!    KM            -         MOMENTUM TRANSFER COEFFICIENT (CU*USTAR)
!    USTAR         -         FRICTION VELOCITY
!    XX            -         PHIM(ZETA) - DIMENSIONLESS WIND SHEAR
!    YY            -         PHIH(ZETA) - DIMENSIONLESS TEMP GRADIENT
!    CU            -         MOMENTUM TRANSPORT COEFFICIENT
!    CT            -         HEAT TRANSPORT COEFFICIENT
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer n,irun,CHOOSEZ0
      real VUS(IRUN),VVS(IRUN),VT1(IRUN),VT2(IRUN),VSH1(IRUN),VSH2(IRUN)
      real VPE(IRUN),VP(IRUN),VZ0(IRUN),LAI(IRUN),VHS(IRUN)
      integer IVWATER(IRUN)
      real VRHO(IRUN)
      real VKM(IRUN),VKH(IRUN),VUSTAR(IRUN),VXX(IRUN)
      real VYY(IRUN),VCU(IRUN),VCT(IRUN),VRIB(IRUN)
      real VZETA(IRUN),VWS(IRUN)
      real, intent(OUT) :: t2m(irun),q2m(irun),u2m(irun),v2m(irun)
      real, intent(OUT) :: t10m(irun),q10m(irun),u10m(irun),v10m(irun)
      real, intent(OUT) :: u50m(irun),v50m(irun)
      LOGICAL LWATER
      integer IVBITRIB(irun)

! Local Variables
      real VHZ(irun),VPSIM(irun),VAPSIM(irun),VPSIG(irun),VPSIHG(irun)
      real VTEMP(irun),VDZETA(irun),VDZ0(irun),VDPSIM(irun)
      real VDPSIH(irun),VZH(irun),VXX0(irun),VYY0(irun)
      real VAPSIHG(irun),VRIB1(irun)
      real VPSIH(irun),VPSIH2(irun),VH0(irun)
      real VX0PSIM(irun),VG(irun),VG0(irun),VR1MG0(irun)
      real VZ2(irun),VDZSEA(irun),VAZ0(irun),VXNUM1(irun)
      real VPSIGB2(irun),VDX(irun),VDXPSIM(irun),VDY(irun)
      real VXNUM2(irun),VDEN(irun),VAWS1(irun),VXNUM3(irun)
      real VXNUM(irun),VDZETA1(irun),VDZETA2(irun)
      real VZCOEF2(irun),VZCOEF1(irun),VTEMPLIN(irun)
      real VDPSIMC(irun),VDPSIHC(irun),VAHS(irun)
      real VTHV1(IRUN),VTHV2(IRUN),VTH1(IRUN),VTH2(IRUN),VPKE(IRUN),VPK(IRUN)

      real vz0h(irun),vh0h(irun),dummy1(irun),dummy2(irun),dummy3(irun),dummy4(irun),dummy5(irun)

! Local Variables
      real USTMX3,USTZ0S,Z0MIN,H0BYZ0,USTH0S,H0VEG,Z0VEGM,PRFAC,Z0MAX
      real XPFAC,DIFSQT
      PARAMETER ( USTMX3 =   0.0632456)
      PARAMETER ( USTZ0S =   0.2030325E-5)
      PARAMETER ( Z0MIN  =  USTZ0S/USTMX3)
      PARAMETER ( Z0MAX  =  USTZ0S/USTMX3)
      PARAMETER ( H0BYZ0 =    30.0    )
      PARAMETER ( USTH0S =  H0BYZ0*USTZ0S )
      PARAMETER ( Z0VEGM =   0.005    )
      PARAMETER ( H0VEG  =  H0BYZ0*Z0VEGM )  !! This prevents discontinuity
      PARAMETER ( PRFAC  = 0.595864   )
      PARAMETER ( XPFAC  = .55        )  
      PARAMETER ( DIFSQT  = 3.872983E-3)

      real psihdiag(irun),psimdiag(irun)
      real rvk,vk2,bmdl(irun)
      integer iwater,itype
      integer i,iter
!
      rvk = 1./MAPL_KARMAN
      vk2 = MAPL_KARMAN*MAPL_KARMAN
      DO I = 1,IRUN
      if( ivwater(i) .eq. 3 ) then 
       BMDL(i)    = 0.
!scale BMDL(i)    = (MAPL_KARMAN * XPFAC * PRFAC / DIFSQT) * exp(-lai(i)*2.)
      else
       BMDL(i)    = (MAPL_KARMAN * XPFAC * PRFAC / DIFSQT)
      endif
      enddo

!     INITIALIZATION 

      DO I = 1,IRUN
       VAHS(I) = 1. / VHS(I)
       VPKE(I) = VPE(I) ** MAPL_KAPPA
       VPK(I) = VP(I) ** MAPL_KAPPA
       VTH1(I) = VT1(I)/VPK(I)
       VTH2(I) = VT2(I)/VPKE(I)
       VTHV1(I) = VTH1(I)*( 1.0 + MAPL_VIREPS*VSH1(I))
       VTHV2(I) = VTH2(I)*( 1.0 + MAPL_VIREPS*VSH2(I))
      ENDDO

!     DETERMINE SURFACE WIND MAGNITUDE AND BULK RICHARDSON NUMBER
!
      DO I = 1,IRUN
       VWS(I) = max(VUS(I)*VUS(I) + VVS(I)*VVS(I),1.e-4)
       VRIB(I) = MAPL_CP*(VPKE(I)-VPK(I))*(VTHV1(I)-VTHV2(I)) / VWS(I)
       VWS(I) = SQRT( VWS(I) )
      ENDDO    
!
!  INITIAL GUESS FOR ROUGHNESS LENGTH Z0 OVER WATER
!
      IWATER = 0
      DO 9002 I = 1,IRUN
       IF (IVWATER(I).EQ.1)  IWATER = IWATER + 1
 9002 CONTINUE
      LWATER = .FALSE.
      IF(IWATER.GE.1)LWATER = .TRUE.
!
      IF(LWATER)THEN
       DO 9004 I = 1,IRUN
        IF (IVWATER(I).EQ.1) VZ0(I) = 0.0003
 9004  CONTINUE
      ENDIF
      do i = 1,irun
       vh0(i) = h0byz0 * vz0(i)
       if(vz0(i).ge.z0vegm)vh0(i) = h0veg
      enddo
       DO I = 1,IRUN
        VZ0H(I) = 0.001
       ENDDO

!     CU AND PSIHG FOR NEUTRALLY STRATIFIED FLOW
!
      DO 9006 I = 1,IRUN
       VHZ(I) = VHS(I) / VZ0(I)
       VPSIM(I) = LOG( VHZ(I) )
       VAPSIM(I) = 1. / VPSIM(I)
       VCU(I) = MAPL_KARMAN * VAPSIM(I)
       VUSTAR(I) = VCU(I) * VWS(I)
!
       VPSIG(I) = BMDL(i)*sqrt(max(VH0(I)*VUSTAR(I)-USTH0S,0.))
       VPSIHG(I) = VPSIM(I) + VPSIG(I)
 9006 CONTINUE
!
!     LINEAR CORRECTION FOR ERROR IN ROUGHNESS LENGTH Z0
!
      IF(LWATER)THEN
       DO 9008 I = 1,IRUN
        VTEMP(I) = 0.
 9008  CONTINUE
       CALL LINADJ(VRIB,VRIB,VWS,VWS,VZ0,VUSTAR,IVWATER,VAPSIM, &
        VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,1,.TRUE.,IRUN,VDZETA, &
        VDZ0,VDPSIM,VDPSIH,IVBITRIB, &
        VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
        VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
        VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,MAPL_KARMAN,bmdl,CHOOSEZ0)
       DO 9010 I = 1,IRUN
        IF ( IVWATER(I).EQ.1 ) THEN
         VCU(I) = VCU(I) * (1. - VDPSIM(I)*VAPSIM(I))
         VZ0(I) = VZ0(I) + VDZ0(I)
         ENDIF 
         IF ( IVWATER(I).EQ.1) THEN
         IF ( VZ0(I) .LE. Z0MIN ) VZ0(I) = Z0MIN 
         vh0(i) = h0byz0 * vz0(i)
         VPSIG(I) = VH0(I) * VCU(I) * VWS(I) - USTH0S
         if(VPSIG(I).lt.0.)  VPSIG(I) = 0.
         VPSIG(I) = SQRT( VPSIG(I) )
         VPSIG(I) = BMDL(i) * VPSIG(I)
         VPSIHG(I) = VPSIM(I) + VDPSIH(I) + VPSIG(I)
        ENDIF  
 9010  CONTINUE
!
      ENDIF
!
!  INITIAL GUESS FOR STABILITY PARAMETER ZETA
!
      DO 9012 I = 1,IRUN
       VZETA(I) = VK2 * VRIB(I) / (VCU(I) * VCU(I) * VPSIHG(I))
 9012 CONTINUE
!
!  RECOMPUTE CU, ESTIMATE PSIHG AND UPDATE ZETA AND Z0
!
      DO 9014 I = 1,IRUN
       VZH(I) = VZ0(I) * VAHS(I)
 9014 CONTINUE
      CALL PSI (VZETA,VZH,VPSIM,VTEMP,IRUN,VXX,VXX0,VYY,VYY0,2)
      DO 9016 I = 1,IRUN
       VCU(I) = MAPL_KARMAN / VPSIM(I)
       VPSIG(I) = VH0(I) * VCU(I) * VWS(I) - USTH0S
       if(VPSIG(I).lt.0.)  VPSIG(I) = 0.
       VPSIG(I) = SQRT(VPSIG(I))
       VPSIG(I) = BMDL(i) * VPSIG(I)
       VPSIHG(I) = VPSIM(I) + VPSIG(I)
       VZETA(I) = VK2 * VRIB(I) / (VCU(I) * VCU(I) * VPSIHG(I))
 9016 CONTINUE
!
!
      IF(LWATER)THEN
       DO 9018 I = 1,IRUN
        IF (IVWATER(I).EQ.1) VUSTAR(I) = VCU(I) * VWS(I)
 9018  CONTINUE
       CALL ZCSUB ( VUSTAR,VHZ,IVWATER,.FALSE.,IRUN,VTEMP,CHOOSEZ0)
       CALL ZCSUB ( VUSTAR,VHZ,IVWATER,.FALSE.,IRUN,vz0h,2)
       DO 9020 I = 1,IRUN
        IF (IVWATER(I).EQ.1 ) then
         VZ0(I) = VTEMP(I)
         IF ( VZ0(I) .LE. Z0MIN ) VZ0(I) = Z0MIN
         IF ( VZ0H(I) .LE. Z0MIN ) VZ0H(I) = Z0MIN
         vh0(i) = h0byz0 * vz0(i)
         vh0h(i) = h0byz0 * vz0h(i)
        endif
 9020  CONTINUE
      ENDIF
!
!  ITERATIVE LOOP - N ITERATIONS
!     COMPUTE CU AND CT
!
      DO 200 ITER = 1,N
       DO 9026 I = 1,IRUN
        VZH(I) = VZ0(I) * VAHS(I)
 9026  CONTINUE
       CALL PSI (VZETA,VZH,VPSIM,VPSIH,IRUN,VXX,VXX0,VYY,VYY0,1)
       DO I = 1,IRUN
        VZH(I) = VZ0H(I) * VAHS(I)
       ENDDO
       if( choosez0.eq.3 .AND. Lwater ) CALL PSI (VZETA,VZH,dummy1,VPSIH,IRUN,dummy2,dummy3,dummy4,dummy5,3)
       DO 9028 I = 1,IRUN
        VCU(I) = MAPL_KARMAN / VPSIM(I)
        VUSTAR(I) = VCU(I) * VWS(I)
!
        VPSIG(I) = VH0(I) * VUSTAR(I) - USTH0S
        if(VPSIG(I).lt.0.)  VPSIG(I) = 0.
        VPSIG(I) = SQRT(VPSIG(I))
        VPSIG(I) = BMDL(i) * VPSIG(I)
        VPSIHG(I) = VPSIH(I) + VPSIG(I)
!
!  LINEAR CORRECTIONS FOR CU, CT, ZETA, AND Z0
!
        VAPSIM(I) = VCU(I) * RVK
        VAPSIHG(I) = 1. / VPSIHG(I)
        VRIB1(I) = VAPSIM(I) * VAPSIM(I) * VPSIHG(I) * VZETA(I)
 9028  CONTINUE
!
       ITYPE = 3
       IF(ITER.EQ.N) ITYPE = 5
!
       CALL LINADJ(VRIB1,VRIB,VWS, &
        VWS,VZ0,VUSTAR,IVWATER, &
        VAPSIM,VAPSIHG,VPSIH, &
        VPSIG,VXX,VXX0, &
        VYY,VYY0,ITYPE,LWATER,IRUN,VDZETA, &
        VDZ0,VDPSIM,VDPSIH, &
        IVBITRIB, &
       VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
       VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
       VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,MAPL_KARMAN,bmdl,CHOOSEZ0)
!
!  UPDATES OF ZETA, Z0, CU AND CT
!
       DO 9032 I = 1,IRUN
        VZETA(I) = VZETA(I) * ( 1. + VDZETA(I) )
        IF (IVBITRIB(I).EQ.1 ) VZETA(I) = VPSIM(I) * VPSIM(I) * VRIB(I) * VAPSIHG(I)
 9032  CONTINUE
!
       IF ( LWATER ) THEN
        DO 9034 I = 1,IRUN
         IF (IVWATER(I).EQ.1 ) then
          VZ0(I) = VZ0(I) * ( 1. + VDZ0(I) )
          VZ0H(I) = VZ0H(I) * ( 1. + VDZ0(I) )
          IF (VZ0(I) .LE. Z0MIN ) VZ0(I) = Z0MIN
          IF (VZ0H(I) .LE. Z0MIN ) VZ0H(I) = Z0MIN
          vh0(i) = h0byz0 * vz0(i)
          vh0h(i) = h0byz0 * vz0h(i)
         endif
 9034   CONTINUE
       ENDIF
!
       IF ( ITER .EQ. N ) THEN
        DO 9036 I = 1,IRUN
         VPSIM(I) = VPSIM(I) + VDPSIM(I)
         VCU(I) = MAPL_KARMAN / VPSIM(I)
         VUSTAR(I) = VCU(I) * VWS(I)
!
         VPSIG(I) = VH0(I) * VUSTAR(I) - USTH0S
         if(VPSIG(I).lt.0.)  VPSIG(I) = 0.
         VPSIG(I) = SQRT(VPSIG(I))
         VPSIG(I) = BMDL(i) * VPSIG(I)
         VPSIHG(I) = VPSIH(I) + VDPSIH(I) + VPSIG(I)
         VCT(I) = MAPL_KARMAN / VPSIHG(I)
 9036   CONTINUE
       ENDIF

!
!  SAVE VALUES OF RIB AND WS
!
        DO 9038 I = 1,IRUN
         VRIB1(I) = VRIB(I)
 9038   CONTINUE
!
 200  CONTINUE
!
!  CALCULATE RHO-SURFACE ( KG / M**3 )
!
       DO I = 1,IRUN
        VTEMP(I) =  10. * VAHS(I) * VZETA(I)
        VZH(I) = VZ0(I) * 0.1
       ENDDO
       CALL PSI (VTEMP,VZH,VHZ,VPSIH2,IRUN,VHZ,VHZ,VHZ,VHZ,3)
       DO I = 1,IRUN
        VTEMP(I) = min(( VPSIH2(I) + VPSIG(I) ) / VPSIHG(I),1.)
        VRHO(I) = VPKE(I)*( VTH2(I) + VTEMP(I) * (VTH1(I)-VTH2(I)) )
        VRHO(I) = VPE(I)*100. / ( MAPL_RGAS * VRHO(I) )
       ENDDO
!
! interpolate uvtq to 2, 10 and 50 meters for diagnostic output
!  use psih and psim which represent non-dim change from ground
!                 to specified level
! and multiply theta by surface p**kappa to get temperatures
!
        do i = 1,irun
         vtemp(i) = 2. * vahs(i) * vzeta(i)
         vzh(i) = min(vz0(i),2.) * 0.5
        enddo
        call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)
        do i = 1,irun
         vtemp(i) = min(( psihdiag(i) + vpsig(i) ) / vpsihg(i),1.)
         t2m(i) = ( (vth2(i) + vtemp(i)* (vth1(i)-vth2(i))) ) * vpke(i)
         q2m(i) = (vsh2(i) + vtemp(i)* (vsh1(i)-vsh2(i)))
         u2m(i) = (psimdiag(i)/vpsim(i) * vus(i))
         v2m(i) = (psimdiag(i)/vpsim(i) * vvs(i))
        enddo

        do i = 1,irun
         vtemp(i) = 10. * vahs(i) * vzeta(i)
         vzh(i) = vz0(i) * 0.1
        enddo
        call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)
        do i = 1,irun
         vtemp(i) = min(( psihdiag(i) + vpsig(i) ) / vpsihg(i),1.)
         t10m(i) = ( (vth2(i) + vtemp(i)* (vth1(i)-vth2(i))) ) * vpke(i)
         q10m(i) = (vsh2(i) + vtemp(i)* (vsh1(i)-vsh2(i)))
         u10m(i) = (psimdiag(i)/vpsim(i) * vus(i))
         v10m(i) = (psimdiag(i)/vpsim(i) * vvs(i))
        enddo

        do i = 1,irun
         vtemp(i) = 50. * vahs(i) * vzeta(i)
         vzh(i) = vz0(i) * 0.02
        enddo
        call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)
        do i = 1,irun
         u50m(i) = (psimdiag(i)/vpsim(i) * vus(i))
         v50m(i) = (psimdiag(i)/vpsim(i) * vvs(i))
        enddo
!
!  EVALUATE TURBULENT TRANSFER COEFFICIENTS
!
      DO 9044 I = 1,IRUN
!!     VKH(I) = VUSTAR(I) * VCT(I)
!!     VKM(I) = VUSTAR(I) * VCU(I)
       VKH(I) = VUSTAR(I) * VCT(I) * VRHO(I)
       VKM(I) = VUSTAR(I) * VCU(I) * VRHO(I)
 9044 CONTINUE

      DO I = 1,IRUN
       VRIB(I) = MAPL_CP*(VPKE(I)-VPK(I))*(VTHV1(I)-VTHV2(I)) /    &
                max(VUS(I)*VUS(I) + VVS(I)*VVS(I),1.e-1)
      ENDDO    
!
end subroutine helfsurface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: phi
! !INTERFACE:
      SUBROUTINE PHI(Z,PHIM,PHIH,IFLAG,N)
!**********************************************************************
!
!  FUNCTION PHI - SOLVES KEYPS EQUATIONS
!               - CALLED FROM PSI
!
!  DESCRIPTION OF PARAMETERS
!     Z     -  INPUTED VALUE OF MONIN- OBUKHOV STABILITY PARAMETER ZETA
!               TIMES APPROPRIATE CONSTANT
!     PHIM  -  OUTPUTED SOLUTION OF KEYPS EQUATION FOR MOMENTUM
!     PHIH  -  OUTPUTED SOLUTION OF KEYPS EQUATION FOR SCALARS
!     IFLAG -  FLAG TO DETERMINE IF X IS NEEDED (IFLAG=2), Y IS NEEDED
!                  (IFLAG=3), OR BOTH (IFLAG=1)
!     N     -  LENGTH OF VECTOR TO BE SOLVED
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer n,iflag
      real PHIM(N),PHIH(N),Z(N)

! Local Variables
      integer I1(N),I2(N)
      real ZSTAR(N),E1(N),E2(N),TEMP1(N)
!
      real PHIM0(385),ZLINM1(75),ZLINM2(75),ZLINM3(36)
      real ZLOGM1(74),ZLOGM2(75),ZLOGM3(50)
      real PHIH0(385),ZLINH1(75),ZLINH2(75),ZLINH3(36)
      real ZLOGH1(74),ZLOGH2(75),ZLOGH3(50)
      EQUIVALENCE (PHIM0(1),ZLINM1(1)),(PHIM0(76),ZLINM2(1))
      EQUIVALENCE (PHIM0(151),ZLINM3(1))
      EQUIVALENCE (PHIM0(187),ZLOGM1(1)),(PHIM0(261),ZLOGM2(1))
      EQUIVALENCE (PHIM0(336),ZLOGM3(1))
      EQUIVALENCE (PHIH0(1),ZLINH1(1)),(PHIH0(76),ZLINH2(1))
      EQUIVALENCE (PHIH0(151),ZLINH3(1))
      EQUIVALENCE (PHIH0(187),ZLOGH1(1)),(PHIH0(261),ZLOGH2(1))
      EQUIVALENCE (PHIH0(336),ZLOGH3(1))
!
       DATA ZLOGM1/ &
                   0.697894,0.678839,0.659598,0.640260, &
        0.620910,0.601628,0.582486,0.563550,0.544877, &
        0.526519,0.508516,0.490903,0.473708,0.456951, &
        0.440649,0.424812,0.409446,0.394553,0.380133, &
        0.366182,0.352695,0.339664,0.327082,0.314938, &
        0.303222,0.291923,0.281029,0.270528,0.260409, &
        0.250659,0.241267,0.232221,0.223509,0.215119, &
        0.207041,0.199264,0.191776,0.184568,0.177628, &
        0.170949,0.164519,0.158331,0.152374,0.146641, &
        0.141123,0.135813,0.130702,0.125783,0.121048, &
        0.116492,0.112107,0.107887,0.103826,0.0999177, &
        0.0961563,0.0925364,0.0890528,0.0857003,0.0824739, &
        0.0793690,0.0763810,0.0735054,0.0707380,0.0680749, &
        0.0655120,0.0630455,0.0606720,0.0583877,0.0561895, &
        0.0540740,0.0520382,0.0500790,0.0481936,0.0463791/
       DATA ZLOGM2/ &
        0.0446330,0.0429526,0.0413355,0.0397792,0.0382816, &
        0.0368403,0.0354533,0.0341185,0.0328340,0.0315978, &
        0.0304081,0.0292633,0.0281616,0.0271013,0.0260809, &
        0.0250990,0.0241540,0.0232447,0.0223695,0.0215273, &
        0.0207168,0.0199369,0.0191862,0.0184639,0.0177687, &
        0.0170998,0.0164560,0.0158364,0.0152402,0.0146664, &
        0.0141142,0.0135828,0.0130714,0.0125793,0.0121057, &
        0.0116499,0.0112113,0.0107892,0.0103830,0.999210E-2, &
        0.961590E-2,0.925387E-2,0.890547E-2,0.857018E-2,0.824752E-2, &
        0.793701E-2,0.763818E-2,0.735061E-2,0.707386E-2,0.680754E-2, &
        0.655124E-2,0.630459E-2,0.606722E-2,0.583880E-2,0.561897E-2, &
        0.540742E-2,0.520383E-2,0.500791E-2,0.481937E-2,0.463792E-2, &
        0.446331E-2,0.429527E-2,0.413355E-2,0.397793E-2,0.382816E-2, &
        0.368403E-2,0.354533E-2,0.341185E-2,0.328340E-2,0.315978E-2, &
        0.304082E-2,0.292633E-2,0.281616E-2,0.271013E-2,0.260809E-2/
       DATA ZLOGM3/ &
        0.250990E-2,0.241541E-2,0.232447E-2,0.223695E-2,0.215273E-2, &
        0.207168E-2,0.199369E-2,0.191862E-2,0.184639E-2,0.177687E-2, &
        0.170998E-2,0.164560E-2,0.158364E-2,0.152402E-2,0.146664E-2, &
        0.141142E-2,0.135828E-2,0.130714E-2,0.125793E-2,0.121057E-2, &
        0.116499E-2,0.112113E-2,0.107892E-2,0.103830E-2,0.999210E-3, &
        0.961590E-3,0.925387E-3,0.890547E-3,0.857018E-3,0.824752E-3, &
        0.793701E-3,0.763818E-3,0.735061E-3,0.707386E-3,0.680754E-3, &
        0.655124E-3,0.630459E-3,0.606722E-3,0.583880E-3,0.561897E-3, &
        0.540742E-3,0.520383E-3,0.500791E-3,0.481937E-3,0.463792E-3, &
        0.446331E-3,0.429527E-3,0.413355E-3,0.397793E-3,0.382816E-3/
       DATA ZLOGH1/ &
                   0.640529,0.623728,0.606937,0.590199, &
        0.573552,0.557032,0.540672,0.524504,0.508553, &
        0.492843,0.477397,0.462232,0.447365,0.432809, &
        0.418574,0.404670,0.391103,0.377878,0.364999, &
        0.352468,0.340284,0.328447,0.316954,0.305804, &
        0.294992,0.284514,0.274364,0.264538,0.255028, &
        0.245829,0.236933,0.228335,0.220026,0.211999, &
        0.204247,0.196762,0.189537,0.182564,0.175837, &
        0.169347,0.163088,0.157051,0.151231,0.145620, &
        0.140211,0.134998,0.129974,0.125133,0.120469, &
        0.115975,0.111645,0.107475,0.103458,0.995895E-1, &
        0.958635E-1,0.922753E-1,0.888199E-1,0.854925E-1,0.822886E-1, &
        0.792037E-1,0.762336E-1,0.733739E-1,0.706208E-1,0.679704E-1, &
        0.654188E-1,0.629625E-1,0.605979E-1,0.583217E-1,0.561306E-1, &
        0.540215E-1,0.519914E-1,0.500373E-1,0.481564E-1,0.463460E-1/
       DATA ZLOGH2/ &
        0.446034E-1,0.429263E-1,0.413120E-1,0.397583E-1,0.382629E-1, &
        0.368237E-1,0.354385E-1,0.341053E-1,0.328222E-1,0.315873E-1, &
        0.303988E-1,0.292550E-1,0.281541E-1,0.270947E-1,0.260750E-1, &
        0.250937E-1,0.241494E-1,0.232405E-1,0.223658E-1,0.215240E-1, &
        0.207139E-1,0.199342E-1,0.191839E-1,0.184618E-1,0.177669E-1, &
        0.170981E-1,0.164545E-1,0.158351E-1,0.152390E-1,0.146653E-1, &
        0.141133E-1,0.135820E-1,0.130707E-1,0.125786E-1,0.121051E-1, &
        0.116494E-1,0.112108E-1,0.107888E-1,0.103826E-1,0.999177E-2, &
        0.961561E-2,0.925360E-2,0.890523E-2,0.856997E-2,0.824733E-2, &
        0.793684E-2,0.763803E-2,0.735048E-2,0.707375E-2,0.680743E-2, &
        0.655114E-2,0.630450E-2,0.606715E-2,0.583873E-2,0.561891E-2, &
        0.540737E-2,0.520379E-2,0.500787E-2,0.481933E-2,0.463789E-2, &
        0.446328E-2,0.429524E-2,0.413353E-2,0.397790E-2,0.382814E-2, &
        0.368401E-2,0.354532E-2,0.341184E-2,0.328338E-2,0.315977E-2, &
        0.304081E-2,0.292632E-2,0.281615E-2,0.271012E-2,0.260809E-2/
       DATA ZLOGH3/ &
        0.250990E-2,0.241540E-2,0.232446E-2,0.223695E-2,0.215273E-2, &
        0.207168E-2,0.199368E-2,0.191862E-2,0.184639E-2,0.177687E-2, &
        0.170997E-2,0.164559E-2,0.158364E-2,0.152402E-2,0.146664E-2, &
        0.141142E-2,0.135828E-2,0.130714E-2,0.125793E-2,0.121057E-2, &
        0.116499E-2,0.112113E-2,0.107892E-2,0.103830E-2,0.999209E-3, &
        0.961590E-3,0.925387E-3,0.890546E-3,0.857018E-3,0.824752E-3, &
        0.793700E-3,0.763818E-3,0.735061E-3,0.707386E-3,0.680754E-3, &
        0.655124E-3,0.630459E-3,0.606722E-3,0.583880E-3,0.561897E-3, &
        0.540742E-3,0.520383E-3,0.500791E-3,0.481937E-3,0.463792E-3, &
        0.446331E-3,0.429527E-3,0.413355E-3,0.397793E-3,0.382816E-3/
 
       DATA ZLINM1/ &
        0.964508,0.962277,0.960062,0.957863,0.955680, &
        0.953512,0.951359,0.949222,0.947100,0.944992, &
        0.942899,0.940821,0.938758,0.936709,0.934673, &
        0.932652,0.930645,0.928652,0.926672,0.924706, &
        0.922753,0.920813,0.918886,0.916973,0.915072, &
        0.913184,0.911308,0.909445,0.907594,0.905756, &
        0.903930,0.902115,0.900313,0.898522,0.896743, &
        0.894975,0.893219,0.891475,0.889741,0.888019, &
        0.886307,0.884607,0.882917,0.881238,0.879569, &
        0.877911,0.876264,0.874626,0.872999,0.871382, &
        0.869775,0.868178,0.866591,0.865013,0.863445, &
        0.861887,0.860338,0.858798,0.857268,0.855747, &
        0.854235,0.852732,0.851238,0.849753,0.848277, &
        0.846809,0.845350,0.843900,0.842458,0.841025, &
        0.839599,0.838182,0.836774,0.835373,0.833980/
       DATA ZLINM2/ &
        0.832596,0.831219,0.829850,0.828489,0.827136, &
        0.825790,0.824451,0.823121,0.821797,0.820481, &
        0.819173,0.817871,0.816577,0.815289,0.814009, &
        0.812736,0.811470,0.810210,0.808958,0.807712, &
        0.806473,0.805240,0.804015,0.802795,0.801582, &
        0.800376,0.799176,0.797982,0.796794,0.795613, &
        0.794438,0.793269,0.792106,0.790949,0.789798, &
        0.788652,0.787513,0.786380,0.785252,0.784130, &
        0.783014,0.781903,0.780798,0.779698,0.778604, &
        0.777516,0.776432,0.775354,0.774282,0.773215, &
        0.772153,0.771096,0.770044,0.768998,0.767956, &
        0.766920,0.765888,0.764862,0.763840,0.762824, &
        0.761812,0.760805,0.759803,0.758805,0.757813, &
        0.756824,0.755841,0.754862,0.753888,0.752918, &
        0.751953,0.750992,0.750035,0.749083,0.748136/
       DATA ZLINM3/ &
        0.747192,0.746253,0.745318,0.744388,0.743462, &
        0.742539,0.741621,0.740707,0.739798,0.738892, &
        0.737990,0.737092,0.736198,0.735308,0.734423, &
        0.733540,0.732662,0.731788,0.730917,0.730050, &
        0.729187,0.728328,0.727472,0.726620,0.725772, &
        0.724927,0.724086,0.723248,0.722414,0.721584, &
        0.720757,0.719933,0.719113,0.718296,0.717483, &
        0.716673/
       DATA ZLINH1/ &
        0.936397,0.932809,0.929287,0.925827,0.922429, &
        0.919089,0.915806,0.912579,0.909405,0.906284, &
        0.903212,0.900189,0.897214,0.894284,0.891399, &
        0.888558,0.885759,0.883001,0.880283,0.877603, &
        0.874962,0.872357,0.869788,0.867255,0.864755, &
        0.862288,0.859854,0.857452,0.855081,0.852739, &
        0.850427,0.848144,0.845889,0.843662,0.841461, &
        0.839287,0.837138,0.835014,0.832915,0.830841, &
        0.828789,0.826761,0.824755,0.822772,0.820810, &
        0.818869,0.816949,0.815050,0.813170,0.811310, &
        0.809470,0.807648,0.805845,0.804060,0.802293, &
        0.800543,0.798811,0.797095,0.795396,0.793714, &
        0.792047,0.790396,0.788761,0.787141,0.785535, &
        0.783945,0.782369,0.780807,0.779259,0.777724, &
        0.776204,0.774696,0.773202,0.771720,0.770251/
       DATA ZLINH2/ &
        0.768795,0.767351,0.765919,0.764499,0.763091, &
        0.761694,0.760309,0.758935,0.757571,0.756219, &
        0.754878,0.753547,0.752226,0.750916,0.749616, &
        0.748326,0.747045,0.745775,0.744514,0.743262, &
        0.742020,0.740787,0.739563,0.738348,0.737141, &
        0.735944,0.734755,0.733574,0.732402,0.731238, &
        0.730083,0.728935,0.727795,0.726664,0.725539, &
        0.724423,0.723314,0.722213,0.721119,0.720032, &
        0.718952,0.717880,0.716815,0.715756,0.714704, &
        0.713660,0.712621,0.711590,0.710565,0.709547, &
        0.708534,0.707529,0.706529,0.705536,0.704549, &
        0.703567,0.702592,0.701623,0.700660,0.699702, &
        0.698750,0.697804,0.696863,0.695928,0.694998, &
        0.694074,0.693155,0.692241,0.691333,0.690430, &
        0.689532,0.688639,0.687751,0.686868,0.685990/
       DATA ZLINH3/ &
        0.685117,0.684249,0.683386,0.682527,0.681673, &
        0.680824,0.679979,0.679139,0.678303,0.677472, &
        0.676645,0.675823,0.675005,0.674191,0.673381, &
        0.672576,0.671775,0.670978,0.670185,0.669396, &
        0.668611,0.667830,0.667054,0.666281,0.665512, &
        0.664746,0.663985,0.663227,0.662473,0.661723, &
        0.660977,0.660234,0.659495,0.658759,0.658027, &
        0.657298/

        integer i
!
      DO 9002 I = 1,N
       ZSTAR(I)    = 100. * Z(I) - 14.
 9002 CONTINUE
!
      DO 9004 I = 1,N
       TEMP1(I) = Z(I)*0.5
       IF( Z(I) .LE. 2. )TEMP1(I) = 1.
       TEMP1(I) = LOG10(TEMP1(I))
       TEMP1(I) = (TEMP1(I) + 9.3) * 20.
       IF( Z(I) .GT. 2. ) ZSTAR(I) = TEMP1(I)
       IF( Z(I).GT.1.78e10 ) ZSTAR(I) = 384.9999
 9004  CONTINUE
!
 60    CONTINUE
!
      DO 9006 I = 1,N
       I1(I) = ZSTAR(I)
       I2(I) = I1(I) + 1
       TEMP1(I) = ZSTAR(I) - I1(I)
!
 9006  CONTINUE
!
      IF( IFLAG .GT. 2 ) GO TO 100
       DO 9008 I = 1,N
       if( z(i).ge.0.15 ) then
       E1(I) = PHIM0( I1(I) )
       E2(I) = PHIM0( I2(I) )
       PHIM(I)  = TEMP1(I) * ( E2(I)-E1(I) )
       PHIM(I)  = PHIM(I) +   E1(I)
       endif
 9008  CONTINUE

  100 CONTINUE
!
      IF( IFLAG .EQ. 2 ) GO TO 200
       DO 9010 I = 1,N
       if( z(i).ge.0.15 ) then
       E1(I) = PHIH0( I1(I) )
       E2(I) = PHIH0( I2(I) )
       PHIH(I)  = TEMP1(I) * ( E2(I)-E1(I) )
       PHIH(I)  = PHIH(I) +   E1(I)
       endif
 9010  CONTINUE

  200 CONTINUE
!
       DO 9012 I = 1,N
       ZSTAR(I) = -Z(I)
 9012  CONTINUE
!
      IF( IFLAG .GT. 2 ) GO TO 300
       DO 9014 I = 1,N
       IF( Z(I) .LT. 0.15 ) PHIM(I) = 1. + ZSTAR(I) &
           *(0.25+ZSTAR(I)*(0.09375+ZSTAR(I)* &
           (0.03125+0.00732422 * ZSTAR(I))))
 9014  CONTINUE
!
  300 CONTINUE
      IF( IFLAG .EQ. 2 ) GO TO 500
       DO 9016 I = 1,N
       IF( Z(I) .LT. 0.15 ) THEN
       PHIH(I) =1.+ Z(I) * (0.5+ZSTAR(I)*(0.375+ZSTAR(I)* &
           (0.5+ZSTAR(I)*(0.8203125+ZSTAR(I)* &
           (1.5+2.93262*ZSTAR(I))))))
       PHIH(I) = 1. / PHIH(I)
      ENDIF
 9016  CONTINUE
!
  500 CONTINUE

end subroutine phi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: psi
! !INTERFACE:
      SUBROUTINE PSI(VZZ,VZH,VPSIM,VPSIH,IRUN,VX,VXS,VY,VYS,IFLAG)
!**********************************************************************
!
!  SUBROUTINE PSI - DETERMINES DIMENSIONLESS WIND AND
!                    SCALAR PROFILES IN SURFACE LAYER
!                 - CALLED FROM helfsurface
!
!  DESCRIPTION OF PARAMETERS
!     ZZ   -  INPUTED VALUE OF MONIN- OBUKHOV STABILITY PARAMETER ZETA
!     ZH   -  INPUTED VALUE OF Z0 DIVIDED BY SFC LAYER HEIGHT
!     PSIM -  OUTPUTED VALUE OF DIMENSIONLESS WIND
!     PSIH -  OUTPUTED VALUE OF DIMENSIONLESS SCALAR
!     X    -  OUTPUTED VALUE OF PHIM(ZETA)
!     XS   -  OUTPUTED VALUE OF PHIM(ZETA0)
!     Y    -  OUTPUTED VALUE OF PHIH(ZETA)
!     YS   -  OUTPUTED VALUE OF PHIH(ZETA0)
!     IFLAG-  FLAG TO DETERMINE IF CU IS NEEDED (IFLAG=2),
!                  IF CT IS NEEDED (IFLAG=3), OR BOTH (IFLAG=1)
!  SUBPROGRAMS NEEDED
!     PHI  -  COMPUTES SIMILARITY FUNCTION FOR MOMENTUM AND SCALARS
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer irun,iflag
      real VZZ(IRUN),VZH(IRUN),VPSIM(IRUN),VPSIH(IRUN), &
           VX(IRUN),VXS(IRUN),VY(IRUN),VYS(IRUN)
 
! Local Variables
      real ZWM,RZWM,Z0M,ZCM,RZCM,CM1,CM2,CM6,CM7,CM8ARG,YCM
      PARAMETER ( ZWM     =    1.    )
      PARAMETER ( RZWM    =  1./ZWM  )
      PARAMETER ( Z0M     =    0.2    )
      PARAMETER ( ZCM     =    42.    )
      PARAMETER ( RZCM    =  1./ZCM  )
      PARAMETER ( CM1     =  1./126. )
      PARAMETER ( CM2     =  1./(6.*CM1)  )
      PARAMETER ( CM6     =  6. / ( 1. + 6.*CM1 )  )
      PARAMETER ( CM7     =  CM2 + ZWM  )
      PARAMETER ( CM8ARG  =  CM7*ZCM*RZWM / (CM2+ZCM)  )
      PARAMETER ( YCM     =  6. / ( 1. + 6.*CM1*ZCM )  )

      integer INTSTB(irun),INTZ0(irun)
      real ZZ0(irun),Z(irun),Z2(irun),Z1(irun),Z0(irun)
      real X0(irun),X1(irun),Y0(irun),Y1(irun)
      real PSI2(irun),TEMP(irun)
      real HZ(irun),ARG0(irun),ARG1(irun),DX(irun)
      real X0NUM(irun),X1NUM(irun),X0DEN(irun)
      real X1DEN(irun),Y1DEN(irun),Z2ZWM(irun)
      real cm3,cm4,cm5,cm8
      integer ibit,indx
      integer i
!
      CM3 =   sqrt( 0.2/CM1-0.01 )
      CM4 =   1./CM3
      CM5 =  (10.-CM1) / (10.*CM1*CM3)
      CM8 =   6. * LOG(CM8ARG)
!
      DO 9000 I = 1,IRUN
       VPSIM(I) = 0.
       VPSIH(I) = 0.
       VX(I) = 0.
       VXS(I) = 0.
       VY(I) = 0.
       VYS(I) = 0.
       ZZ0(I) = VZH(I)*VZZ(I)
 9000 CONTINUE
      IBIT = 0
      DO 9122 I = 1,IRUN
       IF(VZZ(I).LE.-1.e-7)IBIT = IBIT + 1
 9122 CONTINUE
      DO 9022 I = 1,IRUN
       IF(VZZ(I).LE.-1.e-7)THEN
        INTSTB(I) = 1
       ELSE
        INTSTB(I) = 0
       ENDIF
 9022 CONTINUE
!
! ****************************************
! *****    UNSTABLE SURFACE LAYER    *****
! ****************************************
!
      IF(IBIT.LE.0)  GO TO 100
!
      indx = 0
      DO 9002 I = 1,IRUN
       IF (INTSTB(I).EQ.1)THEN
        indx = indx + 1
        Z(indx) = VZZ(I)
        Z0(indx) = ZZ0(I)
       ENDIF
 9002 CONTINUE
!
      DO 9004 I = 1,IBIT
       Z(I) = -18. * Z(I)
       Z0(I) = -18. * Z0(I)
 9004 CONTINUE
 
      CALL PHI( Z,X1,Y1,IFLAG,IBIT )
      CALL PHI( Z0,X0,Y0,IFLAG,IBIT )
 
! ****************************
! *****    COMPUTE PSIM  *****
! ****************************
!
      IF(IFLAG.GE.3) GO TO 75
!
      DO 9006 I = 1,IBIT
       ARG1(I) = 1. - X1(I)
       IF ( Z(I) .LT. 0.013 ) ARG1(I) = Z(I) * ( 0.25 -  0.09375 * Z(I) )
!
       ARG0(I)  = 1. - X0(I)
       IF ( Z0(I) .LT. 0.013 ) ARG0(I) = Z0(I) * ( 0.25 -  0.09375 * Z0(I) )
!
       ARG1(I) = ARG1(I) * ( 1.+X0(I) )
       ARG0(I) = ARG0(I) * ( 1.+X1(I) )
       DX(I) = X1(I) - X0(I)
       ARG1(I) = ARG1(I) / ARG0(I)
       ARG0(I) = -DX(I) / ( 1. + X1(I)*X0(I) )
       ARG0(I) = ATAN( ARG0(I) )
       ARG1(I) = LOG( ARG1(I) )
       PSI2(I) = 2. * ARG0(I) + ARG1(I)
       PSI2(I) = PSI2(I) + DX(I)
 9006 CONTINUE
!
      indx = 0
      DO 9008 I = 1,IRUN
       IF( INTSTB(I).EQ.1 ) THEN
        indx = indx + 1
        VPSIM(I) = PSI2(indx)
        VX(I) = X1(indx)
        VXS(I) = X0(indx)
       ENDIF
 9008 CONTINUE
!
! ****************************
! *****    COMPUTE PSIH  *****
! ****************************
!
      IF(IFLAG.EQ.2) GO TO 100
!
  75  CONTINUE
      DO 9010 I = 1,IBIT
       ARG1(I) = 1. - Y1(I)
       IF( Z(I) .LT. 0.0065 ) ARG1(I) = Z(I) * ( 0.5 -  0.625 * Z(I) )
!
       ARG0(I)  = 1. - Y0(I)
       IF( Z0(I) .LT. 0.0065 ) ARG0(I) = Z0(I) * ( 0.5 -  0.625 * Z0(I) )
!
       ARG1(I) = ARG1(I) * ( 1. + Y0(I) )
       ARG0(I) = ARG0(I) * ( 1. + Y1(I) )
       ARG1(I) = ARG1(I) / ARG0(I)
       PSI2(I) = LOG( ARG1(I) )
       PSI2(I) = PSI2(I) - Y1(I) + Y0(I)
 9010 CONTINUE
!
      indx = 0
      DO 9012 I = 1,IRUN
       IF( INTSTB(I).EQ.1 ) THEN
       indx = indx + 1
       VPSIH(I) = PSI2(indx)
       VY(I) = Y1(indx)
       VYS(I) = Y0(indx)
       ENDIF
 9012 CONTINUE
!
! **************************************
! *****    STABLE SURFACE LAYER    *****
! **************************************
!
  100 CONTINUE
      IBIT = 0
      DO 9114 I = 1,IRUN
       IF(VZZ(I).GT.-1.e-7)THEN
        IBIT = IBIT + 1
       ENDIF
 9114 CONTINUE
      DO 9014 I = 1,IRUN
       IF(VZZ(I).GT.-1.e-7)THEN
        INTSTB(I) = 1
       ELSE
        INTSTB(I) = 0
       ENDIF
 9014 CONTINUE
      IF(IBIT.LE.0)  GO TO 300
      indx = 0
      DO 9016 I = 1,IRUN
       IF (INTSTB(I).EQ.1)THEN
        indx = indx + 1
        Z(indx) = VZZ(I)
        Z0(indx) = ZZ0(I)
        ARG1(indx) = VZH(I)
       ENDIF
 9016 CONTINUE

      DO 9018 I = 1,IBIT
       HZ(I) = 1. / ARG1(I)
       Z1(I) = Z(I)
       Z2(I) = ZWM
!
       IF ( Z(I) .GT. ZWM ) THEN
        Z1(I) = ZWM
        Z2(I) = Z(I)
       ENDIF
!
       IF ( Z0(I) .GT. Z0M ) THEN
        Z0(I) = Z0M
        INTZ0(I) = 1
       ELSE
        INTZ0(I) = 0
       ENDIF
!
       X1NUM(I) = 1. + 5. * Z1(I)
       X0NUM(I) = 1. + 5. * Z0(I)
       X1DEN(I) = 1. / (1. + CM1 * (X1NUM(I) * Z1(I)) )
       X0DEN(I) = 1. + CM1 * (X0NUM(I) * Z0(I))
!
       IF ( (INTZ0(I).EQ.1) .OR. (Z(I).GT.ZWM) ) &
            HZ(I) = Z1(I) / Z0(I)
       ARG1(I) = HZ(I)*HZ(I)*X0DEN(I)*X1DEN(I)
       ARG1(I) = LOG( ARG1(I) )
       ARG1(I) = 0.5 * ARG1(I)
       ARG0(I) = (Z1(I) + 0.1) * (Z0(I) + 0.1)
       ARG0(I) = CM3 + ARG0(I) * CM4
       ARG0(I) = ( Z1(I) - Z0(I) ) / ARG0(I)
       ARG0(I) = ATAN( ARG0(I) )
       TEMP(I) = ARG1(I) + CM5 * ARG0(I)
!
       X0(I) = X0NUM(I) / X0DEN(I)
       IF ( INTZ0(I).EQ.1 ) X0(I) = 0.
       Z2ZWM(I) = Z2(I) * RZWM
 9018 CONTINUE
!
! ****************************
! *****    COMPUTE PSIM  *****
! ****************************
!
      IF( IFLAG.GE.3 ) GO TO 225
!
      DO 9020 I = 1,IBIT
       X1(I) = X1NUM(I) * X1DEN(I)
       ARG1(I) = LOG( Z2ZWM(I) )
       PSI2(I) = TEMP(I) + CM6 * ARG1(I)
 9020 CONTINUE
!
      indx = 0
      DO 9030 I = 1,IRUN
       IF( INTSTB(I).EQ.1 ) THEN
       indx = indx + 1
       VPSIM(I) = PSI2(indx)
       VX(I) = X1(indx)
       VXS(I) = X0(indx)
       ENDIF
 9030 CONTINUE
!
! ****************************
! *****    COMPUTE PSIH  *****
! ****************************
!
       IF(IFLAG.EQ.2)GO TO 300
!
  225 CONTINUE
      DO 9024 I = 1,IBIT
       Y1DEN(I) = 1. + CM1 * ( X1NUM(I) * Z(I) )
       Y1(I) = X1NUM(I) / Y1DEN(I)
       ARG1(I) = CM7 * Z2ZWM(I) / ( CM2 + Z2(I) )
       ARG0(I) = 6.
       IF ( Z2(I) .GT. ZCM ) THEN
        Y1(I) = YCM
        ARG1(I) = Z2(I) * RZCM
        ARG0(I) = YCM
        TEMP(I) = TEMP(I) + CM8
       ENDIF
       ARG1(I) = LOG( ARG1(I) )
       PSI2(I) = TEMP(I) + ARG0(I) * ARG1(I)
 9024 CONTINUE
!
      indx = 0
      DO 9026 I = 1,IRUN
       IF( INTSTB(I).EQ.1 ) THEN
       indx = indx + 1
       VPSIH(I) = PSI2(indx)
       VY(I) = Y1(indx)
       VYS(I) = X0(indx)
       ENDIF
 9026 CONTINUE
!
  300 CONTINUE
!
end subroutine psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: linadj
! !INTERFACE:
      SUBROUTINE LINADJ ( VRIB1,VRIB2,VWS1,VWS2,VZ1,VUSTAR,IWATER, &
       VAPSIM, VAPSIHG,VPSIH,VPSIG,VX,VX0,VY,VY0,ITYPE,LWATER,IRUN, &
       VDZETA,VDZ0,VDPSIM,VDPSIH,INTRIB, &
       VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
       VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
       VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,vk,bmdl,CHOOSEZ0)
!
!**********************************************************************
!
!  ARGUMENTS ::
!
!     INPUT:
!     ------
!    RIB1          -         BULK RICHARDSON NUMBER OF INPUT STATE
!    RIB2          -         DESIRED BULK RICH NUMBER OF OUTPUT STATE
!    WS1           -         SURFACE WIND SPEED OF INPUT STATE
!    WS2           -         DESIRED SURFACE WIND SPEED OF OUTPUT STATE
!    Z1            -         INPUT VALUE OF ROUGHNESS HEIGHT
!    USTAR         -         INPUT VALUE OF CU * WS
!    WATER         -         BIT ARRAY - '1' WHERE OCEAN
!    APSIM         -         (1/PSIM)
!    APSIHG        -         ( 1 / (PSIH+PSIG) )
!    PSIH          -         NON-DIM TEMP GRADIENT
!    PSIG          -         PSIH FOR THE MOLECULAR LAYER
!    X             -         PHIM(ZETA) - DERIVATIVE OF PSIM
!    X0            -         PHIM(ZETA0)
!    Y             -         PHIH(ZETA) - DERIVATIVE OF PSIH
!    Y0            -         PHIH(ZETA0)
!    ITYPE         -         INTEGER FLAG :
!                               1    = NEUTRAL ADJUSTMENT
!                               3, 5 = ADJUSTMENT INSIDE LOOP
!                               5    = ADJUST CU AND CT
!    LWATER        -         LOGICAL - .TRUE. IF THERE ARE WATER POINTS
!    CHOOSEZ0      -         INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                                          1 - Edson Z0 for mom. and heat, high wind limit
!                                          2 - L&P Z0, high wind limit
!                                          3 - Edson Z0 for mom. only, high wind limit
!
!     OUTPUT:
!     -------
!    DZETA         -         D LOG ZETA
!    DZ0           -         D Z0 (ITYPE 1) OR D LOG Z0 (ITYPE 2-5)
!    DPSIM         -         D PSIM
!    DPSIH         -         D PSIH
!    BITRIB        -         BIT ARRAY - '1' WHERE RIB1 = 0
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer irun,itype,CHOOSEZ0
      real VRIB1(IRUN),VRIB2(IRUN)
      real VWS1(IRUN),VWS2(IRUN),VZ1(IRUN),VUSTAR(IRUN)
      integer IWATER(IRUN)
      real VAPSIM(IRUN),VAPSIHG(IRUN)
      real VPSIH(IRUN),VPSIG(IRUN),VX(IRUN)
      real VX0(IRUN),VY(IRUN),VY0(IRUN)
      LOGICAL LWATER
      real VDZETA(IRUN),VDZ0(IRUN),VDPSIM(IRUN)
      real VDPSIH(IRUN)
      integer INTRIB(IRUN)
      real VX0PSIM(irun),VG(irun),VG0(irun),VR1MG0(irun)
      real VZ2(irun),VDZSEA(irun),VAZ0(irun),VXNUM1(irun)
      real VPSIGB2(irun),VDX(irun),VDXPSIM(irun),VDY(irun)
      real VXNUM2(irun),VDEN(irun),VAWS1(irun),VXNUM3(irun)
      real VXNUM(irun),VDZETA1(irun),VDZETA2(irun)
      real VZCOEF2(irun),VZCOEF1(irun),VTEMPLIN(irun)
      real VDPSIMC(irun),VDPSIHC(irun),bmdl(irun)

! Local Variables
      real xx0max,prfac,xpfac,difsqt,ustz0s,h0byz0,usth0s
      PARAMETER ( XX0MAX  =   1.49821 )
      PARAMETER ( PRFAC  = 0.595864   )
      PARAMETER ( XPFAC  = .55        )  
      PARAMETER ( DIFSQT  = 3.872983E-3)
      PARAMETER ( USTZ0S =   0.2030325E-5)
      PARAMETER ( H0BYZ0 =    30.0    )
      PARAMETER ( USTH0S =  H0BYZ0*USTZ0S )

      integer VINT1(irun),VINT2(irun)
      real vk,b2uhs(irun)
      integer i
!
      do i = 1,irun
      B2UHS(i)   = BMDL(i) * BMDL(i) * USTH0S
      enddo

!   COMPUTE X0/PSIM, 1/Z0, G, G0, 1/(1-G0),
!     DEL LOG Z0, D LOG ZO / D USTAR
!
      IF ( (ITYPE.EQ.1) .AND. LWATER ) THEN
       DO 9000 I = 1,IRUN
        IF (IWATER(I).EQ.1) VX0PSIM(I) = VAPSIM(I)
 9000  CONTINUE
      ENDIF
      IF ( ITYPE .GE. 3 ) THEN
       DO 9002 I = 1,IRUN
        VX0PSIM(I) = VX0(I) * VAPSIM(I)
 9002  CONTINUE
      ENDIF
!
       DO 9004 I = 1,IRUN
        VDZ0(I) = 0.
        VG(I) = 0.
        VG0(I) = 0.
        VR1MG0(I) = 1.
 9004  CONTINUE
!
       IF ( LWATER ) THEN
        CALL ZCSUB ( VUSTAR,VDZSEA,IWATER,.TRUE.,IRUN,VZ2,CHOOSEZ0)

        VDZSEA = min( VDZSEA, 0.2*VZ1/VAPSIM ) ! To prevent Divide by Zero as VG0 => 1.0
!
        DO 9006 I = 1,IRUN
         IF ( IWATER(I).EQ.1) THEN
          VAZ0(I) = 1. / VZ1(I)
          VG(I) = VDZSEA(I) * VAZ0(I)
          VG0(I) = VX0PSIM(I) * VG(I)
          VR1MG0(I) = 1. / ( 1. - VG0(I) )
          VDZ0(I) = ( VZ2(I) - VZ1(I) ) * VR1MG0(I)
         ENDIF
 9006   CONTINUE
       ENDIF
!
      IF ( LWATER .AND. (ITYPE.GE.3) ) THEN
       DO 9008 I = 1,IRUN
        IF (IWATER(I).EQ.1) VDZ0(I) = VDZ0(I) * VAZ0(I)
 9008  CONTINUE
      ENDIF
!
!   COMPUTE NUM1,NUM2,NUM3, DEN
!
      IF (ITYPE.GE.3) THEN
       DO 9010 I = 1,IRUN
        VXNUM1(I) = 0.
        IF (VRIB1(I).EQ.0.) THEN
         INTRIB(I) = 1
        ELSE
         INTRIB(I) = 0
        ENDIF
        IF ( INTRIB(I).EQ.0 ) VXNUM1(I) = 1. / VRIB1(I)
        VPSIGB2(I) = 0.
        if(vpsig(i).gt.0.)VPSIGB2(I) = &
              0.5 * ( vpsig(i)*vpsig(i) + b2uhs(i) ) / vpsig(i)
        VDX(I) = VX(I) - VX0(I)
        VDXPSIM(I) = VDX(I) * VAPSIM(I)
        VDY(I) = VY(I) - VY0(I)
        VXNUM3(I) = - VPSIGB2(I)
!
        IF ( LWATER ) THEN
         IF (IWATER(I).EQ.1) THEN
          VDXPSIM(I) = VDXPSIM(I) * VR1MG0(I)
          VXNUM3(I) = VXNUM3(I) + VG(I) * ( VY0(I) - VPSIGB2(I) )
          VXNUM2(I) = VY0(I) - VPSIGB2(I) - VX0PSIM(I) * VPSIGB2(I)
          VXNUM2(I) = (VXNUM2(I) * VAPSIHG(I)) - 2. * VX0PSIM(I)
          VXNUM2(I) = VXNUM2(I) * VDZ0(I)
         ENDIF
        ENDIF
!
        VDEN(I) = VDY(I) + VDXPSIM(I) * VXNUM3(I)
        VDEN(I) = ( 1. + VDEN(I) * VAPSIHG(I) ) - 2. * VDXPSIM(I)
 9010  CONTINUE
      ENDIF
!
      IF (ITYPE.EQ.5) THEN
       DO 9012 I = 1,IRUN
        VAWS1(I) = VR1MG0(I) / VWS1(I)
        VXNUM3(I) = VXNUM3(I) * VAPSIHG(I)
!
        IF ( LWATER ) THEN
         IF(IWATER(I).EQ.1) THEN
          VXNUM3(I) = VXNUM3(I) - 2. * VG0(I)
          VXNUM3(I) = VAWS1(I) * VXNUM3(I)
         ENDIF
        ENDIF
 9012  CONTINUE
      ENDIF
!
!   COMPUTE D LOG ZETA
!
      IF (ITYPE.GE.3) THEN
       DO 9014 I = 1,IRUN
        VXNUM(I) = VRIB2(I) - VRIB1(I)
        IF( (VX0(I).GT.XX0MAX).AND.(VXNUM(I).GE.0.) )VXNUM(I) = 0.
        VXNUM(I) = VXNUM1(I) * VXNUM(I)
 9014  CONTINUE
!
       DO 9018 I = 1,IRUN
        VDZETA1(I) = VXNUM(I)
        IF(LWATER.AND.(IWATER(I).EQ.1)) VXNUM(I) = VXNUM(I) + VXNUM2(I)
        IF ( VDEN(I) .LT.0.1 ) VDEN(I) = 0.1
 9018  CONTINUE
!
       DO 9020 I = 1,IRUN
        VDZETA(I) = VXNUM(I) / VDEN(I)
 9020  CONTINUE
       DO 9022 I = 1,IRUN
        IF((VRIB2(I).EQ.0.).OR.(VDZETA(I).LE.-1.))VDZETA(I) = VDZETA1(I)
 9022  CONTINUE
      ENDIF
!
!   COMPUTE D LOG Z0
!
      IF ( LWATER .AND. (ITYPE.GE.3) )THEN
       DO 9026 I = 1,IRUN
        IF( IWATER(I).EQ.1 ) THEN
         VZCOEF2(I) = VG(I) * VDXPSIM(I)
         VDZ0(I) = VDZ0(I) - VZCOEF2(I) * VDZETA(I)
        ENDIF
 9026  CONTINUE
      ENDIF
!
      IF ( LWATER .AND. (ITYPE.EQ.5) ) THEN
       DO 9028 I = 1,IRUN
        IF(IWATER(I).EQ.1) VZCOEF1(I) = VG(I) * VAWS1(I)
 9028  CONTINUE
      ENDIF
!
!   CALCULATE D PSIM AND D PSIH
!
      IF ( (ITYPE.EQ.1) .AND. LWATER ) THEN
       DO 9032 I = 1,IRUN
        IF (IWATER(I).EQ.1) THEN
         VDPSIM(I) = - VDZ0(I) * VAZ0(I)
         VDPSIH(I) = VDPSIM(I)
        ENDIF
 9032  CONTINUE
      ENDIF
!
      IF (ITYPE.GE.3) THEN
       DO 9034 I = 1,IRUN
        VDPSIM(I) = VDX(I) * VDZETA(I)
        VDPSIH(I) = VDY(I) * VDZETA(I)
        IF ( LWATER ) THEN
         IF (IWATER(I).EQ.1 ) THEN
          VDPSIM(I) = VDPSIM(I) - VX0(I) * VDZ0(I)
          VDPSIH(I) = VDPSIH(I) - VY0(I) * VDZ0(I)
         ENDIF
        ENDIF
 9034  CONTINUE
      ENDIF
!
!   PREVENT OVERCORRECTION OF PSIM OR PSIH FOR UNSTABLE CASE
!
      IF (ITYPE.GE.4) THEN
       DO 9036 I = 1,IRUN
        VDPSIMC(I) = -0.9 - VDPSIM(I) * VAPSIM(I)
        VDPSIHC(I) = -0.9 *  VPSIH(I) - VDPSIH(I)
        IF ( VDPSIMC(I).GT.0.  ) THEN
         VINT1(I) = 1
        ELSE
         VINT1(I) = 0
        ENDIF
        IF ( VDPSIHC(I).GT.0.  ) THEN
         VINT2(I) = 1
        ELSE
         VINT2(I) = 0
        ENDIF
        VDZETA1(I) = 0.
        IF(VINT1(I).EQ.1) VDZETA1(I) = VDPSIMC(I) / VDXPSIM(I)
        IF((VINT1(I).EQ.1).OR.(VINT2(I).EQ.1)) VTEMPLIN(I) = &
              VDY(I) + VY0(I) * VG(I) * VDXPSIM(I)
        IF (VINT2(I).EQ.1) then
             VDZETA2(I) =  VDPSIHC(I) / VTEMPLIN(I)
        IF ( VDZETA2(I).LT.VDZETA1(I) ) VDZETA1(I) = VDZETA2(I)
        endif
        IF((VINT1(I).EQ.1).OR.(VINT2(I).EQ.1)) THEN
         VDZETA(I) = VDZETA1(I) + VDZETA(I)
         VDPSIM(I) = VDPSIM(I) + VDX(I) * VR1MG0(I) * VDZETA1(I)
         VDPSIH(I) = VDPSIH(I) + VTEMPLIN(I) * VDZETA1(I)
         IF ( IWATER(I).EQ.1 ) &
           VDZ0(I) = VDZ0(I) - VG(I) * VDXPSIM(I) * VDZETA1(I)
        ENDIF
 9036  CONTINUE
      ENDIF
!
end subroutine linadj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: zcsub
! !INTERFACE:
      SUBROUTINE ZCSUB (VUSTAR,VDZSEA,IWATER,LDZSEA,IRUN,VZSEA,CHOOSEZ0)
!**********************************************************************
!  FUNCTION ZSEA
!  PURPOSE
!     COMPUTES Z0 AS A FUNCTION OF USTAR OVER WATER SURFACES
!  USAGE
!     CALLED BY helfsurface
!  DESCRIPTION OF PARAMETERS
!     USTAR    -  INPUTED VALUE OF SURFACE-STRESS VELOCITY
!     DZSEA    -  OUTPUTED VALUE OF DERIVATIVE  D(ZSEA)/D(USTAR)
!     WATER    -  INPUTED BIT VECTOR TO DETERMINE WATER POINTS
!     LDZSEA   -  LOGICAL FLAG TO DETERMINE IF DZSEA SHOULD BE COMPUTED
!     ZSEA     -  OUTPUTED VALUE OF ROUGHNESS LENGTH
!     CHOOSEZ0 -  INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                               1 - Edson Z0 for mom. and heat, high wind limit
!                               2 - L&P Z0, high wind limit
!                               3 - Edson Z0 for mom. only, high wind limit
!  SUBPROGRAMS NEEDED
!     NONE
!  RECORD OF MODIFICATIONS
!   Molod 6/8/2011 - Implement new choozez0 options (expand from 0,1 choice)
!  REMARKS:
!        COMPUTE ROUGHNESS LENGTH FOR OCEAN POINTS
!          BASED ON FUNCTIONS OF LARGE AND POND
!          AND OF KONDO --- DESIGNED FOR K = .4
! *********************************************************************
      implicit none 

! Argument List Delcarations
      integer irun, CHOOSEZ0
      real VZSEA(IRUN),VUSTAR(IRUN),VDZSEA(IRUN)
      integer IWATER(IRUN)
      LOGICAL LDZSEA

! Local Variables
      real USTMX1_OLD,USTMX2_OLD
      real USTMX1_NEW,USTMX2_NEW
      real USTMX1,USTMX2,USTMX3

      PARAMETER ( USTMX1_NEW =   0.80 )
      PARAMETER ( USTMX2_NEW =   0.80 )
      PARAMETER ( USTMX1_OLD =   1.1  )
      PARAMETER ( USTMX2_OLD =   0.381844 )
      PARAMETER ( USTMX3     =   0.0632456)

      real AA(IRUN,5),TEMP(IRUN)
      integer INT2(IRUN),INT3(IRUN),INT4(IRUN)
      integer i,k
      real ustloc(irun)

      real AA1(5),AA2(5),AA3(5),AA4(5)
      real AA2_NEW(5),AA3_NEW(5),AA4_NEW(5)
      real AA2_OLD(5),AA3_OLD(5),AA4_OLD(5)

      DATA AA1/.2030325E-5,0.0,0.0,0.0,0.0/

      DATA AA2_NEW/-1.102451E-08,0.1593E-04,0.1E-03,2.918E-03, &
               0.695649E-04/
      DATA AA3_NEW/-1.102451E-08,0.12E-04,0.1E-03,2.918E-03, &
               1.5649E-04/
      DATA AA4_NEW/0.085E-03,1.5E-03,-0.210E-03,0.215E-02, &
               -0.0/

      DATA AA2_OLD/-0.402451E-08,0.239597E-04,0.117484E-03,0.191918E-03, &
               0.395649E-04/
      DATA AA3_OLD/-0.237910E-04,0.228221E-03,-0.860810E-03,0.176543E-02, &
               0.784260E-04/
      DATA AA4_OLD/-0.343228E-04,0.552305E-03,-0.167541E-02,0.250208E-02, &
               -0.153259E-03/

      if( CHOOSEZ0.eq.0 .OR. CHOOSEZ0.eq.2) then
          USTMX1 = USTMX1_OLD
          USTMX2 = USTMX2_OLD
             AA2 =    AA2_OLD
             AA3 =    AA3_OLD
             AA4 =    AA4_OLD
      else
          USTMX1 = USTMX1_NEW
          USTMX2 = USTMX2_NEW
             AA2 =    AA2_NEW
             AA3 =    AA3_NEW
             AA4 =    AA4_NEW
      endif
!
!**********************************************************************
!*****              LOWER CUTOFF CONDITION FOR USTAR                ***
!**********************************************************************
!
      DO 9000 I = 1,IRUN
       IF(VUSTAR(I) .LT. 1.e-6)THEN
        INT3(I) = 1
       ELSE
        INT3(I) = 0
       ENDIF
 9000 CONTINUE
      DO 9002 I = 1,IRUN
       IF(INT3(I).EQ.1) VUSTAR(I) = 1.e-6
 9002 CONTINUE

      ustloc = vustar
!
!***********************************
!*****  LOAD THE ARRAY A(I,K)  *****
!***********************************
!
      DO 9004 I = 1,IRUN
       IF( (ustloc(I) .GT. USTMX1) .AND. (IWATER(I).EQ.1) ) THEN
        if( CHOOSEZ0.gt.0 ) ustloc(i) = ustmx1
        INT4(I) = 1
       ELSE
        INT4(I) = 0
       ENDIF
 9004 CONTINUE
      DO 9006 I = 1,IRUN
       IF(ustloc(I) .GT. USTMX2) THEN
        INT3(I) = 1
       ELSE
        INT3(I) = 0
       ENDIF
 9006 CONTINUE
      DO 9008 I = 1,IRUN
       IF(ustloc(I) .GE. USTMX3) THEN
        INT2(I) = 1
       ELSE
        INT2(I) = 0
       ENDIF
 9008 CONTINUE
!
      DO 100 K=1,5
       DO 9010 I = 1,IRUN
        AA(I,K) = AA1(K)
        IF( INT2(I).EQ.1 )  AA(I,K) = AA2(K)
        IF( INT3(I).EQ.1 )  AA(I,K) = AA3(K)
        IF( INT4(I).EQ.1 )  AA(I,K) = AA4(K)
 9010  CONTINUE
  100 CONTINUE
!
!********************************************************
!*****  EVALUATE THE ENHANCED POLYNOMIAL FOR ZSEA  *****
!********************************************************
!
      DO 9012 I = 1,IRUN
       VDZSEA(I)  =  ( AA(I,4) + AA(I,5) * ustloc(I) ) * ustloc(I)
       VZSEA(I)  =  AA(I,2) + ( AA(I,3) + VDZSEA(I) ) * ustloc(I)
       TEMP(I) = AA(I,1) / ustloc(I)
       VZSEA(I)  =  VZSEA(I) + TEMP(I)
 9012 CONTINUE
!
!**********************************************************************
!*****        EVALUATE THE DERIVATIVE DZSEA IF LDZSEA IS TRUE       ***
!**********************************************************************
!
      IF( LDZSEA ) THEN
       DO 9014 I = 1,IRUN
        VDZSEA(I)  =  3. * VDZSEA(I) -(AA(I,4)*ustloc(I) - AA(I,3))
        VDZSEA(I)  =  VDZSEA(I) * ustloc(I) - TEMP(I)
 9014  CONTINUE
      ENDIF
!
end subroutine zcsub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module clsmf25_surfacelayer
