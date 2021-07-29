!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ######################################################################

MODULE MODI_SURFACE_AERO_COND

CONTAINS

    SUBROUTINE MODI_SURFACE_AERO_COND_SUB(PRI, PZREF, PUREF, PVMOD, PZ0,&
                                     PZ0H, PAC, PRA, PCH    ,HSNOWRES , &
                                     KSIZE1) 
!   ######################################################################
!
!!****  *SURFACE_AERO_COND*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the drag coefficients for heat and momentum near the ground
!         
!     
!!**  METHOD
!!    ------
!
!
!
!    1 and 2 : computation of relative humidity near the ground
!
!    3 : richardson number
!
!    4 : the aerodynamical resistance for heat transfers is deduced
!
!    5 : the drag coefficient for momentum ZCD is computed
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/01/98 
!!                  02/04/01 (P Jabouille) limitation of Z0 with 0.5 PUREF
!!                  05/2016 M. Lafaysse - B. Cluzet : implement Martin and Lejeune 1998 formulation for multiphysics
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XKARMAN
USE MODI_WIND_THRESHOLD
!
USE MODE_THERMOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(1:KSIZE1), INTENT(IN)    :: PRI      ! Richardson number
REAL, DIMENSION(1:KSIZE1), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(1:KSIZE1), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level
REAL, DIMENSION(1:KSIZE1), INTENT(IN)    :: PUREF    ! reference height of the wind
                                              ! NOTE this is different from ZZREF
                                              ! ONLY in stand-alone/forced mode,
                                              ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(1:KSIZE1), INTENT(IN)    :: PZ0      ! roughness length for momentum
REAL, DIMENSION(1:KSIZE1), INTENT(IN)    :: PZ0H     ! roughness length for heat
!
REAL, DIMENSION(1:KSIZE1), INTENT(OUT)   :: PAC      ! aerodynamical conductance
REAL, DIMENSION(1:KSIZE1), INTENT(OUT)   :: PRA      ! aerodynamical resistance
REAL, DIMENSION(1:KSIZE1), INTENT(OUT)   :: PCH      ! drag coefficient for heat
CHARACTER(LEN=3), INTENT(IN)  ::HSNOWRES !surface exchange coefficient option 
INTEGER , INTENT(IN)   :: KSIZE1
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PRI)) :: ZZ0, ZZ0H, ZMU,          &
                               ZFH, ZCHSTAR, ZPH, ZCDN, &
                               ZSTA, ZDI, ZWORK1, ZWORK2, ZWORK3 
REAL, DIMENSION(SIZE(PRI)) :: ZVMOD, ZMARTIN,ZCDN_M98
!
INTEGER                    :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Functions:
REAL :: X, CHSTAR, PH
CHSTAR(X) = 3.2165 + 4.3431*X + 0.5360*X*X - 0.0781*X*X*X
PH    (X) = 0.5802 - 0.1571*X + 0.0327*X*X - 0.0026*X*X*X
!
!-------------------------------------------------------------------------------
!
!*       4.     Surface aerodynamic resistance for heat transfers
!               -------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SURFACE_AERO_COND',0,ZHOOK_HANDLE)
ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
!
DO JJ=1,1 ! SIZE(PRI) we need the first memeber of the array the rest is dummy 
  ZZ0(JJ)  = MIN(PZ0(JJ),PUREF(JJ)*0.5)
  ZZ0H(JJ) = MIN(ZZ0(JJ),PZ0H(JJ))
  ZZ0H(JJ) = MIN(ZZ0H(JJ),PZREF(JJ)*0.5)
!
  ZWORK1(JJ)=LOG( PUREF(JJ)/ZZ0(JJ) )
  ZWORK2(JJ)=PZREF(JJ)/ZZ0H(JJ)
  ZWORK3(JJ)=ZVMOD(JJ)*ZVMOD(JJ)
  ZMU(JJ) = MAX( LOG( ZZ0(JJ)/ZZ0H(JJ) ), 0.0 )
  ZFH(JJ) = ZWORK1(JJ) / LOG(ZWORK2(JJ))
!
  ZCHSTAR(JJ) = CHSTAR(ZMU(JJ))
  ZPH(JJ)     = PH(ZMU(JJ))
!
! 
  ZCDN(JJ) = (XKARMAN/ZWORK1(JJ))**2.
!
!
  ZSTA(JJ) = PRI(JJ)*ZWORK3(JJ)
  ZCDN_M98(JJ)= XKARMAN*XKARMAN/(LOG(PUREF(JJ)/ZZ0(JJ))*LOG(PZREF(JJ)/ZZ0(JJ)))
   
   
  IF(HSNOWRES=='RIL' .OR. HSNOWRES=='DEF') THEN
!
      IF ( PRI(JJ) < 0.0 ) THEN
        ZDI(JJ) = 1. / ( ZVMOD(JJ)                                  &
                       +ZCHSTAR(JJ)*ZCDN(JJ)*15.                         &
                                    *ZWORK2(JJ)**ZPH(JJ)  &
                                    *ZFH(JJ) * SQRT(-ZSTA(JJ))           &
                      ) 
        PAC(JJ) = ZCDN(JJ)*  (  ZVMOD(JJ)-15.* ZSTA(JJ)*ZDI(JJ)  )  *  ZFH(JJ)

      ELSE
        ZDI(JJ) = SQRT(ZWORK3(JJ) + 5. * ZSTA(JJ) )
        PAC(JJ) = ZCDN(JJ)*ZVMOD(JJ)/(1.+15.*ZSTA(JJ)*ZDI(JJ)  &
                 / ZWORK3(JJ) /ZVMOD(JJ) )*ZFH(JJ)    
      ENDIF
    !
      PRA(JJ) = 1. / PAC(JJ)
    !
      PCH(JJ) = 1. / (PRA(JJ) * ZVMOD(JJ))
  ELSE IF (HSNOWRES=='M98')THEN
  ! Martin and Lejeune 1998 ; Cluzet et al 2016
    IF (PRI(JJ)<0.0) THEN
        IF (ZCDN_M98(JJ)==0.) THEN
            ZMARTIN(JJ)=1.
        ELSE
            ZMARTIN(JJ)= 1.  + (7./  ( 0.83*(ZCDN_M98(JJ))**(-0.62) )  ) &
            *LOG(1.-0.83*(ZCDN_M98(JJ))**(-0.62)   *  PRI(JJ))
        ENDIF
           
    ELSE
        IF (PRI(JJ)<0.2) THEN
            ZMARTIN(JJ) = MAX(0.75,( 1. - 5.*PRI(JJ))*(1. - 5. *PRI(JJ)))!
            !Nota B. Cluzet : le min servait à empêcher le CH de remonter pour 
            !des Ri >0.4 c'est une erreur car cela seuille le Ch dès Ri =0 au lieu de le faire dès Ri=0.026
        ELSE
            ZMARTIN(JJ)=MIN(MAX(0.75,( 1. - 5.*PRI(JJ))*(1. - 5. *PRI(JJ))),0.75)
        ENDIF
    ENDIF
    
    PCH(JJ) = ZMARTIN(JJ) * ZCDN_M98(JJ)
    PRA(JJ)=1./(PCH(JJ)*ZVMOD(JJ))! Nota B. Cluzet : checked in noilhan and mahfouf seems ok

    PAC(JJ)=1./(PRA(JJ))
  ENDIF

ENDDO
IF (LHOOK) CALL DR_HOOK('SURFACE_AERO_COND',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MODI_SURFACE_AERO_COND_SUB
END MODULE MODI_SURFACE_AERO_COND
