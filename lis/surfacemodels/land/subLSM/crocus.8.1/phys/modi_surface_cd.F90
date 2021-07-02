!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   #################################################################
MODULE MODI_SURFACE_CD
CONTAINS
    SUBROUTINE MODI_SURFACE_CD_SUB(PRI, PZREF, PUREF, PZ0EFF, PZ0H,   &
                              PCD, PCDN) 
!   #################################################################
!
!!****  *SURFACE_CD*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the drag coefficients for momentum near the ground
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
!!    MODD_GROUND_PAR
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
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XKARMAN
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
REAL, DIMENSION(:), INTENT(IN)    :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
!                                             ! NOTE this is different from ZZREF
!                                             ! ONLY in stand-alone/forced mode,
!                                             ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:), INTENT(IN)    :: PZ0EFF   ! roughness length for momentum
                                              ! with subgrid-scale orography
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H     ! roughness length for heat
!
REAL, DIMENSION(:), INTENT(OUT)   :: PCD      ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN     ! neutral drag coefficient for momentum
!
!*      0.2    declarations of local variables
!
!
REAL                       :: ZZ0EFF, ZZ0H, ZMU,     &
                               ZCMSTAR, ZPM, ZCM, ZFM 
INTEGER                    :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! Functions :
REAL :: X, CMSTAR, PM
CMSTAR(X) = 6.8741 + 2.6933*X - 0.3601*X*X + 0.0154*X*X*X
PM    (X) = 0.5233 - 0.0815*X + 0.0135*X*X - 0.0010*X*X*X

!-------------------------------------------------------------------------------
!
!*       1.     Drag coefficient for momentum transfers
!               ---------------------------------------
!

!
IF (LHOOK) CALL DR_HOOK('SURFACE_CD',0,ZHOOK_HANDLE)
DO JJ=1,1 ! SIZE(PRI) we need the first memeber of the array the rest is dummy 
  ZZ0EFF = MIN(PZ0EFF(JJ),PUREF(JJ)*0.5)
  ZZ0H   = MIN(ZZ0EFF,PZ0H(JJ))
!
  ZMU = LOG( MIN(ZZ0EFF/ZZ0H,200.) )
!
  PCDN(JJ) = (XKARMAN/LOG(PUREF(JJ)/ZZ0EFF))**2

  ZCMSTAR = CMSTAR(ZMU)
  ZPM     = PM(ZMU)
!
  ZCM = 10.*ZCMSTAR*PCDN(JJ)*( PUREF(JJ)/ZZ0EFF )**ZPM
!
  IF ( PRI(JJ) > 0.0 ) THEN
    ZFM = 1. + 10.*PRI(JJ) / SQRT( 1.+5.*PRI(JJ) )
    ZFM = 1. / ZFM
  ELSE
    ZFM = 1. - 10.*PRI(JJ) / ( 1.+ZCM*SQRT(-PRI(JJ)) )
  ENDIF
!
  PCD(JJ) = PCDN(JJ)*ZFM
!
ENDDO
IF (LHOOK) CALL DR_HOOK('SURFACE_CD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MODI_SURFACE_CD_SUB
END MODULE MODI_SURFACE_CD
