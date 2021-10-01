!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
MODULE MODI_WIND_THRESHOLD
CONTAINS
    FUNCTION WIND_THRESHOLD(PWIND,PUREF) RESULT(PWIND_NEW)
!   ############################################################################
!
!!****  *WIND_THRESHOLD*  
!!
!!    PURPOSE
!!    -------
!
!     Set a minimum value to the wind for exchange coefficient computations.
!     This minimum value depends on the forcing height
!     
!!    AUTHOR
!!    ------
!!
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2007 
!-------------------------------------------------------------------------------
!
USE MODD_SURF_ATM, ONLY: XCISMIN, XVMODMIN, LALDTHRES
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN)   :: PWIND      ! wind
REAL, DIMENSION(:), INTENT(IN)   :: PUREF      ! forcing level
!
REAL, DIMENSION(SIZE(PWIND))     :: PWIND_NEW  ! modified wind
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------
!
!  wind gradient
!
IF (LHOOK) CALL DR_HOOK('WIND_THRESHOLD',0,ZHOOK_HANDLE)
IF (.NOT.LALDTHRES) THEN
!        
!  minimum value for exchange coefficients computations : 1m/s / 10m
   PWIND_NEW = MAX(PWIND , 0.1 * MIN(10.,PUREF) )
ELSE
!  minimum value for exchange coefficients computations : 1m/s / 10m
   PWIND_NEW = MAX( XVMODMIN, SQRT( PWIND**2 + (XCISMIN*PUREF)**2 ) )
ENDIF
IF (LHOOK) CALL DR_HOOK('WIND_THRESHOLD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END FUNCTION WIND_THRESHOLD
END MODULE MODI_WIND_THRESHOLD
