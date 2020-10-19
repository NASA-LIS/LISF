!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
MODULE MODI_CLOSE_FILE
CONTAINS
      SUBROUTINE MODI_CLOSE_FILE_SUB(HPROGRAM,KUNIT)
!     #######################################################
!
!!****  *CLOSE_FILE* - generic routine to close a file
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#if defined(SFX_ASC) || defined(SFX_ARO) || defined(SFX_MNH) || defined(SFX_NC)
USE MODI_CLOSE_FILE_ASC
#endif
#ifdef SFX_FA
USE MODI_CLOSE_FILE_FA
#endif
#ifdef SFX_OL
USE MODI_CLOSE_FILE_OL
#endif
#ifdef SFX_LFI
USE MODI_CLOSE_FILE_LFI
#endif
!
#ifdef SFX_MNH
USE MODI_CLOSE_FILE_MNH
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(IN)  :: KUNIT    ! logical unit of file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CLOSE_FILE',0,ZHOOK_HANDLE)
IF (HPROGRAM=='MESONH') THEN
#ifdef SFX_MNH
  CALL CLOSE_FILE_MNH(HPROGRAM,KUNIT)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
#ifdef SFX_OL
  CALL CLOSE_FILE_OL(HPROGRAM,KUNIT)
#endif
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef SFX_FA
  CALL CLOSE_FILE_FA(HPROGRAM,KUNIT)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef SFX_LFI
  CALL CLOSE_FILE_LFI(HPROGRAM,KUNIT)
#endif
ELSE 
#if defined(SFX_ASC) || defined(SFX_ARO) || defined(SFX_MNH) || defined(SFX_NC)
  CALL CLOSE_FILE_ASC(HPROGRAM,KUNIT)
#endif
END IF
IF (LHOOK) CALL DR_HOOK('CLOSE_FILE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MODI_CLOSE_FILE_SUB
END MODULE MODI_CLOSE_FILE
