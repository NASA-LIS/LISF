!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #############################################################
MODULE MODI_ABOR1_SFX
CONTAINS
      SUBROUTINE ABOR1_SFX_SUB(YTEXT)
!     #############################################################
!
!!****  *ABOR1_SFX* - abor1 subroutine
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2008 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC
USE MODD_SURFEX_OMP, ONLY : NBLOCK, NBLOCKTOT
USE MODD_SURF_CONF,  ONLY : CPROGNAME, CSOFTWARE
!
USE MODI_GET_LUOUT
USE MODI_CLOSE_FILE
!      
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=*),  INTENT(IN)  :: YTEXT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=6)  :: YPROGRAM   
 CHARACTER(LEN=20) :: YSTRING
INTEGER           :: ILUOUT         ! logical unit of output file      
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
#ifdef SFX_ARO
#include "abor1.intfb.h"
#endif
!-------------------------------------------------------------------------------
!
!* get output listing file logical unit
!
IF (LHOOK) CALL DR_HOOK('ABOR1_SFX',0,ZHOOK_HANDLE)
YPROGRAM = CPROGNAME
!      
 CALL MODI_GET_LUOUT_SUB(YPROGRAM,ILUOUT)
!
IF (YPROGRAM=='ASCII ' .OR. YPROGRAM=='TEXTE ' .OR. YPROGRAM=='BINARY' .OR. YPROGRAM=='NC    ') THEN
   IF ( NPROC>1 .OR. NBLOCKTOT>1 ) &
     WRITE(*,*)"MPI TASK NUMBER = ",NRANK,", OMP THREAD NUMBER = ",NBLOCK
   WRITE(*,*)YTEXT
   YSTRING='LISTING_'//TRIM(CSOFTWARE)//'.txt'
   WRITE(*,*)'-------------------------------------------------------------------------------'
   WRITE(*,*) 'MORE DETAILS ABOUT THE CRASH IN THE OUTPUT LISTING FILE: ', TRIM(YSTRING)
   WRITE(*,*)'-------------------------------------------------------------------------------'   
ENDIF
!
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '--------------------   FATAL ERROR in SURFEX  -----------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '-                                                                         -'
WRITE(ILUOUT,*)YTEXT
WRITE(ILUOUT,*) '-                                                                         -'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
WRITE(ILUOUT,*) '---------------------------------------------------------------------------'
 CALL MODI_CLOSE_FILE_SUB(YPROGRAM,ILUOUT)
!
#ifdef SFX_ARO
call abor1('abort by abor1_sfx')
#else
 write(0,*) "aborted with text:",trim(ytext),"|"
 CALL ABORT
STOP
#endif
IF (LHOOK) CALL DR_HOOK('ABOR1_SFX',1,ZHOOK_HANDLE)
!
END SUBROUTINE ABOR1_SFX_SUB
END MODULE MODI_ABOR1_SFX
