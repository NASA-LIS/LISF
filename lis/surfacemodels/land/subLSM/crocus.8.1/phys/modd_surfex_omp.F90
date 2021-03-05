!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################
      MODULE MODD_SURFEX_OMP
!     ######################
!
!!****  *MODD_SURFEX_OMP
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      S. Faroux   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       26/06/12
!!      Modified    11/2013 by J.Escobar :add !$ to inhibit completly omp
!!                                 dependency
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
! MN 
!#ifdef AIX64 
! USE OMP_LIB
!#endif
!
IMPLICIT NONE
!
!#ifndef AIX64
!  INCLUDE 'omp_lib.h'
!#endif
 ! MN
!
!RJ: this broke non openmp version before
!RJ: OMP_GET_THREAD_NUM() returns 0 for first omp thread
!RJ: OMP_GET_NUM_THREADS() returns 1 for omp thread count
#ifdef RJ_OFIX
INTEGER :: NBLOCKTOT = 1
INTEGER :: NBLOCK = 0
#else
INTEGER :: NBLOCKTOT = 1
INTEGER :: NBLOCK = 1
#endif
!$OMP THREADPRIVATE(NBLOCK)
INTEGER :: IDC = 0
!
END MODULE MODD_SURFEX_OMP

