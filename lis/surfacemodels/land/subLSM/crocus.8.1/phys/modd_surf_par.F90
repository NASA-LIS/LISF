!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!####################
MODULE MODD_SURF_PAR
!####################
!
!!****  *MODD_SURF_PAR - declaration of surface parameters
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
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
!!      V. Masson *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       02/2004
!!      J.Escobar     06/2013  for REAL4/8 add EPSILON management
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
!-----------------------------------------------------------------------------------------------------
INTEGER :: NVERSION  ! surface version
INTEGER :: NBUGFIX   ! bugfix number of this version
!
#ifndef SFX_MNH
REAL,    PARAMETER :: XUNDEF = 1.E+20
#else 
#ifdef MNH_MPI_DOUBLE_PRECISION
REAL,    PARAMETER :: XUNDEF = 1.E+20! HUGE(XUNDEF) ! Z'7FFFFFFFFFFFFFFF' !  undefined value
#else
REAL,    PARAMETER :: XUNDEF = 1.E+6 ! HUGE(XUNDEF) ! Z'7FBFFFFF' ! undefined value
#endif
#endif
INTEGER, PARAMETER :: NUNDEF = 1E+9   !  HUGE(NUNDEF) !  undefined value
REAL,    PARAMETER :: XSURF_EPSILON = EPSILON(XSURF_EPSILON)  ! minimum 
REAL,    PARAMETER :: XSURF_HUGE    = HUGE(XSURF_HUGE) 
REAL,    PARAMETER :: XSURF_TINY    = TINY(XSURF_TINY) 
!-----------------------------------------------------------------------------------------------------
!
END MODULE MODD_SURF_PAR
