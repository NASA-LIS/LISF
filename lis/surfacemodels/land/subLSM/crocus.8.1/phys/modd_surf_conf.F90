!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #####################
      MODULE MODD_SURF_CONF
!     #####################
!
!!****  *MODD_SURF_CONF - surfex configuration
!!
!!    PURPOSE
!!    -------
!     Declaration of program name
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
!!      P. Le Moigne *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/2008
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
!-----------------------------------------------------------------------------------------------------
 CHARACTER(LEN=6) :: CPROGNAME
 CHARACTER(LEN=7) :: CSOFTWARE="       " ! software used: 'PGD    ','PREP  ','OFFLINE','      '
!
END MODULE MODD_SURF_CONF
