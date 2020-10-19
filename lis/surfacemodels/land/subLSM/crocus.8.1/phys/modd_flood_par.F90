!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################
      MODULE MODD_FLOOD_PAR
!     ######################
!
!!****  *MODD_FLOOD_PAR* - declaration of parameters related
!!                          to the flood parameterization
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     parameters related to the surface parameterization of flood.
!
!!
!!      
!!
!!    AUTHOR
!!    ------
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       02/2010                     
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL, SAVE       :: XCFFV
!                   Coefficient for calculation of floodplain fraction over vegetation
!
REAL, SAVE       :: XZ0FLOOD
!                   Roughness length for flood (m)
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_FLOOD_PAR












