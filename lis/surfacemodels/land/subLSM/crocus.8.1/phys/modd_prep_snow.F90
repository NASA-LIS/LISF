!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################
      MODULE MODD_PREP_SNOW
!     ################
!
!!****  *MODD_PREP - declaration for field interpolations
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
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!      Modif M Lafaysse 04/2014 : LSNOW_PREP_PERM
!!       Modified by F. Tuzet (06/2016): Add of a new dimension for impurity: The type of impurity
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_SNOW_PAR
!
IMPLICIT NONE
!
!--------------------------------------------------------------------------
!
!* climatological gradient for vertical extrapolations of snow content
!  a rate of 8cm of snow per degree below 0 C is chosen for these mountain tops
! (climatology from Etchevers 2000 in the Alps and the Jura mountains).
!
REAL, PARAMETER  :: XWSNOW_CLIM_GRAD = - 0.08 * 300.     * (-0.0065)
!
!--------------------------------------------------------------------------
! Parameters for snow field uniforn initialization
!
LOGICAL :: LSNOW_FRAC_TOT
INTEGER, PARAMETER :: NSNOW_LAYER_MAX = 50
LOGICAL :: LSNOW_PREP_PERM ! activate or disactivate initialization over permanent ice areas
INTEGER, PARAMETER               :: NIMPUR_MAX = 5                     !Maximum number of impurity types
 CHARACTER (len=4), DIMENSION(NIMPUR_MAX),PARAMETER ::IMPTYP=(/'Soot','Dust','OrgM','Othr','....'/)
INTEGER                          :: NIMPUR                             ! number of impurity types    
REAL, DIMENSION (NIMPUR_MAX) , PARAMETER            ::SCAVEN_COEF= (/0.0,0.0,0.,0.,0./)           !!Scavenging efficiency of the differrent impurities                      
!
!--------------------------------------------------------------------------
!
!* normalized dimensions for interpolation grids for soil
INTEGER, PARAMETER           :: NGRID_LEVEL = 40
REAL, DIMENSION(NGRID_LEVEL) :: XGRID_SNOW = &
(/0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.10,0.12,0.14,&
  0.16,0.18,0.21,0.25,0.30,0.35,0.40,0.45,0.50,0.55,&
  0.60,0.65,0.70,0.75,0.80,0.85,0.87,0.88,0.89,0.90,&
  0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00/)

!--------------------------------------------------------------------------
!
END MODULE MODD_PREP_SNOW
