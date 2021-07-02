!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################
      MODULE MODD_REPROD_OPER
!     ######################
!
!!****  *MODD_REPROD_OPER* - declaration of ISBA parameters
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify temporary
!       old parameters related to the surface parameterization ISBA 
!       to ensure reproductibility with previous oper cycle
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
!!      Original       04/2013                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!------------------------------------------------------------------------------------
! Old global ISBA param temporary activated in NAM_SURF_REPROD_OPER (for reproductibility)
!------------------------------------------------------------------------------------
!
! * Tropical evergreen forest parameter
!
!XEVERG_RSMIN : old = 250. (Manzi 1993) but observations range 
!               from 140 to 180. According to Delire et al. (1997) and 
!               new tests over 6 local sites, 175. is recommended
!               Should be the default after check with AROME/ALADIN
!
REAL             :: XEVERG_RSMIN
!
!XEVERG_VEG : old = 0.99 (Manzi 1993) but according to Delire et al. (1997) and 
!             new tests over 6 local sites, 1.0 is recommended because 0.99
!             induces unrealistic bare soil evaporation for Tropical forest
!             Should be the default after check with AROME/ALADIN
!
REAL             :: XEVERG_VEG
!
! * Soil depth average
!
!CDGAVG : old         = 'ARI' Arithmetic average for all depths 
!         recommended = 'INV' Harmonic average for all depths
!
CHARACTER(LEN=3) :: CDGAVG
!
! * Soil depth with ISBA-DF
!
!CDGDIF : old         = 'SOIL' Total soil depth (d3) in Ecoclimap
!         recommended = 'ROOT' Root depth (d2) in Ecoclimap
!
CHARACTER(LEN=4) :: CDGDIF
!
! * wind implicitation
!
CHARACTER(LEN=3) :: CIMPLICIT_WIND ! wind implicitation option
!                                  ! 'OLD' = direct
!                                  ! 'NEW' = Taylor serie, order 1 (recommended)
!
! * qsat computation
!
CHARACTER(LEN=3) :: CQSAT ! qsat computation option
!                         ! 'OLD' = do not depend on temperature
!                         ! 'NEW' = qsat and qsati merged (recommended)
!
! * Charnock parameter
!
CHARACTER(LEN=3) :: CCHARNOCK ! Charnock parameter option
!                             ! 'OLD' = constant equal to XVCHRNK
!                             ! 'NEW' = vary between 0.011 et 0.018 according
!                             !         to Chris Fairall's data as in coare3.0
!                             !         (recommended)
!
!--------------------------------------------------------------------------------
!
END MODULE MODD_REPROD_OPER












