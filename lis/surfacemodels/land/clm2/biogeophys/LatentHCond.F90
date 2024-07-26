!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

subroutine LatentHCond (raw,   rbw,   rdw,   rpp,   wtaq,  &
                        wtlq,  wtgq,  wtaq0, wtlq0, wtgq0, &
                        wtalq, wtgaq, wtglq, clm    )

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Provides dimensional and non-dimensional latent heat 
! conductances for canopy and soil flux calculations.  Latent fluxes 
! differs from the sensible heat flux due to stomatal resistance.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: LatentHCond.F90,v 1.6 2004/11/24 22:56:30 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

  real(r8), intent(in) :: raw    ! aerodynamical resistance [s/m]
  real(r8), intent(in) :: rbw    ! leaf boundary layer resistance [s/m]
  real(r8), intent(in) :: rdw    ! latent heat resistance between ground and bottom 
                                 ! of canopy
  real(r8), intent(in) :: rpp    ! fraction of potential evaporation from leaf [-]

  real(r8), intent(out) :: wtaq  ! latent heat conduactance for air [m s-1]
  real(r8), intent(out) :: wtlq  ! latent heat conduactance for leaf [m s-1]
  real(r8), intent(out) :: wtgq  ! latent heat conduactance for ground [m s-1]
  real(r8), intent(out) :: wtaq0 ! normalized latent heat conduactance for air [-]
  real(r8), intent(out) :: wtlq0 ! normalized latent heat conduactance for leaf [-]
  real(r8), intent(out) :: wtgq0 ! normalized heat conduactance for ground [-]
  real(r8), intent(out) :: wtalq ! normalized latent heat cond. for air and leaf [-]
  real(r8), intent(out) :: wtglq ! normalized latent heat cond. for leaf and ground [-]
  real(r8), intent(out) :: wtgaq ! normalized latent heat cond. for air and ground [-]

!----Local Variables----------------------------------------------------

  real(r8) wtsqi                 ! latent heat resistance for air, grd and leaf [-]

!----End Variable List--------------------------------------------------

  wtaq  = clm%frac_veg_nosno/raw                                ! air
  wtlq  = clm%frac_veg_nosno*(clm%elai+clm%esai)/rbw * rpp      ! leaf
  wtgq  = clm%frac_veg_nosno/rdw                                ! ground
  wtsqi = 1./(wtaq+wtlq+wtgq)

  wtgq0 = wtgq*wtsqi                    ! ground
  wtlq0 = wtlq*wtsqi                    ! leaf
  wtaq0 = wtaq*wtsqi                    ! air

  wtglq = wtgq0+wtlq0                   ! ground + leaf
  wtgaq = wtaq0+wtgq0                   ! air + ground
  wtalq = wtaq0+wtlq0                   ! air + leaf

end subroutine LatentHCond
