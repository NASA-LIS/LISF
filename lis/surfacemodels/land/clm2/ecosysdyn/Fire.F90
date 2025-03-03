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

subroutine Fire(fire_length, litter_ag , present , &
                tree       , resist    , nind    , &
                lm_ind     , sm_ind    , hm_ind  , &
                rm_ind     , afire_frac, fpc_grid, &
                acflux_fire) 

!----------------------------------------------------------------------- 
! 
! Purpose:
! 
! Method: Called once per year
! 
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine fire)
! 
!-----------------------------------------------------------------------
! $Id: Fire.F90,v 1.6 2004/11/24 22:56:54 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  implicit none

! ----------------------------- arguments ------------------------------
  real(r8), intent(out)   :: afire_frac
  real(r8), intent(out)   :: acflux_fire
  real(r8), intent(inout) :: litter_ag
  real(r8), intent(inout) :: nind
  real(r8), intent(in) :: fire_length
  real(r8), intent(in) :: resist
  real(r8), intent(in) :: lm_ind
  real(r8), intent(in) :: sm_ind
  real(r8), intent(in) :: hm_ind
  real(r8), intent(in) :: rm_ind
  real(r8), intent(in) :: fpc_grid
  logical , intent(in) :: present
  logical , intent(in) :: tree
! ----------------------------------------------------------------------

! --------------------------- local variables --------------------------
  real(r8), parameter :: minfuel = 200.0 !fuel threshold to carry a fire (gC/m2)
  real(r8) :: fire_index
  real(r8) :: fire_term
  real(r8) :: disturb
  real(r8) :: fuel                       !=litter_ag_total per grid cell
! ----------------------------------------------------------------------

! slevis: Orig. had a daily loop to calculate fire_length
!         Now fire_length comes from subroutine FireSeason

! Calculate annual fire index

  fire_index = fire_length / 365.0

! Calculate the fraction of the grid cell affected by fire

  fire_term = fire_index - 1.0
  afire_frac = max(fire_index *                                    &
                   exp(fire_term /                                 &
                       (-0.13*fire_term**3 + 0.6*fire_term**2 +    &
                          0.8*fire_term    + 0.45              )), &
                   0.001)

! Reduce fraction of patch affected by fire when fuel
! becomes limiting (reduced carrying capacity)

  if (litter_ag < minfuel * fpc_grid) afire_frac = 0.001

! Implement the effect of the fire on vegetation structure and litter
! in the disturbed fraction.

! Each PFT is assigned a resistance to fire, representing the fraction of
! the PFT which survives a fire. Grasses assumed already to have completed
! their life cycle and thus are not affected by fire, giving them
! a competitive advantage against woody PFTs.

  if (present .and. tree) then

! Calculate the fraction of individuals in grid cell which die
! (slevis: 'in grid cell' because nind is grid average)

     disturb = (1.0-resist) * afire_frac

! Calculate carbon flux to atmosphere (gC/m2) due to burnt biomass

     acflux_fire = disturb * nind * (lm_ind+sm_ind+hm_ind+rm_ind)

! Update the individual density

     nind = nind * (1.0-disturb)

! Add combusted litter to carbon flux to atmosphere term

     acflux_fire = acflux_fire + afire_frac * litter_ag

  else
     acflux_fire =               afire_frac * litter_ag
  endif

! Update the above ground litter term

  litter_ag = (1.0 - afire_frac) * litter_ag

  return
end subroutine Fire
