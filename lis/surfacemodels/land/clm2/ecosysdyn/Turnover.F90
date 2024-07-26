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

subroutine Turnover (pftpar, litter_ag, litter_bg, &
                     lm_ind, sm_ind   , hm_ind   , &
                     rm_ind, nind     , present  , turnover_ind)

!----------------------------------------------------------------------- 
! 
! Purpose: Turnover of PFT-specific fraction from each living C pool
!          Leaf and root C transferred to litter, sapwood C to heartwood
! 
! Method: Called once per year
! 
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. turnover)
! 
!-----------------------------------------------------------------------
! $Id: Turnover.F90,v 1.6 2004/11/24 22:57:02 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_varpar, ONLY: npftpar
  implicit none

! ----------------------------- arguments ------------------------------
  real(r8), intent(out)   :: turnover_ind
  real(r8), intent(inout) :: litter_ag
  real(r8), intent(inout) :: litter_bg
  real(r8), intent(inout) :: lm_ind
  real(r8), intent(inout) :: sm_ind
  real(r8), intent(inout) :: hm_ind
  real(r8), intent(inout) :: rm_ind
  real(r8), intent(in)    :: nind
  real(r8), intent(in)    :: pftpar(npftpar)
  logical , intent(in)    :: present
! ----------------------------------------------------------------------

! --------------------------- local variables --------------------------
  real(r8) :: l_torate
  real(r8) :: s_torate
  real(r8) :: r_torate
  real(r8) :: lm_turn
  real(r8) :: sm_turn
  real(r8) :: rm_turn
! ----------------------------------------------------------------------

! Turnover rates are reciprocals of tissue longevity

  l_torate = 1.0 / pftpar(9)
  s_torate = 1.0 / pftpar(11)
  r_torate = 1.0 / pftpar(12)

  if (present) then

! Calculate the biomass turnover in this year

     lm_turn = lm_ind * l_torate
     sm_turn = sm_ind * s_torate
     rm_turn = rm_ind * r_torate

! Update the pools

     lm_ind = lm_ind - lm_turn
     sm_ind = sm_ind - sm_turn
     rm_ind = rm_ind - rm_turn

! Convert sapwood to heartwood

     hm_ind = hm_ind + sm_turn

! Transfer to litter pools

     litter_ag = litter_ag + lm_turn * nind
     litter_bg = litter_bg + rm_turn * nind

! Record total turnover

     turnover_ind = lm_turn + sm_turn + rm_turn

  endif

  return
end subroutine Turnover
