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

subroutine Kill (bm_inc, litter_ag, litter_bg, &
                 lm_ind, sm_ind   , hm_ind   , &
                 rm_ind, nind     , present  , tree)

!----------------------------------------------------------------------- 
! 
! Purpose: Removal of PFTs with negative annual C increment
!          NB: PFTs newly beyond their bioclimatic limits are removed in
!          subroutine establishment
! 
! Method: Called once per year
! 
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. kill)
! 
!-----------------------------------------------------------------------
! $Id: Kill.F90,v 1.6 2004/11/24 22:56:56 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  implicit none

! ----------------------------- arguments ------------------------------
  logical , intent(inout) :: present
  real(r8), intent(inout) :: litter_ag
  real(r8), intent(inout) :: litter_bg
  real(r8), intent(in)    :: bm_inc
  real(r8), intent(in)    :: lm_ind
  real(r8), intent(in)    :: sm_ind
  real(r8), intent(in)    :: hm_ind
  real(r8), intent(in)    :: rm_ind
  real(r8), intent(in)    :: nind
  logical , intent(in)    :: tree
! ----------------------------------------------------------------------

  if (present) then

     if (bm_inc < 0.0) then !negative C increment this year

        present = .false.   !remove PFT

! Transfer killed biomass to litter

        if (tree) then !redundant if block? (slevis)
           litter_ag = litter_ag + (lm_ind + sm_ind + hm_ind) * nind
        else !if grass
           litter_ag = litter_ag + lm_ind * nind
        endif

        litter_bg = litter_bg + rm_ind * nind

     endif

  endif

  return
end subroutine Kill
