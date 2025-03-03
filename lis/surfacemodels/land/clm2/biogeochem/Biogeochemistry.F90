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

subroutine Biogeochemistry (clm)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! surface biogeochemical fluxes
! 
! Method: 
! 
! Author: Gordon Bonan and Sam Levis
! 
!-----------------------------------------------------------------------
! $Id: Biogeochemistry.F90,v 1.6 2004/11/24 22:56:12 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
#if (defined DGVM)
  use clm2_shr_const_mod, only: SHR_CONST_CDAY
  use clm2_varcon  , only : tfrz
  use clm2_varpar  , only : nlevsoi
  use pft_varcon  , only : pftpar
  use LIS_timeMgrMod, only : get_step_size
#endif
  use histFileMod, only : histslf, histmlf 
  use DustEmissionMod, only : Dust
  implicit none

! ------------------------ arguments ------------------------------
  type (clm1d), intent(inout) :: clm       !CLM 1-D Module
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  integer  i,j           !indices
  real(r8) tf            !temperature factor
  real(r8) dmcf          !co2-to-biomass conversion (ug biomass / umol CO2)
#if (defined DGVM)
  real(r8), parameter :: k = 0.0548 / SHR_CONST_CDAY !from [/day] to [/second]
  real(r8) respcoeff     !respiration coefficient (LPJ)
  real(r8) l_cton        !c/n for leaves (LPJ)
  real(r8) s_cton        !c/n for stems (LPJ)
  real(r8) r_cton        !c/n for roots (LPJ)
  real(r8) tsoi,dep
#endif
! -----------------------------------------------------------------
  print*,' Biogeochem not supposed to be called...'
#if 0 
  dmcf = 28.5

! determine vegetation type

  i = clm%itypveg

! total photosynthesis

  clm%fpsn = clm%psnsun*clm%laisun + clm%psnsha*clm%laisha

#if (!defined DGVM)

  clm%frmf = 0.
  clm%frms = 0.
  clm%frmr = 0.

#else

! maintenance respiration: LPJ equations w/ units of [gC m-2 gridcell s-1]
! converted to LSM units of [umol CO2 m-2 patch s-1]

! --------------------------------------------------------------------
! Soil temperature to a depth of 0.25 m.
! --------------------------------------------------------------------

  tsoi = 0.
  dep = 0.
  do j = 1, nlevsoi
     if (clm%z(j)+0.5*clm%dz(j) <= 0.25) then
        tsoi = tsoi + clm%t_soisno(j)*clm%dz(j)
        dep = dep + clm%dz(j)
     end if
  end do
  if (dep /= 0.) then
     clm%tsoi25 = tsoi/dep
  else
     clm%tsoi25 = clm%t_soisno(1)
  end if

  respcoeff = pftpar(i,5)
  l_cton    = pftpar(i,13)
  s_cton    = pftpar(i,14)
  r_cton    = pftpar(i,15)

  if (i > 0 .and. clm%fpcgrid > 0.0) then

     if (clm%t_veg >= tfrz-40.) then
        tf = exp(308.56 * (1.0/56.02 - 1.0/(clm%t_veg-227.13)))
     else
        tf = 0.0
     end if

     clm%frmf = respcoeff * k * clm%lm_ind * clm%nind / l_cton * tf * clm%dphen * 2.0e6 / dmcf / clm%fpcgrid
     clm%frms = respcoeff * k * clm%sm_ind * clm%nind / s_cton * tf * &
  2.0e6 / dmcf / clm%fpcgrid

     if (clm%tsoi25 >= tfrz-40.) then
        tf = exp(308.56 * (1.0/56.02 - 1.0/(clm%tsoi25-227.13)))
     else
        tf = 0.0
     end if

     clm%frmr = respcoeff * k * clm%rm_ind * clm%nind / r_cton * tf * clm%dphen * 2.0e6 / dmcf / clm%fpcgrid

  else
     clm%frmf = 0.0
     clm%frms = 0.0
     clm%frmr = 0.0
  end if
#endif

  clm%frm  = clm%frmf + clm%frms + clm%frmr          

! growth respiration and production

#if (defined DGVM)
  clm%frg = 0.25 * max(clm%fpsn - clm%frm, 0.0)      !changed to match LPJ
  clm%dmi = (clm%fpsn - clm%frm - clm%frg) * dmcf
#else
  clm%frg = 0.                     
  clm%dmi = 0.     
#endif

#if (defined DGVM)
  clm%bm_inc = clm%bm_inc + clm%dmi * get_step_size() * 0.5e-6 !bm_inc=[gC/m2 patch area] from dmi=[ug dry matter/m2 patch area/s]

! microbial respiration

! DGVM calculates clm%fmicr in LitterSOM; problem with units in relation
! to history grid averages; {fmicr}=[gC/m2 gridcell vegetated area] in
! LPJ calculation => sum(fmicr) over a gridcell would give the total for
! the gridcell; it would be wrong to convert to [/m2 patch area] because
! soil carbon continues to exist even when the plants above it die and
! fpcgrid goes to 0; how to reconcile with history calculation which
! will take the following values and weight them by the corresponding
! weights of their patches? Could chg. soil carbon and plant litter to
! gridcell level pools; this will affect fco2, as well; could treat this
! issue simultaneously with soil water, which we also want converted to
! the gridcell level for plant competition. (slevis)
  clm%afmicr = clm%afmicr + clm%fmicr      ![gC/m2 gridcell vegetated area]
#else
  clm%fmicr = 0.
#endif

! net CO2 flux

!  clm%fco2 = -clm%fpsn + clm%frm + clm%frg + clm%fmicr

! call dust mobilization routine (C. Zender's modified codes)

  call Dust(clm)

! call VOC emission routine (A. Gunther's model)

  call VOCEmission(clm)
#endif
  return
end subroutine Biogeochemistry


