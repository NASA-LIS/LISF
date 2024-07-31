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

subroutine EcosystemDyn (clm, doalb, endofyr)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Ecosystem dynamics: phenology, vegetation
! 
! Method: 
! Calling sequence:
!    EcosystemDynamics:             !ecosystem dynamics driver
! This subroutine calculates:
!    o leaf areas (tlai, elai)
!    o stem areas (tsai, esai)
!    o height (htop)
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: EcosystemDyn.F90,v 1.7 2004/11/24 23:24:17 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
#if (defined DGVM)
  use clm2_shr_const_mod, only : SHR_CONST_CDAY
  use HistFileMod ,  only : histslf, histmlf
  use pft_varcon  ,  only : sla
#endif
  use LIS_timeMgrMod,  only : LIS_get_curr_date
  use LIS_coreMod, only : LIS_rc
  implicit none

! ------------------------ arguments ------------------------------
  type (clm1d), intent(inout) :: clm   !CLM 1-D Module
  logical     , intent(in)    :: doalb !true = surface albedo calculation time step
  logical     , intent(in)    :: endofyr
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  real(r8) :: ol      !thickness of canopy layer covered by snow (m)
  real(r8) :: fb      !fraction of canopy layer covered by snow
#if (defined DGVM)
  integer  :: yr      !year (0 -> ...)
  integer  :: mon     !month (1 -> 12)
  integer  :: day     !day (1 -> 31)
  integer  :: ncsec   !seconds of current date
#endif
! -----------------------------------------------------------------

#if (defined DGVM)

! today's phenology: returns dphen
! with small chg in subr, could be called every tstep if preferred

  if (clm%nstep==0 .or. &
       mod(clm%nstep-1, nint(SHR_CONST_CDAY/clm%dtime)) == 0. .or. endofyr) then
     call Phenology (clm)
  end if

! fire season; returns firelength for use at end of yr in subr. Fire
! litter and soil decomposition; returns litterag,bg and fmicr

  if (clm%nstep > 1) then
     call FireSeason (clm)
     call get_curr_date (LIS_rc%t, yr, mon, day, ncsec) 
     call LitterSOM (clm, yr)
  end if

! obtain vegetation structure here (tlai+sai,hbot+top)

  if (doalb .or. endofyr) then

     if (clm%itypveg > 0 .and. clm%fpcgrid > 0.0) then
        if (clm%itypveg <= 11) then          !woody vegetation (following IBIS)
           clm%tlai = clm%lai_ind * clm%dphen
           clm%tsai = 0.250 * clm%lai_ind
           clm%hbot = max(0.0, min(3.00, clm%htop-1.00))
        else                             !grasses (following IBIS)
           clm%tlai = clm%lai_ind * clm%dphen / clm%fpcgrid
           clm%tsai = 0.050 * clm%tlai
           clm%htop = max(0.25, clm%tlai * 0.25)
           clm%hbot = max(0.0, min(0.05, clm%htop-0.20))
        end if
     else
        clm%tlai = 0.0
        clm%tsai = 0.0
        clm%htop = 0.0
        clm%hbot = 0.0
     end if

! adjust lai and sai for burying by snow. if exposed lai and sai
! are less than 0.05 but non zero, set to 0.05 to prevent numerical
! problems associated with very small lai and sai.

     ol = min( max(clm%snowdp-clm%hbot,0._r4), clm%htop-clm%hbot)
     fb = 1. - ol / max(1.e-06,clm%htop-clm%hbot)
     clm%elai = clm%tlai*fb
     clm%esai = clm%tsai*fb
     if (clm%elai > 0.0 .and. clm%elai < 0.05) clm%elai = 0.05
     if (clm%esai > 0.0 .and. clm%esai < 0.05) clm%esai = 0.05

! Fraction of vegetation free of snow

     if ((clm%elai + clm%esai) >= 0.05) then
        clm%frac_veg_nosno_alb = 1
     else
        clm%frac_veg_nosno_alb = 0
     endif
  
  end if

#else

  if (doalb) then

! need to update elai and esai only every albedo time step so do not 
! have any inconsistency in lai and sai between SurfaceAlbedo calls (i.e., 
! if albedos are not done every time step).

! leaf phenology
! Set leaf and stem areas based on day of year
! Interpolate leaf area index, stem area index, and vegetation heights
! between two monthly 
! The weights below (timwt(1) and timwt(2)) were obtained by a call to 
! routine InterpMonthlyVeg in subroutine NCARlsm. 
!                 Field   Monthly Values
!                -------------------------
! leaf area index LAI  <- mlai1 and mlai2
! leaf area index SAI  <- msai1 and msai2
! top height      HTOP <- mhvt1 and mhvt2
! bottom height   HBOT <- mhvb1 and mhvb2

!     clm%tlai = timwt(1)*mlai1(clm%kpatch) + timwt(2)*mlai2(clm%kpatch)
!     clm%tsai = timwt(1)*msai1(clm%kpatch) + timwt(2)*msai2(clm%kpatch)
!     clm%htop = timwt(1)*mhvt1(clm%kpatch) + timwt(2)*mhvt2(clm%kpatch)
!     clm%hbot = timwt(1)*mhvb1(clm%kpatch) + timwt(2)*mhvb2(clm%kpatch)

! adjust lai and sai for burying by snow. if exposed lai and sai
! are less than 0.05, set equal to zero to prevent numerical 
! problems associated with very small lai and sai.

     ol = min( max(clm%snowdp-clm%hbot,0._r4), clm%htop-clm%hbot)
     fb = 1. - ol / max(1.e-06, clm%htop-clm%hbot)
     clm%elai = clm%tlai*fb
     clm%esai = clm%tsai*fb
     if (clm%elai < 0.05) clm%elai = 0._r4
     if (clm%esai < 0.05) clm%esai = 0._r4

! Fraction of vegetation free of snow


     if ((clm%elai + clm%esai) >= 0.05) then
        clm%frac_veg_nosno_alb = 1
     else
        clm%frac_veg_nosno_alb = 0
     endif

     if (clm%itypveg .eq.LIS_rc%bareclass &
          .and. (clm%elai + clm%esai) >= 0.05) then
       clm%frac_veg_nosno_alb = 0
     endif

     if (clm%itypveg.eq.0) & 
          clm%frac_veg_nosno_alb = 0 
!      write(*,*) clm%frac_veg_nosno_alb,clm%elai, clm%esai
  
  endif

#endif

  return
end subroutine EcosystemDyn
