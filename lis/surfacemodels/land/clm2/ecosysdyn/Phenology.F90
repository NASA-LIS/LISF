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

subroutine Phenology (clm)

!----------------------------------------------------------------------- 
! 
! Purpose: Summer and drought phenology
! 
! Method: Called once per day
! 
! Author: Sam Levis (adapted from Jon Foley's IBIS subroutine pheno)
! 
!-----------------------------------------------------------------------
! $Id: Phenology.F90,v 1.6 2004/11/24 22:57:00 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, ONLY : tfrz
  use pft_varcon, ONLY : tree, summergreen, raingreen, pftpar
  implicit none

! ------------------------ arguments ------------------------------
  type (clm1d), intent(inout) :: clm   !CLM 1-D Module
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  real(r8), parameter :: ddfacu = 1.0 / 15.0 !'drop day factor' causes pheno-
  real(r8), parameter :: ddfacl = 1.0 /  5.0 !logy to switch in 15 or 5 days
  real(r8)            :: tthreshold
! -----------------------------------------------------------------
  print*,'shouldnt be in phenology..'
#if 0 
  clm%t10min = min (clm%t10min, clm%t10) !reset to 1.0e+36 once per year

  if (tree(clm%itypveg) .and. .not.summergreen(clm%itypveg) .and. .not.raingreen(clm%itypveg)) then
     clm%dphen = 1.0
  else if (tree(clm%itypveg) .and. summergreen(clm%itypveg)) then
! ---------------------------------------------------------------------
! * * * upper canopy winter phenology * * *
! ---------------------------------------------------------------------

! temperature threshold for budburst and senescence

! temperature threshold is assumed to be tfrz (= 273.16 K)
! or 5 degrees warmer than the coldest monthly temperature
! slevis: change that to be coldest 10-day running mean temperature

     tthreshold = max (tfrz, clm%t10min + 5.0)

! determine if growing degree days are initiated
! slevis:*t10 = 10-day running mean air temperature (K)

! determine leaf display

     if (clm%t10 < tthreshold) then
       clm%dphen = max (0.0, clm%dphen - ddfacu)
     else
       print*,'not supposed to be here.. agdd0 removed..'
!       clm%dphen = min (1.0, max (0.0, clm%agdd0 - 100.0) / 50.0)
     endif

  else if (clm%itypveg > 0 .and. .not.tree(clm%itypveg)) then !NB: grass has no specific phenology
! ---------------------------------------------------------------------
! * * * lower canopy phenology * * *
! ---------------------------------------------------------------------

! temperature threshold for budburst and senescence

! temperature threshold is assumed to be tfrz

    tthreshold = tfrz

! determine leaf display

    if      (clm%t10 < tthreshold) then        !cold phenology for grasses
      clm%dphen = max (0.0, clm%dphen - ddfacl)!slevis: made ddfacl=1/5
    else if (clm%fnpsn10 < 0.0) then           !drought phenology for grasses
      clm%dphen = max (0.1, clm%dphen - ddfacl)!slevis: made ddfacl=1/5
    else
!     clm%dphen = min (1.0, max (0.0, clm%agdd5 - 150.0) / 50.0)
      clm%dphen = min (1.0, clm%dphen + ddfacl)!slevis: try this line instead
    endif

  end if

  if (tree(clm%itypveg) .and. raingreen(clm%itypveg)) then
! ---------------------------------------------------------------------
! * * * upper canopy drought phenology * * *
! ---------------------------------------------------------------------

    if (clm%fnpsn10 <  0.0) clm%dphen = max (0.1, clm%dphen - ddfacu)
    if (clm%fnpsn10 >= 0.0) clm%dphen = min (1.0, clm%dphen + ddfacu)

! trying out enforced drought phenology

    if (clm%dphen > 0.95) clm%leafon = clm%leafon + 1.0
    if (clm%leafon >= 365.0*pftpar(clm%itypveg,10)) then
       clm%dphen = 0.1
       clm%leafof = clm%leafof + 1.0
       if (clm%leafof >= 365.0*pftpar(clm%itypveg,10)) then
          clm%dphen = 1.0
          clm%leafof = 0.0
          clm%leafon = 1.0
       end if
    end if

  end if
#endif
  return
end subroutine Phenology














