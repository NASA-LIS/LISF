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

subroutine Allocation (pftpar   , nind     , bm_inc   , &
                       tree     , height   , sla      , &
                       lm_ind   , sm_ind   , hm_ind   , &
                       rm_ind   , present  , crownarea, &
                       litter_ag, litter_bg, lai_ind  , &
                       fpc_grid , fpc_inc  , wscal    )

!----------------------------------------------------------------------- 
! 
! Purpose:
! 
! Method: Called once per year
! 
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. allocation)
! 
!-----------------------------------------------------------------------
! $Id: Allocation.F90,v 1.6 2004/11/24 22:56:52 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_shr_const_mod, ONLY: SHR_CONST_PI
  use clm2_varpar, ONLY: npftpar
  use pft_varcon, ONLY: reinickerp, allom1, allom2, allom3, latosa, wooddens
  implicit none

! ----------------------------- arguments ------------------------------
  real(r8), intent(out)   :: lai_ind
  real(r8), intent(out)   :: fpc_inc
  real(r8), intent(inout) :: fpc_grid
  real(r8), intent(inout) :: crownarea
  real(r8), intent(inout) :: height
  real(r8), intent(inout) :: lm_ind
  real(r8), intent(inout) :: sm_ind
  real(r8), intent(inout) :: hm_ind
  real(r8), intent(inout) :: rm_ind
  real(r8), intent(inout) :: litter_ag
  real(r8), intent(inout) :: litter_bg
  real(r8), intent(in)    :: bm_inc
  real(r8), intent(in)    :: nind
  real(r8), intent(in)    :: pftpar(npftpar)
  real(r8), intent(in)    :: sla
  logical , intent(in)    :: present
  logical , intent(in)    :: tree
  real(r8), intent(in)    :: wscal
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
  integer , parameter :: nseg = 20
  real(r8), parameter :: xacc = 0.1     !threshold x-axis and threshold
  real(r8), parameter :: yacc = 1.0e-10 !y-axis precision of allocation soln
  real(r8) :: bm_inc_ind
  real(r8) :: lmtorm
  real(r8) :: crownarea_max
  real(r8) :: lminc_ind_min
  real(r8) :: rminc_ind_min
  real(r8) :: lminc_ind
  real(r8) :: rminc_ind
  real(r8) :: sminc_ind
  real(r8) :: fpc_ind
  real(r8) :: fpc_grid_old
  real(r8) :: x1, x2, dx
  real(r8) :: sap_xsa
  real(r8) :: stemdiam
  real(r8) :: fx1
  real(r8) :: fmid
  real(r8) :: xmid
  real(r8) :: sign
  real(r8) :: rtbis
  real(r8) :: pi
! ----------------------------------------------------------------------

  pi = SHR_CONST_PI

  if (present) then

     bm_inc_ind = bm_inc / nind

! calculate this year's leaf to fine root mass ratio from mean annual
! water scalar and pft specific parameter
! slevis: in lpj wscal=awscal(pft)/aleafdays(pft), awscal=SUM(dwscal) (dphen>0)
!         dwscal=min(supply(pft)/demandpot(pft),1) or =1 when demand(pft)=0 etc
!         here wscal=annpsn/annpsnpot

     if (tree) then
        lmtorm = pftpar(18) * wscal

        crownarea_max = pftpar(20)

! TREE ALLOCATION

! Allocation of this year's biomass increment (bm_inc_ind) to the
! three living carbon pools, such that the basic allometric
! relationships (A-C below) are always satisfied.

! (A) (leaf area) = latosa * (sapwood xs area)
!       (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
! (B) (leaf mass) = lmtorm * (root mass)
! (C) height = allom2 * (stem diameter)**allom3 (source?)
! (D) (crown area) = min (allom1 * (stem diameter)**reinickerp, crownarea_max)

! Mathematical derivation:

!   (1) bm_inc_ind = lminc_ind + sminc_ind + rminc_ind
!   (2) leaf_area_new = latosa * sap_xsa_new   [from (A)]
!   (3) leaf_area_new = (lm_ind + lminc_ind) * sla
! from (2) & (3),
!   (4) (lm_ind + lminc_ind) * sla = latosa * sap_xsa_new
! from (4),
!   (5) sap_xsa_new = (lm_ind + lminc_ind) * sla / latosa
!   (6) (lm_ind + lminc_ind) = lmtorm * (rm_ind + rminc_ind)
!         [from (B)]
!   (7) height_new = allom2 * stemdiam_new**allom3  [from (C)]
! from (1),
!   (8) sminc_ind = bm_inc_ind - lminc_ind - rminc_ind
! from (6),
!   (9) rminc_ind=((lm_ind + lminc_ind) / lmtorm) - rm_ind
! from (8) & (9),
!  (10) sminc_ind = bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind)  / lmtorm) + rm_ind
!  (11) wooddens = (sm_ind + sminc_ind + hm_ind) / stemvolume_new
!  (12) stemvolume_new = height_new * pi * stemdiam_new**2 / 4
! from (10), (11) & (12)
!  (13) stemdiam_new = [ ((sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
!         / wooddens) / (height_new * pi / 4) ]**(1/2)
! combining (7) and (13),
!  (14) height_new = allom2 * [ ((sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
!         / wooddens) / (height_new * pi / 4) ]**(1/2 * allom3)
! from (14),
!  (15) height_new**(1 + 2 / allom3) = allom2**(2 / allom3)
!         * ((sm_ind + bm_inc_ind - lminc_ind - ((lm_ind + lminc_ind)
!         / lmtorm) + rm_ind + hm_ind) / wooddens) / (pi / 4)
!  (16) wooddens = (sm_ind + sminc_ind) / sapvolume_new
! from (10) and (16),
!  (17) wooddens = (sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind) / sapvolume_new
!  (18) sapvolume_new = height_new * sap_xsa_new
! from (17) and (18),
!  (19) sap_xsa_new = (sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind)
!         / (height_new * wooddens)
! from (19),
!  (20) height_new = (sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
!         / (sap_xsa_new * wooddens)
! from (5) and (20),
!  (21) height_new**(1 + 2 / allom3) = [ (sm_ind + bm_inc_ind
!         - lminc_ind - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
!         / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
!         **(1 + 2 / allom3)
! -------------------------------------------------------------------
!  (15) and (21) are two alternative expressions for
!       height_new**(1 + 2 / allom3). Combining these,

!  (22) allom2**(2 / allom3) * ((sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
!         / wooddens) / (pi / 4) - [ (sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
!         / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
!         **(1 + 2 / allom3)
!         = 0

! Equation (22) can be expressed in the form f(lminc_ind)=0.

! Numerical methods are used to solve the equation for the
! unknown lminc_ind.
! -------------------------------------------------------------------

! Work out minimum leaf production to maintain current sapmass

!  (23) sap_xsa = sm_ind / wooddens / height
! from (A) and (23),
!  (24) leaf_mass * sla = latosa * sap_mass / wooddens / height
! from (24),
!  (25) leaf_mass = latosa * sap_mass / (wooddens * height * sla)
! from (25), assuming sminc_ind=0,
!  (26) lm_ind + lminc_ind_min = latosa * sm_ind
!         / (wooddens * height * sla)
! from (26),
!  (27) lminc_ind_min = latosa * sm_ind / (wooddens * height * sla)
!         - lm_ind

        lminc_ind_min = latosa*sm_ind/(wooddens*height*sla) - lm_ind !eqn (27)

! Work out minimum root production to support this leaf mass
! (i.e. lm_ind + lminc_ind_min)
! May be negative following a reduction in soil water limitation
! (increase in lmtorm) relative to last year.

! from (B) and (25),
!  (28) root_mass = latosa * sap_mass / (wooddens * height * sla)
!         / lmtorm
! from (28), assuming sminc_ind=0,
!  (29) rm_ind + rminc_ind_min = latosa * sm_ind
!         / (wooddens * height * sla * lmtorm)
! from (29),
!  (30) rminc_ind_min = latosa * sm_ind
!         / (wooddens * height * sla * lmtorm) - rm_ind

        rminc_ind_min = latosa*sm_ind/(wooddens*height*sla*lmtorm) - rm_ind !eqn(30)

        if (rminc_ind_min > 0.0 .and. lminc_ind_min > 0.0 .and. &
            rminc_ind_min + lminc_ind_min <= bm_inc_ind          ) then

! Normal allocation (positive increment to all living C compartments)

! Calculation of leaf mass increment (lminc_ind) that satisfies
! Eqn (22) using Bisection Method (Press et al 1986, p 346)

! Seeking a root for non-negative lminc_ind, rminc_ind and
! sminc_ind.  There should be exactly one (no proof presented, but
! Steve has managed one) and it should lie between x1=0 and
! x2=(bm_inc_ind-(lm_ind/lmtorm-rm_ind))/(1+1/lmtorm).

           x1 = 0.0
           x2 = (bm_inc_ind - (lm_ind/lmtorm-rm_ind)) / (1.0 + 1.0 / lmtorm)
           dx = (x2-x1) / real(nseg)

           if (lm_ind == 0.0) x1 = x1 + dx !to avoid division by zero

! evaluate f(x1)=LHS of eqn (22) at x1

           fx1 = allom2**(2.0/allom3) *                                       &
        (sm_ind + bm_inc_ind - x1 - (lm_ind + x1)/lmtorm + rm_ind + hm_ind) / &
        wooddens / (pi/4.0) -                                                 &
        ( (sm_ind + bm_inc_ind - x1 - (lm_ind + x1)/lmtorm + rm_ind) /        &
          ((lm_ind + x1)*sla*wooddens/latosa) )**(1.0+2.0/allom3)

! Find approximate location of leftmost root on the interval
! (x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
! change in sign of f(xmid) relative to f(x1).

           fmid = fx1
           xmid = x1

           do while (fmid*fx1 > 0.0 .and. xmid < x2)

              xmid = xmid + dx
              fmid = allom2**(2.0/allom3) *                                   &
      (sm_ind + bm_inc_ind - xmid - (lm_ind+xmid)/lmtorm + rm_ind + hm_ind) / &
      wooddens / (pi/4.0) -                                                   &
      ( (sm_ind + bm_inc_ind - xmid - (lm_ind+xmid)/lmtorm + rm_ind) /        &
        ((lm_ind + xmid)*sla*wooddens/latosa) )**(1.0+2.0/allom3)

           enddo

           x1 = xmid - dx
           x2 = xmid

! Apply bisection method to find root on the new interval (x1,x2)

           fx1 = allom2**(2.0/allom3) *                                       &
        (sm_ind + bm_inc_ind - x1 - (lm_ind + x1)/lmtorm + rm_ind + hm_ind) / &
        wooddens / (pi/4.0) -                                                 &
        ( (sm_ind + bm_inc_ind - x1 - (lm_ind + x1)/lmtorm + rm_ind) /        &
          ((lm_ind + x1)*sla*wooddens/latosa) )**(1.0+2.0/allom3)

           if (fx1 >= 0.0) then
              sign=-1.0
           else
              sign=1.0
           endif

           rtbis=x1
           dx=x2-x1

! Bisection loop
! Search iterates on value of xmid until xmid lies within
! xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0

           fmid=1.0  !dummy value to guarantee entry to loop

           do while (dx >= xacc .and. abs(fmid) > yacc)

              dx=dx*0.5
              xmid=rtbis+dx

! calculate fmid=f(xmid) [eqn (22)]

              fmid = allom2**(2.0/allom3) *                                   &
      (sm_ind + bm_inc_ind - xmid - (lm_ind+xmid)/lmtorm + rm_ind + hm_ind) / &
      wooddens / (pi/4.0) -                                                   &
      ( (sm_ind + bm_inc_ind - xmid - (lm_ind+xmid)/lmtorm + rm_ind) /        &
        ((lm_ind + xmid)*sla*wooddens/latosa) )**(1.0+2.0/allom3)

              if (fmid*sign <= 0.0) rtbis=xmid

           enddo

! Now rtbis contains numerical solution for lminc_ind given eqn (22)

           lminc_ind=rtbis

! Calculate increments in other compartments using allometry relationships

           rminc_ind = (lm_ind + lminc_ind) / lmtorm - rm_ind !eqn (9)
           sminc_ind = bm_inc_ind - lminc_ind - rminc_ind     !eqn (1)

        else

! Abnormal allocation: reduction in some C compartment(s) to satisfy allometry

! Attempt to distribute this year's production among leaves and roots only
!  (31) bm_inc_ind = lminc_ind + rminc_ind
! from (31) and (9),
!  (32) bm_inc_ind = lminc_ind + ((lm_ind + lminc_ind) / lmtorm)
!         - rm_ind
! from (32)
!  (33) lminc_ind = (bm_inc_ind - lm_ind / lmtorm + rm_ind) /
!         (1 + 1 / lmtorm)

           lminc_ind = (bm_inc_ind - lm_ind/lmtorm + rm_ind) / &
                       (1.0 + 1.0/lmtorm) !eqn (33)

           if (lminc_ind >= 0.0) then

! Positive allocation to leafmass

              rminc_ind = bm_inc_ind - lminc_ind  !eqn (31)

! Add killed roots (if any) to below-ground litter

              if (rminc_ind < 0.0) then
                 lminc_ind = bm_inc_ind
                 rminc_ind = (lm_ind + lminc_ind) / lmtorm - rm_ind
                 litter_bg = litter_bg - rminc_ind * nind
              endif
           else

! Negative allocation to leaf mass

              rminc_ind = bm_inc_ind
              lminc_ind = (rm_ind + rminc_ind) * lmtorm - lm_ind !from eqn (9)

! Add killed leaves to litter

              litter_ag = litter_ag - lminc_ind * nind

           endif

! Calculate sminc_ind (must be negative)

! from (25),
!  (34) lm_ind + lminc_ind = latosa * (sm_ind + sminc_ind)
!         / (wooddens * height * sla)
! from (34),
!  (35) sminc_ind = (lm_ind + lminc_ind) * wooddens * height * sla
!         / latosa - sm_ind

           sminc_ind = (lm_ind + lminc_ind)*wooddens*height*sla / latosa - &
                       sm_ind !eqn (35)

! Convert killed sapwood to heartwood

           hm_ind = hm_ind - sminc_ind

        endif

! Increment C compartments

        lm_ind = lm_ind + lminc_ind
        rm_ind = rm_ind + rminc_ind
        sm_ind = sm_ind + sminc_ind

! Calculate new height, diameter and crown area

        sap_xsa = lm_ind * sla / latosa  !eqn (5)
        height = sm_ind / sap_xsa / wooddens
        stemdiam = (height/allom2)**(1.0/allom3) !eqn (C)
        crownarea = min(allom1*stemdiam**reinickerp, crownarea_max) !eqn (D)

     else !grasses
        lmtorm = pftpar(18) !slevis: eliminate influence from soil H2O

! GRASS ALLOCATION
! Distribute this year's production among leaves and fine roots
! according to leaf to rootmass ratio [eqn (33)]
! Relocation of C from one compartment to the other not allowed:
! negative increment in either compartment transferred to litter

        lminc_ind = (bm_inc_ind - lm_ind/lmtorm + rm_ind) / (1.0 + 1.0/lmtorm)
        rminc_ind = bm_inc_ind - lminc_ind

        if (lminc_ind >= 0.0) then

! Add killed roots (if any) to below-ground litter

! CHECK: take out if statement because if rminc is negative than
! root mass has been translocated to the leaves, therefore mass balance
! problem since this carbon stays in the vegetation but is in addition
! added to the litter pool. ALLOW translocation from roots to leaves
! i.e. assume carbon stores in the roots which can be delivered
! to the leaves under times of stress.

!          if (rminc_ind < 0.0) litter_bg = litter_bg -rminc_ind * nind

        else

! Negative allocation to leaf mass

           rminc_ind = bm_inc_ind
           lminc_ind = (rm_ind + rminc_ind)*lmtorm - lm_ind !from eqn (9)

! Add killed leaves to litter

           litter_ag = litter_ag - lminc_ind * nind

        endif

! Increment C compartments

        lm_ind = lm_ind + lminc_ind
        rm_ind = rm_ind + rminc_ind

     endif

! Update LAI and FPC

     if (crownarea > 0.0) then
        lai_ind = lm_ind * sla / crownarea
     else
        lai_ind = 0.0
     endif

     fpc_ind = 1.0 - exp(-0.5*lai_ind)
     fpc_grid_old = fpc_grid
     fpc_grid = crownarea * nind * fpc_ind
     fpc_inc = max(0.0, fpc_grid - fpc_grid_old)

! diagnostic (slevis)
!    write(15,*)lminc_ind/bm_inc_ind,rminc_ind/bm_inc_ind,sminc_ind/bm_inc_ind
! end diagnostic
  endif

  return
end subroutine Allocation
