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

subroutine LitterSOM (clm, kyr)

!----------------------------------------------------------------------- 
! 
! Purpose: Litter and soil decomposition
! 
! Method: Incorporates analytical solution for soil pool sizes
!         once litter inputs are (assumed to be) at equilibrium,
!         reducing spin-up time for carbon fluxes due to soil respiration.
! 
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine littersom)
! 
!-----------------------------------------------------------------------
! $Id: LitterSOM.F90,v 1.6 2004/11/24 22:56:58 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : tfrz
  use clm2_shr_const_mod, only : SHR_CONST_CDAY
  implicit none

! ------------------------ arguments ------------------------------
  type (clm1d), intent(inout) :: clm   !CLM 1-D Module
  integer     , intent(in) :: kyr      !year (0, ...) for nstep+1
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  integer , parameter :: soil_equil_year = 400     !number of years until pool sizes for soil decomposition solved analytically
  real(r8), parameter :: k_litter10 = 0.5          !litter decomp. rate at 10 deg C (/year)
  real(r8), parameter :: k_soil_fast10 = 0.03      !fast pool decomp. rate at 10 deg C (/year)
  real(r8), parameter :: k_soil_slow10 = 0.001     !slow pool decomp. rate at 10 deg C (/year)
  real(r8), parameter :: atmfrac = 0.7             !fraction of litter decomp. going directly into the atmosphere
  real(r8), parameter :: soilfrac = 1.0 - atmfrac  !fraction of litter decomp. going to soil C pools
  real(r8), parameter :: fastfrac = 0.985          !fraction of litter entering fast soil decomposition pool
  real(r8), parameter :: slowfrac = 1.0 - fastfrac !fraction of litter entering slow soil decomposition pool
  real(r8) :: temp_resp          !temperature response of decomposition
  real(r8) :: moist_resp         !moisture response of decomposition
  real(r8) :: k_litter           !litter decomposition rate (/tstep)
  real(r8) :: k_fast             !fast pool decomposition rate (/tstep)
  real(r8) :: k_slow             !slow pool decomposition rate (/tstep)
  real(r8) :: litter_decom       !litter decomposition
  real(r8) :: litter_decom_ag    !above-ground litter decomposition
  real(r8) :: litter_decom_bg    !below-ground litter decomposition
  real(r8) :: cflux_litter_soil  !litter decomposition flux to soil
  real(r8) :: cflux_litter_atmos !litter decomposition flux to atmosphere
  real(r8) :: cflux_fast_atmos   !soil fast pool decomposition flux to atmos.
  real(r8) :: cflux_slow_atmos   !soil slow pool decomposition flux to atmos.
! real(r8) :: tsoi,dep
! integer  :: j
! -----------------------------------------------------------------

! Temperature response function is a modified Q10 relationship
! (Lloyd & Taylor 1994)
! slevis: Original code used monthly avg soil temp (K); I use tstep value
  print*,'LitterSOM not supposed to be called..'
#if 0 
  if (clm%tsoi25 <= tfrz - 40.0) then !avoid division by zero
     temp_resp=0.0
  else ! Lloyd & Taylor 1994
     temp_resp=exp(308.56*((1.0/56.02)-(1.0/(clm%tsoi25-227.13))))
  endif

! Moisture response based on soil layer 1 moisture content (Foley 1995)
! slevis: Orig. code used monthly soil water in upper 0.5 m (fraction of whc)
!         I use the tstep value

  moist_resp = 0.25 + 0.75 * clm%wf

! Original divided by 12 to get monthly decomposition rates (k, /month)
! as a function of temperature and moisture
! slevis: make rates /tstep by dividing by the number of tsteps per year

  k_litter = k_litter10    * temp_resp * moist_resp * clm%dtime / (SHR_CONST_CDAY * 365.)
  k_fast   = k_soil_fast10 * temp_resp * moist_resp * clm%dtime / (SHR_CONST_CDAY * 365.)
  k_slow   = k_soil_slow10 * temp_resp * moist_resp * clm%dtime / (SHR_CONST_CDAY * 365.)

! Calculate monthly litter decomposition using equation
!   (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate
! from (1),
!   (2) c = c0*exp(-kt) where c0=initial pool size
! from (2), decomposition in any month given by
!   (3) delta_c = c0 - c0*exp(-k)
! from (3)
!   (4) delta_c = c0*(1.0-exp(-k))

  litter_decom_ag = clm%litterag * (1.0-exp(-k_litter))  !eqn 4
  litter_decom_bg = clm%litterbg * (1.0-exp(-k_litter))
  litter_decom    = litter_decom_ag + litter_decom_bg

! Update the litter pools

  clm%litterag = clm%litterag - litter_decom_ag
  clm%litterbg = clm%litterbg - litter_decom_bg

! Calculate carbon flux to atmosphere and soil

  cflux_litter_atmos = atmfrac  * litter_decom
  cflux_litter_soil  = soilfrac * litter_decom

! Further subdivide soil fraction between fast and slow soil pools

  clm%cpool_fast = clm%cpool_fast + fastfrac * cflux_litter_soil
  clm%cpool_slow = clm%cpool_slow + slowfrac * cflux_litter_soil

! Calculate monthly soil decomposition to the atmosphere

  cflux_fast_atmos = clm%cpool_fast * (1.0-exp(-k_fast))  !eqn 4
  cflux_slow_atmos = clm%cpool_slow * (1.0-exp(-k_slow))  !eqn 4

! Update the soil pools

  clm%cpool_fast = clm%cpool_fast - cflux_fast_atmos
  clm%cpool_slow = clm%cpool_slow - cflux_slow_atmos

! Calculate heterotrophic respiration (in LSM referred to as microbial)

  clm%fmicr = cflux_litter_atmos + cflux_fast_atmos + cflux_slow_atmos

! Empty soil pools below a minimum threshold

  if (clm%cpool_fast < 1.0e-5) clm%cpool_fast = 0.0
  if (clm%cpool_slow < 1.0e-5) clm%cpool_slow = 0.0

  if (kyr <= soil_equil_year) then

! Update running average respiration rates and litter input
! slevis: had to multiply the denominator to chg units from years to tsteps

     clm%k_fast_ave       = clm%k_fast_ave       + k_fast / &
        (real(soil_equil_year) * 365. * SHR_CONST_CDAY / clm%dtime)
     clm%k_slow_ave       = clm%k_slow_ave       + k_slow / &
        (real(soil_equil_year) * 365. * SHR_CONST_CDAY / clm%dtime)
     clm%litter_decom_ave = clm%litter_decom_ave + litter_decom / &
        (real(soil_equil_year) * 365. * SHR_CONST_CDAY / clm%dtime)

! SOIL DECOMPOSITION EQUILIBRIUM CALCULATION
! Analytical solution of differential flux equations for fast and slow
! soil carbon pools.  Implemented after (soil_equil_year) simulation
! years, when annual litter inputs should be close to equilibrium.  Assumes
! average climate (temperature and soil moisture) from all years up to
! soil_equil_year.
! slevis: next could be done once

  else if (kyr == soil_equil_year+1) then

! Analytically calculate pool sizes this year only

! Rate of change of soil pool size = litter input - decomposition
!   (5) dc/dt = litter_decom - kc
! At equilibrium,
!   (6) dc/dt = 0
! From (5) & (6),
!   (7) c = litter_decom / k

    clm%cpool_fast=soilfrac*fastfrac*clm%litter_decom_ave/clm%k_fast_ave !eqn 7
    clm%cpool_slow=soilfrac*slowfrac*clm%litter_decom_ave/clm%k_slow_ave !eqn 7

  endif
#endif
  return
end subroutine LitterSOM
