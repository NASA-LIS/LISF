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

subroutine VOCEmission(clm)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Volatile organic compound emission
! 
! Method: 
! This code simulates volatile organic compound emissions
! following the algorithm presented in Gunther, A., 1999: Modeling
! Biogenic Volatile Organic Compound Emissions to the Atmosphere. In
! Reactive Hydrocarbons in the Atmosphere, Ch. 3
!
! This model relies on the assumption that 90% of isoprene and monoterpene
! emissions originate from canopy foliage:
!
! E = epsilon * gamma * density * delta
!
! The factor delta (longterm activity factor) applies to isoprene emission
! from deciduous plants only. We neglect this factor at the present time.
! This factor is discussed in Gunther (1997).
!
! Subroutine written to operate at the patch level.
!
! IN FINAL IMPLEMENTATION, REMEMBER:
!
! 1. may wish to call this routine only as freq. as rad. calculations
! 2. may wish to place epsilon values directly in pft-physiology file
!
! Output:
!
! vocflx(begpatch:endpatch,nvoc)) !VOC flux [ug C m-2 h-1]
! 
! Author: Sam Levis
! 
!-----------------------------------------------------------------------
! $Id: VOCEmission.F90,v 1.6 2004/11/24 22:56:14 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varpar, only : numpft	
  use pft_varcon, only : sla	
  use clm2_shr_const_mod, only : SHR_CONST_RGAS
  implicit none

! ----------------------- arguments -------------------------------
  type (clm1d), intent(inout) :: clm       !CLM 1-D Module
! -----------------------------------------------------------------

! ----------------------- local variables -------------------------
  integer n        !loop index

  real(r8) epsilon !emission factor [ug g-1 h-1]
  real(r8) gamma   !activity factor (instantaneous light and temp. condns)
  real(r8) density !source density factor [g dry wgt foliar mass/m2 ground]

  real(r8) cl
  real(r8) ct
  real(r8) par
  real(r8) reciprod

! constants

  real(r8), parameter :: alpha = 0.0027 !empirical coefficient
  real(r8), parameter :: cl1 = 1.066    !empirical coefficient
  real(r8), parameter :: ct1 = 95000.0  !empirical coefficient [J mol-1]
  real(r8), parameter :: ct2 = 230000.0 !empirical coefficient [J mol-1]
  real(r8), parameter :: ct3 = 0.961    !empirical coefficient
  real(r8), parameter :: tm  = 314.0    !empirical coefficient [K]
  real(r8), parameter :: R   = SHR_CONST_RGAS*0.001 !univ. gas constant [J K-1 mol-1]
  real(r8), parameter :: tstd = 303.0   !std temperature [K]
  real(r8), parameter :: bet = 0.09     !beta empirical coefficient [K-1]

! Specific leaf areas [m2 leaf g-1 carbon]
! These are the values from my version of genesis-ibis / 1000.
! With DGVM defined, use LPJ's sla [m2 leaf g-1 carbon]
! Divide by 2 in the equation to get dry weight foliar mass from grams carbon

#if (!defined DGVM)
  sla( 0) = 0.
  sla( 1) = 0.0125 !needleleaf
  sla( 2) = 0.0125 !Gordon Bonan suggests NET = 0.0076
  sla( 3) = 0.0125 !Gordon Bonan suggests NDT = 0.0200
  sla( 4) = 0.0250 !broadleaf
  sla( 5) = 0.0250 !Gordon Bonan suggests BET = 0.0178
  sla( 6) = 0.0250 !Gordon Bonan suggests BDT = 0.0274
  sla( 7) = 0.0250
  sla( 8) = 0.0250
  sla( 9) = 0.0250
  sla(10) = 0.0250
  sla(11) = 0.0250
  sla(12) = 0.0200 !grass
  sla(13) = 0.0200
  sla(14) = 0.0200
  sla(15) = 0.0200
  sla(16) = 0.0200 !numpft = 16
#endif

! -----------------------------------------------------------------

  do n = 1, nvoc

! epsilon: use values from table 3 in Gunther (1997) which originate in
! -------  Gunther et al. (1995). In the comments below, I mention the pft
!          category as described in table 3. Some values were taken directly
!          from Gunther et al. (1995). Units: [ug g-1 h-1]

! isoprenes:

     if (n == 1) then
        if (clm%itypveg >= 1 .and. clm%itypveg <= 3) then        !coniferous pfts
           epsilon = 18.
        else if (clm%itypveg == 4) then                      !trop. rain forest
           epsilon = 24.
        else if (clm%itypveg >= 5 .and. clm%itypveg <= 11) then  !deciduous pfts
           epsilon = 45.
        else if (clm%itypveg == 12) then                     !alpine grassland
           epsilon = 16.
        else if (clm%itypveg >= 13 .and. clm%itypveg <= 14) then !other grassland
           epsilon = 5.
        else if (clm%itypveg >= 15 .and. clm%itypveg <= 17) then !crops
           epsilon = 5.
        else
           epsilon = 0.
        end if

! monoterpenes:

     else if (n == 2) then
        if (clm%itypveg >= 1 .and. clm%itypveg <= 3) then        !coniferous pfts
           epsilon = 2.3
        else if (clm%itypveg == 4) then                      !trop. rain forest
           epsilon = 0.4
        else if (clm%itypveg >= 5 .and. clm%itypveg <= 11) then  !deciduous pfts
           epsilon = 0.8
        else if (clm%itypveg == 12) then                     !alpine grassland
           epsilon = 0.8
        else if (clm%itypveg >= 13 .and. clm%itypveg <= 14) then !other grassland
           epsilon = 0.2
        else if (clm%itypveg >= 15 .and. clm%itypveg <= 17) then !crops
           epsilon = 0.2
        else
           epsilon = 0.
        end if

! other VOCs (OVOCs)

     else if (n == 3) then
        epsilon = 1.5                 !Gunther (1995)

! other reactive VOCs (ORVOCs)

     else if (n == 4) then
        epsilon = 1.5                 !Gunther (1995)

     end if

! gamma: Activity factor. Units [dimensionless]
! -----

! isoprenes:

     if (n == 1) then

        par = clm%forc_solad(1) + clm%forc_solai(1)
        cl = alpha * cl1 * par * (1. + alpha * alpha * par * par)**(-0.5)

        reciprod = 1. / (R * clm%t_veg * tstd) !using t_veg for tleaf; ok?
        ct = exp(ct1 * (clm%t_veg - tstd) * reciprod) / &
             (ct3 + exp(ct2 * (clm%t_veg - tm) * reciprod))

        gamma = cl * ct !gamma = 1 under std temp & light condns

! monoterpenes, OVOCs, and ORVOCs (Gunther, 1999 and 1995):

     else

        gamma = exp(bet * (clm%t_veg - tstd)) !using t_veg for tleaf; ok?

     end if

! density: Source density factor [g dry weight foliar mass m-2 ground]
! -------

     if (clm%itypveg > 0) then
        density = clm%tlai / (sla(clm%itypveg) * 0.5)
     else
        density = 0.
     end if

! calculate the voc flux
! ----------------------

!     clm%vocflx(n) = epsilon * gamma * density

  end do

  return
end subroutine VOCEmission
