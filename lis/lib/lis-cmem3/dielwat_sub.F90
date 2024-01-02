!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!!!! Insert banner here !!!!
!SUBROUTINE DIEL_WAT 

! Purpose : 
!   Calculate dielectric constant of water in three different media : 
!   pure water, sea water, soil water

! Reference:
!  Dielectric constant of pure water
!   Ulaby p 2020
!  Dielectric constant of saline water
!   1) Stogryn, A. (1971): Equations for calculating the dielectric constant of
!    saline water, IEEE Transactions on Microwave Theory and Techniques,
!    Vol. MTT-19, 733-736.
!   2) Klein, L. A. and C. T. Swift (1977): An improved model
!    for the dielectric constant of sea water at microwave
!    frequencies, IEEE Transactions on  Antennas and Propagation,
!    Vol. AP-25, No. 1, 104-111.
!  Dielectric constant of soil water
!   1) Dobson '85. Modified Debye expression
!         Stern_Gouy double layer theory
!   2) Ulaby p 2024

! Input Variables:
!   fghz         frequency (GHz)
!   T 		 temperature of water (C)
!   sal 	 water salinity (psu = ppt(weight) )
!   wc           volumetric soil water content (cm3/cm3)
!   sand         clay fraction (% of dry weight of soil)
!   clay         clay fraction (% of dry weight of soil)
!   medium 	 pure water(0) sea water(1) soil water(2)
!   isal 	 Stogryn (1) Klein and Swift (2)

! Output variables:
!   ew           dielectric constant of water 

! local variables :
!   rho_s        rho_s : soil specific density (g/cm3), SMOS ATBD 2.664
!   rho_b        soil bulk density (g/cm3) ~1.3
!  N : normality from salinity (Stogryn, modified by Klein and Swift 1977)
!   ew  : dielectric constant of water
!  sal : water salinity (psu = ppt(weight) )
!  eps_w0 : static dielectric constant of pure water (Klein and Swift 1977) 
!  eps_sw0 : Static dielectric constant of soil water
!  tau_w : relaxation time of pure water (stogryn, 1970)
!  tau_sw : relaxation time of saline water
!  sigma : ionic conductivity
!  sigma_eff : effective conductivity of water (S/m)
!---------------------------------------------------------------------------

SUBROUTINE DIEL_WAT (fghz, T, sal, wc, sand, clay, medium, isal, ew)

IMPLICIT NONE

INTEGER :: medium, isal
real :: fghz, T, sal, wc, sand, clay
complex :: ew

real :: f, pi, eps_0, eps_winf, rho_b, rho_s
REAL :: N, omega
REAL :: sigma, sigma_eff
REAL :: tau_w, tau_sw
REAL :: eps_w0, eps_sw0, a, bb
COMPLEX :: j
!---------------------------------------------------------------------------

f = fghz * 1.0e9
pi = acos(-1.0)
eps_0 = 8.854e-12
eps_winf = 4.9
omega = 2.0 * pi * f

rho_s = 2.66
rho_b = (sand*1.6 +  clay*1.1 + (100.-sand-clay)*1.2)/100.

j = (0.,1.)

tau_w = 1.768e-11  - 6.068e-13 * T  + 1.104e-14 * T**2  - 8.111e-17 * T**3
! same as:
!tau_w = 1./(2.*pi) * (1.1109e-10 - 3.824e-12 * T + 6.938e-14 * T**2  &
!     - 5.096e-16 * T**3)

SELECT CASE (isal)
  
  CASE ( 1 )
    N = 0.9141 * sal * (1.707e-2 + 1.205e-5 * sal + 4.058e-9 * sal**2)
    
    eps_sw0 = 87.74  - 0.4008 * T  + 9.398e-4 * T**2  + 1.410e-6 * T**3
    a = 1.0  - 0.2551 * N  + 5.151e-2 * N**2  - 6.889e-3 * N**3
    eps_sw0 = eps_sw0 * a

    bb = 1.0 - 0.04896 * N - 0.02967 * N**2 + 5.644e-3 &
       & * N**3 + 0.1463e-2 * N * T
    tau_sw  = tau_w * bb
   
  CASE ( 2 )
  
    eps_sw0 = 87.134  - 1.949e-1 * T  - 1.276e-2 * T**2  + 2.491e-4 * T**3
    a = 1.000  + 1.613e-5 * sal * T  - 3.656e-3 * sal  + 3.210e-5 * sal**2  &
             -  4.232e-7 * sal**3
    eps_sw0 = eps_sw0 * a
    
    bb = 1.000  + 2.282e-5 * sal * T  - 7.638e-4 * sal  - 7.760e-6 * sal**2  &
               + 1.105e-8 * sal**3
    tau_sw  = tau_w * bb
   
END SELECT

SELECT CASE (medium)

  CASE ( 0 ) ! pure water
    eps_w0 = 88.045 - 0.4147 * T + 6.295e-4 * T**2 + 1.075e-5 * T**3 
    ew = eps_winf + (eps_w0 - eps_winf) / (1. - j * omega * tau_w)
    
  CASE ( 1 )
    CALL ION_CONDUCT (T,sal,sigma)
    ! Debye expression '41
    ew = eps_winf + (eps_sw0 - eps_winf) / (1. - j * omega * tau_sw)  &
         + j * sigma / (omega * eps_0)
  
  CASE ( 2 )
    sigma_eff = -1.645  + 1.939 * rho_b  - 0.02256 * sand  + 0.01594 * clay
    !  sigma_eff gets negative for very sandy soils
    !   with low bulk densities. If this happens set sigma_eff to zero.
    IF (sigma_eff < 0.) sigma_eff = 0.
    
    ! Modified Debye expression, Dobson '85
    !IF (wc < 0.001) wc = 0.001     ! to avoid dividing by zero
     wc = MAX(0.001, wc)  ! to avoid dividing by zero
    ew = eps_winf + (eps_sw0 - eps_winf) / (1. - j * omega * tau_sw)  &
         + j * sigma_eff / (omega * eps_0) * (rho_s - rho_b) / (rho_s * wc)

END SELECT

END SUBROUTINE DIEL_WAT

!===========================================================================
! Internal subroutines
!===========================================================================

SUBROUTINE ION_CONDUCT (T,sal,sigma)

! Purpose :
!  Calculate ionic conductivity for saline water
! Reference:
!  Stogryn, 1971

! local variables :
!  T : temperature of water (C)
!  sigma25 : ionic conductivity for sea water at 25 degree C (S/m)
!  sigma : ionic conductivity for saline water (S/m)
!---------------------------------------------------------------------------

IMPLICIT NONE

REAL :: T
REAL :: sal
REAL :: alphac
REAL :: delta
REAL :: sigma25
REAL :: sigma
!---------------------------------------------------------------------------

delta  =  25. - T

alphac  =  2.033e-2  +  1.266e-4 * delta  +  2.464e-6 * delta ** 2.0   &
          - sal * (1.849e-5  -  2.551e-7 * delta  +  2.551e-8 * delta ** 2.0)

sigma25  =  sal * (0.182521  -  1.46192e-3 * sal        &
                             +  2.09324e-5 * sal ** 2.0   &
                             -  1.28205e-7 * sal ** 3.0)

sigma   = sigma25 * exp(-1. * delta * alphac)

END SUBROUTINE ION_CONDUCT
