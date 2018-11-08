! cable_constants.f90
!
! Source file containing constants for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach, 
!
! bugs to bernard.pak@csiro.au
!
! This file contains modules:
!  math_constants
!  physical_constants
!  other_constants
!  photosynthetic_constants
!  spatial_heterogeneity

MODULE cable_photosynthetic_constants
  USE cable_dimensions, ONLY : i_d, r_1
  IMPLICIT NONE
  PRIVATE i_d, r_1
  INTEGER(i_d), PARAMETER :: maxiter=20 ! max # interations for leaf temperature
  ! a1c3 is defined inside cable_canopy.f90
  REAL(r_1), PARAMETER :: a1c4 = 4.0
  REAL(r_1), PARAMETER :: a1c3_default = 9.0
  REAL(r_1), PARAMETER :: d0c3_hawkesbury = 5000.0 ! Empirical coef for vpd sensitivity of stomata
  REAL(r_1), PARAMETER :: d0c3_default = 1500.0 ! Empirical coef for vpd sensitivity of stomata
  REAL(r_1), PARAMETER :: d0c4 = 1500.0 ! Empirical coef for vpd sensitivity of stomata
  REAL(r_1), PARAMETER :: alpha3 = 0.200
  REAL(r_1), PARAMETER :: alpha4  = 0.05
  REAL(r_1), PARAMETER :: cfrd3  = 0.015
  REAL(r_1), PARAMETER :: cfrd4  = 0.025
  REAL(r_1), PARAMETER :: conkc0 = 302.0E-6  !mol mol^-1
  REAL(r_1), PARAMETER :: conko0 = 256.0E-3  !mol mol^-1
  REAL(r_1), PARAMETER :: convx3 = 1.0E-2
  REAL(r_1), PARAMETER :: convx4 = 0.8
  REAL(r_1), PARAMETER :: ekc = 59430.0  !J mol^-1
  REAL(r_1), PARAMETER :: eko = 36000.0  !J mol^-1
  REAL(r_1), PARAMETER :: gam0 = 28.0E-6  !mol mol^-1 @ 20C = 36.9 @ 25C
  REAL(r_1), PARAMETER :: gam1 = 0.0509
  REAL(r_1), PARAMETER :: gam2 = 0.0010
  REAL(r_1), PARAMETER :: gsw03  = 0.01
  REAL(r_1), PARAMETER :: gsw04  = 0.04
  REAL(r_1), PARAMETER :: rgbwc  = 1.32
  REAL(r_1), PARAMETER :: rgswc  = 1.57
  REAL(r_1), PARAMETER :: tmaxj  = 45.0
  REAL(r_1), PARAMETER :: tmaxv  = 45.0
  REAL(r_1), PARAMETER :: tminj  = -5.0
  REAL(r_1), PARAMETER :: tminv  = -5.0
  REAL(r_1), PARAMETER :: toptj  = 20.0
  REAL(r_1), PARAMETER :: toptv  = 20.0
  REAL(r_1), PARAMETER :: trefk= 298.2  !reference temperature K
END MODULE cable_photosynthetic_constants
