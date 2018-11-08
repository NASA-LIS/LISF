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

MODULE cable_math_constants
  USE cable_dimensions, ONLY : i_d, r_1
  IMPLICIT NONE
  PRIVATE i_d, r_1
  REAL(r_1), PARAMETER :: pi = 3.141592653589793238462643383279502884197
  REAL(r_1), PARAMETER :: pi180 = pi / 180.0 ! radians / degree
  REAL(r_1), PARAMETER :: two_pi = 2.0 * pi
END MODULE cable_math_constants
