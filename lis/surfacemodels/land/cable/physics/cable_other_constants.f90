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

!=========================================================================
!=========================================================================
MODULE cable_other_constants
  USE cable_dimensions, ONLY : i_d, r_1, nrb
  IMPLICIT NONE
  PRIVATE i_d, r_1, nrb
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: gauss_w=(/0.308,0.514,0.178/) ! Gaussian integ. weights
  ! values in refl and taul are slightly modified since Oct07 and Mar08 (YP)
  ! leaf reflectance
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: refl = (/ 0.1, 0.425, 0.02 /) ! mar08
  ! leaf transmittance
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: taul = (/ 0.1, 0.425, 0.02 /) ! mar08
  INTEGER(i_d), PARAMETER :: istemp = 4 ! soil temp:	 1,2,3,4 = FR,kf,mrr,mrrkf
  INTEGER(i_d), PARAMETER :: ismois = 2 ! soil moist:  1,2,3	 = MP84,NP89,Richards
  INTEGER(i_d), PARAMETER :: isinf  = 2 ! soil infilt: 1,2	 = MP84, FC96
  INTEGER(i_d), PARAMETER :: isevap = 2 ! soil evap: 1,2,3 = alfa,beta,threshold
  INTEGER(i_d), PARAMETER :: itherm = 1 ! VW or KGK algorithm for hconds,rkapps
  INTEGER(i_d), PARAMETER :: irktem = 5 ! RK steps in soil temp schemes
  INTEGER(i_d), PARAMETER :: irkmoi = 5 ! RK steps in soil moisture schemes
  ! soil water parameters:
  REAL(r_1), PARAMETER :: etarct = 0.7  ! rel soil moisture for finding zst1,zst2
  REAL(r_1), PARAMETER :: dbde = 1.3333 ! d(beta)/d(etar): 4/3, 8 for D78,KGK91
  INTEGER(i_d), PARAMETER :: istsw = 1  !
  INTEGER(i_d), PARAMETER :: iresp=0    ! unscaled (iresp=0) or scaled (iresp=1) respiration
END MODULE cable_other_constants
