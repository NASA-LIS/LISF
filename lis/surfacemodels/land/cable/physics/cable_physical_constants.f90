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
MODULE cable_physical_constants
  USE cable_dimensions, ONLY : i_d, r_1
  IMPLICIT NONE
  PRIVATE i_d, r_1
  REAL(r_1), PARAMETER :: capp   = 1004.64  ! air spec. heat capacity (J/kg/K)
  REAL(r_1), PARAMETER :: dheat  = 21.5E-6  ! molecular diffusivity for heat
  REAL(r_1), PARAMETER :: grav   = 9.80     ! gravity acceleration (m/s2)
  REAL(r_1), PARAMETER :: rgas   = 8.3143   ! universal gas const  (J/mol/K)
  REAL(r_1), PARAMETER :: rmair  = 0.02897  ! molecular wt: dry air (kg/mol)
  REAL(r_1), PARAMETER :: rmh2o  = 0.018016 ! molecular wt: water	(kg/mol)
  REAL(r_1), PARAMETER :: sboltz = 5.67e-8  ! Stefan-Boltz. constant (W/m2/K4)
  REAL(r_1), PARAMETER :: tfrz   = 273.16   ! Temp (K) corresp. to 0 C
  ! Teten coefficients
  REAL(r_1), PARAMETER :: tetena = 6.106  ! ??? refs?
  REAL(r_1), PARAMETER :: tetenb = 17.27
  REAL(r_1), PARAMETER :: tetenc = 237.3
  ! Aerodynamic parameters, diffusivities, water density:
  REAL(r_1), PARAMETER :: vonk   = 0.40 ! von Karman constant
  REAL(r_1), PARAMETER :: a33    = 1.25 ! inertial sublayer sw/us
  REAL(r_1), PARAMETER :: csw    = 0.50 ! canopy sw decay (Weil theory)
  REAL(r_1), PARAMETER :: ctl    = 0.40 ! Wagga wheat (RDD 1992, Challenges)
  REAL(r_1), PARAMETER :: apol   = 0.70 ! Polhausen coeff: single-sided plate
  REAL(r_1), PARAMETER :: prandt = 0.71 ! Prandtl number: visc/diffh
  REAL(r_1), PARAMETER :: schmid = 0.60 ! Schmidt number: visc/diffw
  REAL(r_1), PARAMETER :: diffwc = 1.60 ! diffw/diffc = H2O/CO2 diffusivity
  REAL(r_1), PARAMETER :: rhow   = 1000.0 ! liquid water density   [kg/m3]
  REAL(r_1), PARAMETER :: emleaf = 1.0  ! leaf emissivity
  REAL(r_1), PARAMETER :: emsoil = 1.0  ! soil emissivity
  REAL(r_1), PARAMETER :: cr = 0.3	! element drag coefficient
  REAL(r_1), PARAMETER :: cs = 0.003    ! substrate drag coefficient
  REAL(r_1), PARAMETER :: beta  = cr/cs ! ratio cr/cs
  REAL(r_1), PARAMETER :: ccd   = 15.0  ! constant in d/h equation
  REAL(r_1), PARAMETER :: ccw   = 2.0   ! ccw=(zw-d)/(h-d)
  REAL(r_1), PARAMETER :: usuhm = 0.3   ! (max of us/uh)
  ! Turbulence parameters:
  INTEGER(i_d), PARAMETER :: niter = 4  ! number of iterations for za/L
  REAL(r_1), PARAMETER :: zetmul = 0.4  ! if niter=2, final zeta=zetmul*zetar(2)
  REAL(r_1), PARAMETER :: zeta0  = 0.0  ! initial value of za/L
  REAL(r_1), PARAMETER :: zetneg = -10.0 ! negative limit on za/L when niter>=3
  REAL(r_1), PARAMETER :: zetpos = 0.5  ! positive limit on za/L when niter>=3
  REAL(r_1), PARAMETER :: zdlin  = 1.0  ! height frac of d below which TL linear
  REAL(r_1), PARAMETER :: umin   = 1.0
END MODULE cable_physical_constants
