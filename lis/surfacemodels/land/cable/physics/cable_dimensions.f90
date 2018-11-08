! cable_dimensions.f90
!
! Source file containing dimensions variables for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach

MODULE cable_dimensions
  INTEGER, PARAMETER    :: mp = 1  ! # land grid cells
  INTEGER, PARAMETER    :: mp_patch = 1 ! # patches(>=mp)one land grid pt can contain many patches
  INTEGER, PARAMETER    :: max_vegpatches = 1 ! The maximum # of patches in any grid cell
  INTEGER, PARAMETER    :: mf = 2  ! # leaves (sunlit, shaded)
  INTEGER, PARAMETER    :: nrb = 3 ! # radiation bands
  INTEGER, PARAMETER    :: ms = 6  ! # soil layers
  INTEGER, PARAMETER    :: msn = 3 ! max # snow layers
  INTEGER, PARAMETER    :: ncp = 3 ! # vegetation carbon stores
  INTEGER, PARAMETER    :: ncs = 2 ! # soil carbon stores
  ! i_d is default kind for representing integer values.
  INTEGER, PARAMETER :: i_d = KIND(9)
  ! r_1 is default kind for representing REAL values (typically 32 or 64 bits).
  INTEGER, PARAMETER :: r_1  = KIND(1.0)
  !INTEGER, PARAMETER :: r_1  = SELECTED_REAL_KIND(12, 50)
  ! r_2 is kind for representing REAL values with at least 10-digit precision
  ! (typically 64 bits).
  INTEGER, PARAMETER :: r_2  = SELECTED_REAL_KIND(12, 50)
  INTEGER, PARAMETER :: ncstringlength = 25 ! max length of string variable in netcdf file
  ! define single precision for netcdf in- and output
  INTEGER, PARAMETER :: SP  = KIND(1.0)
  INTEGER, PARAMETER :: LGT = KIND(.true.)

!  !$OMP THREADPRIVATE(mp,mp_patch)

END MODULE cable_dimensions
