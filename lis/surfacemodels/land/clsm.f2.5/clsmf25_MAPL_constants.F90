!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module clsmf25_MAPL_constants

!  $Id: MAPL_Constants.F90,v 1.2 2009-01-14 14:10:21 f4mjs Exp $

implicit none
private


!=============================================================================
!BOP

! !MODULE: -- A container module for global constants

! !PUBLIC VARIABLES:

!EOP

!BOC
real(kind=8), parameter, public :: MAPL_PI_R8     = 3.14159265358979323846
real, parameter, public :: MAPL_PI     = MAPL_PI_R8
real, parameter, public :: MAPL_GRAV   = 9.80                   ! m^2/s
real, parameter, public :: MAPL_RADIUS = 6376.0E3               ! m
real, parameter, public :: MAPL_OMEGA  = 2.0*MAPL_PI/86164.0    ! 1/s
real, parameter, public :: MAPL_ALHL   = 2.4665E6               ! J/kg @15C
real, parameter, public :: MAPL_ALHF   = 3.3370E5               ! J/kg
real, parameter, public :: MAPL_ALHS   = MAPL_ALHL+MAPL_ALHF    ! J/kg
real, parameter, public :: MAPL_STFBOL = 5.6734E-8              ! W/(m^2 K^4)
real, parameter, public :: MAPL_AIRMW  = 28.97                  ! kg/Kmole
real, parameter, public :: MAPL_H2OMW  = 18.01                  ! kg/Kmole
real, parameter, public :: MAPL_O3MW   = 47.9982                ! kg/Kmole
real, parameter, public :: MAPL_RUNIV  = 8314.3                 ! J/(Kmole K)
real, parameter, public :: MAPL_KAPPA  = 2.0/7.0                ! --
real, parameter, public :: MAPL_RVAP   = MAPL_RUNIV/MAPL_H2OMW  ! J/(kg K)
real, parameter, public :: MAPL_RGAS   = MAPL_RUNIV/MAPL_AIRMW  ! J/(kg K)
real, parameter, public :: MAPL_CP     = MAPL_RGAS/MAPL_KAPPA   ! J/(kg K)
real, parameter, public :: MAPL_P00    = 100000.0               ! Pa
real, parameter, public :: MAPL_CAPICE = 2000.                  ! J/(K kg)
real, parameter, public :: MAPL_CAPWTR = 4218.                  ! J/(K kg)
real, parameter, public :: MAPL_RHOWTR = 1000.                  ! kg/m^3
real, parameter, public :: MAPL_NUAIR  = 1.533E-5               ! m^2/S (@ 18C)
real, parameter, public :: MAPL_TICE   = 273.16                 ! K
real, parameter, public :: MAPL_SRFPRS = 98470                  ! Pa
real, parameter, public :: MAPL_KARMAN = 0.40                   ! --
real, parameter, public :: MAPL_USMIN  = 1.00                   ! m/s
real, parameter, public :: MAPL_VIREPS = MAPL_AIRMW/MAPL_H2OMW-1.0   ! --
real, parameter, public :: MAPL_AVOGAD = 6.023E26               ! 1/kmol

integer,parameter, public :: MAPL_R8 = selected_real_kind(12) ! 8 byte real
integer,parameter, public :: MAPL_R4 = selected_real_kind( 6) ! 4 byte real
integer,parameter, public :: MAPL_RN = kind(1.0)              ! native real
integer,parameter, public :: MAPL_I8 = selected_int_kind (13) ! 8 byte integer
integer,parameter, public :: MAPL_I4 = selected_int_kind ( 6) ! 4 byte integer
integer,parameter, public :: MAPL_IN = kind(1)                ! native integer
!EOC

end module clsmf25_MAPL_constants

