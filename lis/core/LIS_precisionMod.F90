!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_precisionMod
!BOP
!
! !MODULE: LIS_precisionMod
! 
! !REVISION HISTORY:
! 
! 14 Nov 2002; Sujay Kumar  Initial Specification 
!
! !ARGUMENTS:
  integer, parameter :: r4 = selected_real_kind(5)
  integer, parameter :: r8 = selected_real_kind(6)
  integer, parameter :: i8 = selected_int_kind(13)

#if (defined DOUBLE_PRECISION)
  real(r8), parameter :: inf = O'777600000000000000000'
  real(r8), parameter :: nan = O'777677777777777777777'
#else

#if (defined SYSLINUX)
  real(r8), parameter :: inf = Z'7F800000'
  real(r8), parameter :: nan = Z'7FC00000'
#else
!  real(r8), parameter :: inf = O'17740000000'
!  real(r8), parameter :: nan = O'17757777777'
  real(r8), parameter :: inf = 0
  real(r8), parameter :: nan = 0
#endif


#endif
!  integer, parameter  :: bigint = 100000000
  integer, parameter  :: bigint = 100
!
! !DESCRIPTION: 
!  Define the precision to use for floating point and integer operations
!  throughout the model.
!
!EOP

end module LIS_precisionMod
 
