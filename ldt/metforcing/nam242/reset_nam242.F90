!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: reset_nam242
! \label{reset_nam242}
! 
! !REVISION HISTORY: 
!     Sep 2012: NOHRSC/NOAA: Initial specification
! 
! !INTERFACE:
subroutine reset_nam242
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use nam242_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for NAM forcing. 
!
!EOP  
  implicit none

  real :: gridDesci(20)
  integer :: n 
  
  do n=1,LDT_rc%nnest
     gridDesci = 0
     gridDesci(1) = 5
     gridDesci(2) = 553
     gridDesci(3) = 425
     gridDesci(4) = 30
     gridDesci(5) = -173
     gridDesci(6) = 0       ! not used
     gridDesci(7) = -135    ! Dagang question
     gridDesci(8) = 11.25
     gridDesci(9) = 11.25
     gridDesci(10) = 60     ! Dagang question
     gridDesci(11) = -135   ! Dagang question
     gridDesci(20) = 0      ! Dagang question
     nam242_struc(n)%nc = 553
     nam242_struc(n)%nr = 425
     nam242_struc(n)%mi    = nam242_struc(n)%nc*nam242_struc(n)%nr
  enddo
end subroutine reset_nam242
