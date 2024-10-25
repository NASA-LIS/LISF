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
!BOP
! !MODULE: reset_chirps2
! \label{reset_chirps2}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 09Oct2015; Kristi Arsenault, Added to CHIRPS2 reader
! 
! !INTERFACE:
subroutine reset_chirps2
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use chirps2_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for CHIRPS v2 forcing. 
!
!EOP  
  implicit none

  integer :: n

  do n=1,LIS_rc%nnest
     chirps2_struc(:)%chirpstime1 = 0.
     chirps2_struc(:)%chirpstime2 = 0.
     chirps2_struc(:)%reset_flag = .true.
  enddo

end subroutine reset_chirps2
