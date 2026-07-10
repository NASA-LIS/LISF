!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: reset_chirps
! \label{reset_chirps}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 09Oct2015; Kristi Arsenault, Added to CHIRPS2 reader
! 29Oct2025; J Nattala, Updated to support both CHIRPS v2.0 and v3.0
! 
! !INTERFACE:
subroutine reset_chirps
! !USES:
  use LDT_coreMod,    only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use chirps_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for CHIRPS forcing. 
!
!EOP  
  implicit none

  integer :: n

  do n=1,LDT_rc%nnest
     chirps_struc(:)%chirpstime1 = 0.
     chirps_struc(:)%chirpstime2 = 0.
     chirps_struc(:)%reset_flag = .true.
  enddo

end subroutine reset_chirps
