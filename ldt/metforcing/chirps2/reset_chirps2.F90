!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
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
  use LDT_coreMod,    only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use chirps2_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for CHIRPS v2 forcing. 
!
!EOP  
  implicit none

  integer :: n 

  do n=1,LDT_rc%nnest
     chirps2_struc(:)%chirpstime1 = 0.
     chirps2_struc(:)%chirpstime2 = 0.
     chirps2_struc(:)%reset_flag = .true.
  enddo

end subroutine reset_chirps2
