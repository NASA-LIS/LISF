!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: reset_cmap
! \label{reset_cmap}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_cmap
! !USES:
  use LDT_coreMod,    only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use cmap_forcingMod
!
! !DESCRIPTION:
!  Routine to reset gridchange fields for CMAP forcing. 
!
!EOP  
  implicit none

  integer :: n 
  
  do n=1,LDT_rc%nnest

     cmap_struc(n)%gridchange1 = .true.
     cmap_struc(n)%gridchange2 = .true.
     cmap_struc(n)%gridchange3 = .true.
     cmap_struc(n)%gridchange4 = .true.
     cmap_struc(n)%gridchange5 = .true.

  enddo

end subroutine reset_cmap
