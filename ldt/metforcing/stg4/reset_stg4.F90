!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: reset_stg4
!  \label{reset_stg4}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_stg4()

! !USES:
  use LDT_coreMod, only : LDT_rc
  use stg4_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup Stage IV forcing related memory allocations.   
! 
!EOP
  implicit none
  integer   :: n

  do n=1,LDT_rc%nnest  
     stg4_struc(n)%stg4time = 0.0
  enddo

end subroutine reset_stg4
