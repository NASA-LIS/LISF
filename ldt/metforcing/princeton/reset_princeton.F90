!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit System (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: reset_princeton
!  \label{reset_princeton}
!
! !REVISION HISTORY: 
! 08/21/2014 Bailing Li;initial creation 
! 
! !INTERFACE:
subroutine reset_princeton()
! !USES:
  use LDT_coreMod,          only : LDT_rc
  use princeton_forcingMod, only : princeton_struc
!
! !DESCRIPTION:
!  Routine to reset nldas2 forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LDT_rc%nnest
     princeton_struc(n)%princetontime1 = 3000.0
     princeton_struc(n)%princetontime2 = 0.0
  enddo

end subroutine reset_princeton
