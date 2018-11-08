!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_TRMM3B42V6
! \label{finalize_TRMM3B42V6}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine finalize_TRMM3B42V6(findex)

! !USES:
  use TRMM3B42V6_forcingMod, only : TRMM3B42V6_struc
! !DESCRIPTION: 
!  Routine to cleanup TRMM3B42V6 forcing related memory allocations.   
!
!EOP
  implicit none
  integer :: findex
  
  deallocate(TRMM3B42V6_struc)

end subroutine finalize_TRMM3B42V6
