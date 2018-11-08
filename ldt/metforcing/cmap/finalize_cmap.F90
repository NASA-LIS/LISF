!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_cmap
! \label{finalize_cmap}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine finalize_cmap(findex)

! !USES:
  use cmap_forcingMod, only : cmap_struc

! !DESCRIPTION: 
!  Routine to cleanup CMAP forcing related memory allocations.   
!
!EOP
  implicit none
  
  integer, intent(IN) :: findex

  deallocate(cmap_struc)

end subroutine finalize_cmap
