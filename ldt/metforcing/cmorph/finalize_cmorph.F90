!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_cmorph
! \label{finalize_cmorph}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine finalize_cmorph(findex)

! !USES:
  use cmorph_forcingMod, only : cmorph_struc
! !DESCRIPTION: 
!  Routine to cleanup CMORPH forcing related memory allocations.   
!
!EOP
  implicit none
  integer, intent(IN) :: findex  

  deallocate(cmorph_struc)

end subroutine finalize_cmorph
