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
