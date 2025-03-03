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
  use LIS_coreMod,     only : LIS_rc
  use cmap_forcingMod, only : cmap_struc
! !DESCRIPTION: 
!  Routine to cleanup CMAP forcing related memory allocations.   
!
!EOP
  implicit none
  
  integer, intent(in) :: findex

  integer :: n
  integer :: rc

  do n=1,LIS_rc%nnest
     deallocate(cmap_struc(n)%n111, stat=rc)
     deallocate(cmap_struc(n)%n121, stat=rc)
     deallocate(cmap_struc(n)%n211, stat=rc)
     deallocate(cmap_struc(n)%n221, stat=rc)
     deallocate(cmap_struc(n)%w111, stat=rc)
     deallocate(cmap_struc(n)%w121, stat=rc)
     deallocate(cmap_struc(n)%w211, stat=rc)
     deallocate(cmap_struc(n)%w221, stat=rc)

     deallocate(cmap_struc(n)%n112, stat=rc)
     deallocate(cmap_struc(n)%n122, stat=rc)
     deallocate(cmap_struc(n)%n212, stat=rc)
     deallocate(cmap_struc(n)%n222, stat=rc)
     deallocate(cmap_struc(n)%w112, stat=rc)
     deallocate(cmap_struc(n)%w122, stat=rc)
     deallocate(cmap_struc(n)%w212, stat=rc)
     deallocate(cmap_struc(n)%w222, stat=rc)
  enddo

  deallocate(cmap_struc)

end subroutine finalize_cmap
