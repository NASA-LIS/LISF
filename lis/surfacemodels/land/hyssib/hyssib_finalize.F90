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
! !ROUTINE: hyssib_finalize
! \label{hyssib_finalize}
!
! !REVISION HISTORY:
!
! 28 Apr 2002: Sujay Kumar; Initial code
! 05 Sep 2007: Sujay Kumar; Addition of Hyssib to LIS 5.0
! 
! !INTERFACE:
subroutine hyssib_finalize()
! !USES:
  use LIS_coreMod, only : LIS_rc
  use hyssib_lsmMod
  use hyssibveg_module
  use hyssibalb_module
!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in Hyssib
!  
!EOP
  implicit none
  integer :: n
  
  do n=1,LIS_rc%nnest
     deallocate(hyssib_struc(n)%hyssib)
  enddo
  deallocate(hyssib_struc)
  !deallocate veg and albedo structures
  deallocate(hyssibveg)
  deallocate(hyssibalb)

end subroutine hyssib_finalize
