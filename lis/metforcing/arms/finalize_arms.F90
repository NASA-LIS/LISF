!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_arms
! \label{finalize_arms}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine finalize_arms(findex)

! !USES:
  use arms_forcingMod, only : arms_struc
! !DESCRIPTION: 
!  Routine to cleanup ARMS forcing related memory allocations.   
!
!EOP
  implicit none
  integer, intent(in) :: findex

  deallocate(arms_struc)

end subroutine finalize_arms
