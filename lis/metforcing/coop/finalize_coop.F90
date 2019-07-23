!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_coop
! \label{finalize_coop}
! 
! !REVISION HISTORY: 
! 08Jun2011: Yuqiong Liu; Initial Specification
! 
! !INTERFACE:
subroutine finalize_coop(findex)

! !USES:
  use coop_forcingMod, only : coop_struc
! !DESCRIPTION: 
!  Routine to cleanup COOP forcing related memory allocations.   
!
!EOP
  implicit none
  integer, intent(IN) :: findex

  deallocate(coop_struc)

end subroutine finalize_coop
