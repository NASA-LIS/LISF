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
! !ROUTINE: finalize_snotel
! \label{finalize_snotel}
! 
! !REVISION HISTORY: 
! 08Jun2011: Yuqiong Liu; Initial Specification
! 
! !INTERFACE:
subroutine finalize_snotel(findex)

! !USES:
  use snotel_forcingMod, only : snotel_struc
! !DESCRIPTION: 
!  Routine to cleanup SCAN forcing related memory allocations.   
!
!EOP
  implicit none
  integer, intent(IN) :: findex

  deallocate(snotel_struc)

end subroutine finalize_snotel
