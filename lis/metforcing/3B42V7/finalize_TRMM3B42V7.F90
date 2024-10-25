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
! !ROUTINE: finalize_TRMM3B42V7
! \label{finalize_TRMM3B42V7}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine finalize_TRMM3B42V7(findex)

! !USES:
  use TRMM3B42V7_forcingMod, only : TRMM3B42V7_struc
! !DESCRIPTION: 
!  Routine to cleanup TRMM3B42V7 forcing related memory allocations.   
!
!EOP
  implicit none
  integer :: findex
  
  deallocate(TRMM3B42V7_struc)

end subroutine finalize_TRMM3B42V7
