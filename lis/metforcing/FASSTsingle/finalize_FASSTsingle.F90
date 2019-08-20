!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_FASSTsingle
! \label{finalize_FASSTsingle}
! 
! !REVISION HISTORY: 
! 08 Dec 2004: Sujay Kumar; Initial Specification
! 05 Oct 2010: David Mocko, Updated for FASST single-point test case
!
! !INTERFACE:
subroutine finalize_FASSTsingle(findex)

  ! !USES:
  use FASSTsingle_forcingMod, only : FASSTsingle_struc
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  ! !DESCRIPTION:
  !  Routine to cleanup FASSTsingle forcing related memory allocations.
  !
  !EOP
  implicit none

  integer, intent(IN) :: findex

  !      write(LIS_logunit,*) 'starting finalize_FASSTsingle'
  deallocate(FASSTsingle_struc)

end subroutine finalize_FASSTsingle

