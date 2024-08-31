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
! !ROUTINE: finalize_Bondville
! \label{finalize_Bondville}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine finalize_Bondville(findex)

  ! !USES:
  use Bondville_forcingMod, only : Bondville_struc
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  ! !DESCRIPTION:
  !  Routine to cleanup Bondville forcing related memory allocations.
  !
  !EOP
  implicit none

  integer, intent(IN) :: findex

  !      write(LIS_logunit,*) 'starting finalize_Bondville'
  deallocate(Bondville_struc)

end subroutine finalize_Bondville

