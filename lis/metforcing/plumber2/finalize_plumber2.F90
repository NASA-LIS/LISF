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
! !ROUTINE: finalize_plumber2
! \label{finalize_plumber2}
! 
! !REVISION HISTORY: 
! 16 Sep 2021: Mark Beauharnois, based on Bondville reader
! 
! !INTERFACE:
subroutine finalize_plumber2(findex)

  ! !USES:
  use plumber2_forcingMod, only : plumber2_struc
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  ! !DESCRIPTION:
  !  Routine to cleanup PLUMBER2 forcing related memory allocations.
  !
  !EOP
  implicit none

  integer, intent(IN) :: findex

  !      write(LIS_logunit,*) 'starting finalize_plumber2'
  deallocate(plumber2_struc)

end subroutine finalize_plumber2
