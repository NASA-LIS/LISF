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
! !ROUTINE: finalize_Loobos
! \label{finalize_Loobos}
! 
! !REVISION HISTORY: 
! 08Dec2004: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine finalize_Loobos(findex)

  ! !USES:
  use Loobos_forcingMod, only : Loobos_struc
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  ! !DESCRIPTION:
  !  Routine to cleanup Loobos forcing related memory allocations.
  !
  !EOP
  implicit none

  integer, intent(IN) :: findex

  !      write(LIS_logunit,*) 'starting finalize_Loobos'
  deallocate(Loobos_struc)

end subroutine finalize_Loobos

