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
! !ROUTINE: reset_nldas20
!  \label{reset_nldas20}
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from reset_nldas2.F90)
!
! !INTERFACE:
subroutine reset_nldas20()
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use nldas20_forcingMod, only : nldas20_struc

  implicit none
!
! !DESCRIPTION:
!  Routine to reset NLDAS-2 forcing related memory allocations.
!
!EOP
  integer   :: n
  integer   :: findex

  do n = 1,LIS_rc%nnest
     nldas20_struc(n)%nldas20time1 = 3000.0
     nldas20_struc(n)%nldas20time2 = 0.0
  enddo

end subroutine reset_nldas20

