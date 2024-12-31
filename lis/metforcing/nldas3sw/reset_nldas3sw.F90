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
! !ROUTINE: reset_nldas3sw
!  \label{reset_nldas3sw}
!
! !REVISION HISTORY:
! 27 Dec 2024: David Mocko, Initial Specification
!                           (derived from reset_nldas20.F90)
!
! !INTERFACE:
subroutine reset_nldas3sw()
! !USES:
  use LIS_coreMod,         only : LIS_rc
  use nldas3sw_forcingMod, only : nldas3sw_struc

  implicit none
!
! !DESCRIPTION:
!  Routine to reset CERES SWdown forcing related memory allocations.
!
!EOP
  integer   :: n
  integer   :: findex

  do n = 1,LIS_rc%nnest
     nldas3sw_struc(n)%nldas3swtime1 = 3000.0
     nldas3sw_struc(n)%nldas3swtime2 = 0.0
  enddo

end subroutine reset_nldas3sw

