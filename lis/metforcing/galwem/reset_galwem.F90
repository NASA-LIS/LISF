!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: reset_galwem
!  \label{reset_galwem}
!
! !REVISION HISTORY:
! 05 Apr 2022; Yeosang Yoon, Initial Code
!
! !INTERFACE:
subroutine reset_galwem()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use galwem_forcingMod, only : galwem_struc
!
! !DESCRIPTION:
!  Routine to reset GALWEM forcing related memory allocations.
!
!EOP
  implicit none

  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     galwem_struc(n)%fcsttime1 = 3000.0
     galwem_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine reset_galwem
