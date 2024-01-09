!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: reset_mogrepsg
!  \label{reset_mogrepsg}
!
! !REVISION HISTORY:
! 26 Jan 2023; Yeosang Yoon, Initial Code
!
! !INTERFACE:
subroutine reset_mogrepsg()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use mogrepsg_forcingMod, only : mogrepsg_struc
!
! !DESCRIPTION:
!  Routine to reset GALWEM-GE forcing related memory allocations.
!
!EOP
  implicit none

  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     mogrepsg_struc(n)%fcsttime1 = 3000.0
     mogrepsg_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine reset_mogrepsg
