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
! !MODULE: reset_gddp
!  \label{reset_gddp}
!
! !REVISION HISTORY: 
! 03 Feb 2022; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_gddp()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use gddp_forcingMod, only : gddp_struc
!
! !DESCRIPTION:
!  Routine to reset GDDP forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     gddp_struc(n)%gddptime1 = 3000.0
     gddp_struc(n)%gddptime2 = 0.0
  enddo

end subroutine reset_gddp
