!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: reset_nldas2
!  \label{reset_nldas2}
!
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_nldas2()
! !USES:
  use LDT_coreMod,       only : LDT_rc
  use nldas2_forcingMod, only : nldas2_struc
!
! !DESCRIPTION:
!  Routine to reset nldas2 forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LDT_rc%nnest
     nldas2_struc(n)%nldas2time1 = 3000.0
     nldas2_struc(n)%nldas2time2 = 0.0
  enddo

end subroutine reset_nldas2
