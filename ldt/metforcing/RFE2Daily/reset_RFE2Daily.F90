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
! !MODULE: reset_RFE2Daily
! \label{reset_RFE2Daily}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_RFE2Daily
! !USES:
  use LDT_coreMod,    only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use RFE2Daily_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for RFE2DAILY forcing. 
!
!EOP  
  implicit none

  real :: gridDesci(20)
  integer :: n 

  do n=1,LDT_rc%nnest
     RFE2Daily_struc(n)%RFE2DailyEndTime = LDT_rc%udef  
  enddo

end subroutine reset_RFE2Daily
