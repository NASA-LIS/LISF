!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: reset_merraland
! \label{reset_merraland}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_merraland
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use merraland_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for MERRALAND forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LIS_rc%nnest
     merraland_struc(n)%startFlag = .true. 
     merraland_struc(n)%dayFlag = .true. 
     merraland_struc(n)%merralandtime1 = 3000.0
     merraland_struc(n)%merralandtime2 = 0.0
  enddo
end subroutine reset_merraland
