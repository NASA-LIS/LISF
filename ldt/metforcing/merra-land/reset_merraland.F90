!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: reset_merraland
! \label{reset_merraland}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 22Jan2015: KR Arsenault, added to LDT
! 
! !INTERFACE:
subroutine reset_merraland
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use merraland_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for MERRALAND forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LDT_rc%nnest
     merraland_struc(n)%startFlag = .true. 
     merraland_struc(n)%dayFlag = .true. 
     merraland_struc(n)%merralandtime1 = 3000.0
     merraland_struc(n)%merralandtime2 = 0.0
  enddo
end subroutine reset_merraland
