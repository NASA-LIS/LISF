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
! !MODULE: reset_merra2
! \label{reset_merra2}
! 
! !REVISION HISTORY: 
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
! 
! !INTERFACE:
subroutine reset_merra2
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use merra2_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for merra2 forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LDT_rc%nnest
     merra2_struc(n)%startFlag = .true. 
     merra2_struc(n)%dayFlag = .true. 
     merra2_struc(n)%merra2time1 = 3000.0  ! original value
!     merra2_struc(n)%merra2time1 = 0.0  ! alternate value (KRA)
     merra2_struc(n)%merra2time2 = 0.0
     merra2_struc(n)%ringtime = 0.0
     merra2_struc(n)%reset_flag = .true.
  enddo
end subroutine reset_merra2
