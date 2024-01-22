!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: reset_merra2_ac
! \label{reset_merra2_ac}
! 
! !REVISION HISTORY: 
! 01 Jun 2022: Michel Bechtold, initial code (based on merra-2 data preprocessed
! to daily data)
! 17 Jan 2024: Louise Busschaert, AC71 implementation in NASA master
! 
! !INTERFACE:
subroutine reset_merra2_ac
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use merra2_ac_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for merra2_ac forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LIS_rc%nnest
     merra2_ac_struc(n)%startFlag = .true. 
     merra2_ac_struc(n)%dayFlag = .true. 
     merra2_ac_struc(n)%merra2time1 = 3000.0
     merra2_ac_struc(n)%merra2time2 = 0.0
     merra2_ac_struc(n)%ringtime = 0.0
     merra2_ac_struc(n)%reset_flag = .true.
  enddo
end subroutine reset_merra2_ac
