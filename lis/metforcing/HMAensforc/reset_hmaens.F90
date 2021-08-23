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
! !MODULE: reset_hmaens
! \label{reset_hmaens}
! 
! !REVISION HISTORY: 
! 23 dec 2019: Sujay Kumar, initial code 
! 
! !INTERFACE:
subroutine reset_hmaens
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use hmaens_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for hmaens forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LIS_rc%nnest
     hmaens_struc(n)%startFlag = .true. 
     hmaens_struc(n)%dayFlag = .true. 
     hmaens_struc(n)%hmaenstime1 = 3000.0
     hmaens_struc(n)%hmaenstime2 = 0.0
     hmaens_struc(n)%ringtime = 0.0
     hmaens_struc(n)%reset_flag = .true.
  enddo
end subroutine reset_hmaens

