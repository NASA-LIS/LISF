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
! !MODULE: reset_era5
! \label{reset_era5}
! 
! !REVISION HISTORY: 
! 23 dec 2019: Sujay Kumar, initial code 
! 
! !INTERFACE:
subroutine reset_era5
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use era5_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for era5 forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LDT_rc%nnest
     era5_struc(n)%startFlag = .true. 
     era5_struc(n)%dayFlag = .true. 
     era5_struc(n)%era5time1 = 3000.0
     era5_struc(n)%era5time2 = 0.0
     era5_struc(n)%ringtime = 0.0
     era5_struc(n)%reset_flag = .true.
  enddo
end subroutine reset_era5
