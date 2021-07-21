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
! !ROUTINE: RAPID_routing_readrst
! \label{RAPID_routing_readrst}
!
! !REVISION HISTORY:
! 18 Mar 2021: Yeosang Yoon;  Initial implementation

subroutine RAPID_routing_output(n)
  
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_histDataMod
  use LIS_historyMod
  use LIS_fileIOMod
  use RAPID_routingMod

  implicit none
  
  integer, intent(in)   :: n 
  
  character(len=12)     :: cdate1
  integer               :: iret
  character*100         :: filename
  character*100         :: name
  integer               :: ftn
  integer               :: mo, da
  logical               :: open_stats
  logical               :: alarmCheck

! use rapid_write_Qout_file.F90


end subroutine RAPID_routing_output
