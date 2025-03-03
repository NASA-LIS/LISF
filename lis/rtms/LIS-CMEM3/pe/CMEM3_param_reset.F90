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
! !ROUTINE: CMEM3_param_reset
!  \label{CMEM3_param_reset}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine CMEM3_param_reset()
!!$! !USES:
!!$  use ESMF
!!$  use LIS_coreMod, only : LIS_rc
!!$  use LIS_timeMgrMod, only : LIS_calendar, LIS_clock
!!$  use LIS_logMod,  only : LIS_verify
!!$  use CMEM3_Mod !,      only : CMEM3_data, CMEM3_ctl
!!$
!!$  implicit none
!!$! !ARGUMENTS: 
!!$!
!!$! !DESCRIPTION:
!!$!  
!!$
!!$! 
!!$!EOP
!!$  type(ESMF_TimeInterval) :: alarmInterval
!!$  type(ESMF_Time)         :: time2, alarmTime, currTime
!!$  type(ESMF_TimeInterval) :: deltat
!!$  integer                 :: m, n 
!!$  integer                 :: rc
!!$
!!$
!!$  n = 1
!!$
!!$  CMEM3_ctl%ls_status = 0 
!!$  CMEM3_ctl%ls_accum_status = 0 
!!$
!!$!need to reinitialize all alarms: 
!!$  call ESMF_TimeIntervalSet(alarmInterval, h=CMEM3_ctl%dt, rc=rc)
!!$  call LIS_verify(rc, 'Time interval set failed in CMEM3_init')
!!$  
!!$  call ESMF_TimeSet(time2, yy = LIS_rc%yr, &
!!$       mm = LIS_rc%mo, &
!!$       dd = LIS_rc%da, &
!!$       h  = 0, &
!!$       m  = 0, & 
!!$       s  = 0, &
!!$       calendar = LIS_calendar,&
!!$       rc = rc)
!!$  call LIS_verify(rc, 'Error in time2 set in CMEM3_init')
!!$  
!!$  alarmTime = time2 + alarmInterval
!!$  
!!$  CMEM3_ctl%dt_alarm = ESMF_AlarmCreate(clock=LIS_clock, &
!!$       ringTime = alarmTime, &
!!$       ringInterval = alarmInterval, rc=rc)
!!$  call LIS_verify(rc, 'Alarm Create failed in CMEM3_init')
!!$  
!!$  !accumulation alarm 1
!!$  call ESMF_TimeIntervalSet(alarmInterval, h=CMEM3_ctl%acc_interval1, rc=rc)
!!$  call LIS_verify(rc, 'Time interval set failed in CMEM3_init')
!!$  
!!$  call ESMF_TimeSet(time2, yy = LIS_rc%yr, &
!!$       mm = LIS_rc%mo, &
!!$       dd = LIS_rc%da, &
!!$       h  = 0, &
!!$       m  = 0, & 
!!$       s  = 0, &
!!$       calendar = LIS_calendar,&
!!$       rc = rc)
!!$  call LIS_verify(rc, 'Error in time2 set in CMEM3_init')
!!$  
!!$  alarmTime = time2 + alarmInterval
!!$  
!!$  CMEM3_ctl%acc1_alarm = ESMF_AlarmCreate(clock=LIS_clock, &
!!$       ringTime = alarmTime, &
!!$       ringInterval = alarmInterval, rc=rc)
!!$  call LIS_verify(rc, 'Alarm Create failed in CMEM3_init')
!!$  
!!$  !accumulation alarm 2
!!$  call ESMF_TimeIntervalSet(alarmInterval, h=CMEM3_ctl%acc_interval2, rc=rc)
!!$  call LIS_verify(rc, 'Time interval set failed in CMEM3_init')
!!$  
!!$  call ESMF_TimeSet(time2, yy = LIS_rc%yr, &
!!$       mm = LIS_rc%mo, &
!!$       dd = LIS_rc%da, &
!!$       h  = 0, &
!!$       m  = 0, & 
!!$       s  = 0, &
!!$       calendar = LIS_calendar,&
!!$       rc = rc)
!!$  call LIS_verify(rc, 'Error in time2 set in CMEM3_init')
!!$       
!!$  alarmTime = time2 + alarmInterval
!!$  
!!$  call ESMF_TimeIntervalSet(alarmInterval, h=CMEM3_ctl%acc_interval1, rc=rc)
!!$  call LIS_verify(rc, 'Time interval set failed in CMEM3_init')
!!$  
!!$  CMEM3_ctl%acc2_alarm = ESMF_AlarmCreate(clock=LIS_clock, &
!!$       ringTime = alarmTime, &
!!$       ringInterval = alarmInterval, rc=rc)
!!$  call LIS_verify(rc, 'Alarm Create failed in CMEM3_init')
!!$       
!!$  !accumulation alarm 3
!!$  call ESMF_TimeIntervalSet(alarmInterval, h=CMEM3_ctl%acc_interval3, rc=rc)
!!$  call LIS_verify(rc, 'Time interval set failed in CMEM3_init')
!!$  
!!$  call ESMF_TimeSet(time2, yy = LIS_rc%yr, &
!!$       mm = LIS_rc%mo, &
!!$       dd = LIS_rc%da, &
!!$       h  = 0, &
!!$       m  = 0, & 
!!$       s  = 0, &
!!$       calendar = LIS_calendar,&
!!$       rc = rc)
!!$  call LIS_verify(rc, 'Error in time2 set in CMEM3_init')
!!$  
!!$  alarmTime = time2 + alarmInterval
!!$       
!!$  call ESMF_TimeIntervalSet(alarmInterval, h=CMEM3_ctl%acc_interval1, rc=rc)
!!$  call LIS_verify(rc, 'Time interval set failed in CMEM3_init')
!!$  
!!$       
!!$  CMEM3_ctl%acc3_alarm = ESMF_AlarmCreate(clock=LIS_clock, &
!!$       ringTime = alarmTime, &
!!$       ringInterval = alarmInterval, rc=rc)
!!$  call LIS_verify(rc, 'Alarm Create failed in CMEM3_init')
!!$
end subroutine CMEM3_param_reset

