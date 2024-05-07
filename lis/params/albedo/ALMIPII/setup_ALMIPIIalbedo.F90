!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine setup_ALMIPIIalbedo(n)
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_albedoMod,   only : LIS_alb

  integer, intent(in)  :: n 
  
  integer              :: rc

  LIS_alb(n)%albInterval  = 864000
  LIS_alb(n)%albIntervalType = "10-day" !10-day

  call LIS_registerAlarm("LIS alb read alarm",LIS_rc%ts, &
       LIS_alb(n)%albInterval,&
       intervalType=LIS_alb(n)%albIntervalType)

  call ESMF_ConfigFindLabel(LIS_config,&
       "ALMIPII albedo data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config, &
       LIS_alb(n)%albfile,rc=rc)
  call LIS_verify(rc,'ALMIPII albedo data directory: not defined')

end subroutine setup_ALMIPIIalbedo



