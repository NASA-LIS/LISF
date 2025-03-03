!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine setup_ALMIPIIroughness(n)

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_vegDataMod,   only : LIS_roughness

  integer, intent(in)  :: n 
  integer              :: rc

  LIS_roughness(n)%roughnessInterval  = 864000
  LIS_roughness(n)%roughnessIntervalType = "10-day" !10-day

  call LIS_registerAlarm("LIS roughness read alarm",LIS_rc%ts, &
       LIS_roughness(n)%roughnessInterval,&
       intervalType=LIS_roughness(n)%roughnessIntervalType) 

  call ESMF_ConfigFindLabel(LIS_config,&
       "ALMIPII roughness data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config, &
       LIS_roughness(n)%roughnessfile,rc=rc)
  call LIS_verify(rc,'ALMIPII roughness data directory: not defined')
end subroutine setup_ALMIPIIroughness



