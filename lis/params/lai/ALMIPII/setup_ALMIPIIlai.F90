!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine setup_ALMIPIIlai(n)
  
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_vegDataMod,   only : LIS_lai

  integer, intent(in)  :: n 
  integer              :: rc

  LIS_lai(n)%laiInterval  = 864000
  LIS_lai(n)%laiIntervalType = "10-day" !10-day

  call LIS_registerAlarm("LIS lai read alarm",LIS_rc%ts, &
       LIS_lai(n)%laiInterval,&
       intervalType=LIS_lai(n)%laiIntervalType)
 
  call ESMF_ConfigFindLabel(LIS_config,&
       "ALMIPII LAI data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config, &
       LIS_lai(n)%laifile,rc=rc)
  call LIS_verify(rc,'ALMIPII LAI data directory: not defined')

end subroutine setup_ALMIPIIlai



