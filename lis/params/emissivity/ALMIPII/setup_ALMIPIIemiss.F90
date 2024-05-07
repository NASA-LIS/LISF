!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine setup_ALMIPIIemiss(n)

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_emissMod,   only : LIS_emiss

  integer, intent(in)  :: n 
  integer              :: rc

  LIS_emiss(n)%emissInterval  = 864000
  LIS_emiss(n)%emissIntervalType = "10-day" !10-day

  call LIS_registerAlarm("LIS emiss read alarm",LIS_rc%ts, &
       LIS_emiss(n)%emissInterval,&
       intervalType=LIS_emiss(n)%emissIntervalType) 

  call ESMF_ConfigFindLabel(LIS_config,&
       "ALMIPII emissivity data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config, &
       LIS_emiss(n)%emissfile,rc=rc)
  call LIS_verify(rc,'ALMIPII emissivity data directory: not defined')
end subroutine setup_ALMIPIIemiss



