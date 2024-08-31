!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine setup_NESDISgfrac(n)
  use ESMF
  use LIS_coreMod
  use LIS_logMod,     only : LIS_verify
  use LIS_vegDataMod, only : LIS_gfrac
  use LIS_timeMgrMod

  integer, intent(in)  :: n 

  integer              :: rc
  real :: gridDesci(50)

  LIS_gfrac(n)%gfracIntervalType = "weekly" !weekly
  LIS_gfrac(n)%gfracInterval = 86400 ! 7*86400
  !
  ! This is the time that NESDIS data starts. 
  !
!  call ESMF_TimeSet(refTime, yy=1981,mm=8,dd=24, &
!       calendar=LIS_Calendar,rc=rc)
!  call LIS_verify(rc,'reftimeset:setup_nesdisgfrac')

  ! Attach the reference time to the gfrac alarm

  gridDesci = 0
  gridDesci(1) = 0
  gridDesci(2) = 2500
  gridDesci(3) = 1250
  gridDesci(4) = -89.9280
  gridDesci(5) = -179.9280
  gridDesci(6) = 128         
  gridDesci(7) = 89.9280
  gridDesci(8) = 179.9280
  gridDesci(9) = 0.144
  gridDesci(10) = 0.144
  gridDesci(20) = 64

  allocate(LIS_gfrac(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  allocate(LIS_gfrac(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  ! EMK...Fixed argument list
  call bilinear_interp_input(n, gridDesci, &
       LIS_gfrac(n)%n111, LIS_gfrac(n)%n121, &
       LIS_gfrac(n)%n211, LIS_gfrac(n)%n221, &
       LIS_gfrac(n)%w111, LIS_gfrac(n)%w121, &
       LIS_gfrac(n)%w211, LIS_gfrac(n)%w221 )

  call LIS_registerAlarm("LIS gfrac read alarm",LIS_rc%ts, &
       LIS_gfrac(n)%gfracInterval,&
       intervalType=LIS_gfrac(n)%gfracIntervalType) 


  call ESMF_ConfigFindLabel(LIS_config,&
       "NESDIS greenness data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config, &
       LIS_gfrac(n)%gfracfile,rc=rc)
  call LIS_verify(rc,'NESDIS greenness data directory: not defined')


end subroutine setup_NESDISgfrac



