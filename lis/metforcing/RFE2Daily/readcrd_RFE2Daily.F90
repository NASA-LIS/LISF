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
!
! !ROUTINE: readcrd_RFE2Daily
! \label{readcrd_RFE2Daily}
!
! !REVISION HISTORY:
!  26 MAY 2010; Soni Yatheendradas, Initial LIS version for FEWSNET
!  20 Mar 2013; KR Arsenault, Cleaned up code and added adjustable hour offset
!
! !INTERFACE:
subroutine readcrd_RFE2Daily()
! !USES:
  use ESMF
  use LIS_coreMod,         only : LIS_rc, LIS_config
  use LIS_logMod,          only : LIS_logunit, LIS_endrun, LIS_verify
  use LIS_timeMgrMod,      only : LIS_calendar, LIS_date2time 
  use RFE2Daily_forcingMod, only: RFE2Daily_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to the CPC-hosted 
!  USAID/FEWS-NET RFE2.0 forcing, from the LIS configuration file.
!
!EOP

  implicit none
  integer :: n, rc
  integer :: doy
  real    :: gmt
  integer :: status
  type(ESMF_TimeInterval) :: modelTimeStep 
  type(ESMF_Time)         :: LISstartTime
  type(ESMF_Time)         :: LISendTime
! _________________________________________________

  RFE2Daily_struc%hour_offset = 6

! Read in RFE2 Daily PPT directory:
  call ESMF_ConfigFindLabel(LIS_config,"RFE2Daily forcing directory:",rc=rc)
  call LIS_verify(rc,"LISconfig: RFE2Daily forcing directory not in lis.config")
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,RFE2Daily_struc(n)%RFE2DailyDir,rc=rc)
     call LIS_verify(rc,"LISconfig: RFE2Daily forcing dir value not correct/given")
     write(LIS_logunit,*) 'For nest ', n, ', RFE2Daily forcing directory: ', &
                             trim(RFE2Daily_struc(n)%RFE2DailyDir)
  enddo
  write(LIS_logunit,*) 'Using RFE2Daily forcing'

! Read in RFE2 Daily PPT time-based hour-offset (should be 6Z; but 0Z for WRSI model):
  call ESMF_ConfigFindLabel(LIS_config, "RFE2Daily time offset:",rc=rc) 
  call LIS_verify(rc,"LISconfig: RFE2Daily time offset not in lis.config")
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,RFE2Daily_struc(n)%hour_offset,rc=rc)
     write(LIS_logunit,*) 'For nest ', n, ', RFE2Daily hour offset: ', &
                             RFE2Daily_struc(n)%hour_offset
  enddo

  do n=1,LIS_rc%nnest

   ! Data interval is daily
     call ESMF_TimeIntervalSet(RFE2Daily_struc(n)%timeStep, &
                               s=24*60*60,rc=status)
     call ESMF_TimeIntervalSet(modelTimeStep, &
                               s=nint(LIS_rc%nts(n)),rc=status)
     
     IF (modelTimeStep .GT. RFE2Daily_struc(n)%timeStep) THEN
     ! Confirm if model run time step is really less than RFE2Daily forcing time step
       WRITE(LIS_logunit,*) 'Model run time step should be sub-daily!!'
       WRITE(LIS_logunit,*) 'Program stopping ... '
       CALL LIS_endrun()
     ENDIF 

     ! Start date of RFE2Daily availability: October 30, 2000 ! SY: (Changed to 
     ! reflect data representative period of 0600 - 0600 GMT as per Ronald W Lietzow & 
     ! James Rowland, USGS, E-mail Communication 01/07/2011)
     call ESMF_TimeSet(RFE2Daily_struc(n)%startTime, yy=2000, &
           mm = 10, dd = 30, h=RFE2Daily_struc(n)%hour_offset, &
           m = 0, s=0, calendar=LIS_calendar, rc=status)
     call LIS_verify(status, &
                 'readcrd_RFE2Daily: ESMF_TimeSet RFE2Daily_struc%startTime')
     ! SY: Changed here also to reflect 0600 - 0600 GMT as for RFE2Daily_struc(n)%startTime above
     call LIS_date2time(RFE2Daily_struc(n)%st_real, doy, gmt, &
          2000, 10, 30, RFE2Daily_struc(n)%hour_offset, 0, 0)
     call ESMF_TimeSet(LISstartTime, yy=LIS_rc%syr, &
           mm = LIS_rc%smo, dd = LIS_rc%sda, h=LIS_rc%shr, &
           m = LIS_rc%smn, s=LIS_rc%sss, calendar=LIS_calendar, &
           rc=status)
     call LIS_verify(status, 'readcrd_RFE2Daily: ESMF_TimeSet LISstartTime')
     call ESMF_TimeSet(LISendTime, yy=LIS_rc%eyr1, &
           mm = LIS_rc%emo1, dd = LIS_rc%eda1, h=LIS_rc%ehr1, &
           m = LIS_rc%emn1, s=LIS_rc%ess1, calendar=LIS_calendar, &
           rc=status)
     call LIS_verify(status, 'readcrd_RFE2Daily: ESMF_TimeSet LISendTime')
     IF (RFE2Daily_struc(n)%startTime .GE. LISendTime) THEN
     ! Confirm if at least some part of the run time duration uses RFE2Daily forcing
       WRITE(LIS_logunit,*) &
                  'No part of run time duration has RFE2Daily availability!!'
       WRITE(LIS_logunit,*) 'RFE2Daily availability is October 30, 2000 0600 GMT onwards'
       WRITE(LIS_logunit,*) 'Please modify lis.config'
       WRITE(LIS_logunit,*) 'Program stopping ... '
       CALL LIS_endrun()
     ENDIF 

     IF (LISstartTime .LT. RFE2Daily_struc(n)%startTime) THEN
     ! Give warning that some part of the run time duration does not use RFE2Daily forcing 
       WRITE(LIS_logunit,*) &
             'NOTE: RFE2Daily data not available till October 30, 2000 0600 GMT'
     ENDIF 

!------------------------------------------------------------------------
! Setting global observed precip times to LIS_rc%udef to ensure
! data is read in during first time step
!------------------------------------------------------------------------
     RFE2Daily_struc(n)%RFE2DailyEndTime = LIS_rc%udef

  enddo

end subroutine readcrd_RFE2Daily

