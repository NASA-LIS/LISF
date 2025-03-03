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
!  26 MAY 2010; Soni Yatheendradas, Initial LDT version for FEWSNET
!  20 Mar 2013; KR Arsenault, Cleaned up code and added adjustable hour offset
!
! !INTERFACE:
subroutine readcrd_RFE2Daily()

! !USES:
  use ESMF
  use LDT_coreMod,         only : LDT_rc, LDT_config
  use LDT_logMod,          only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_timeMgrMod,      only : LDT_calendar, LDT_date2time 
  use RFE2Daily_forcingMod,only : RFE2Daily_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to the CPC-hosted 
!  USAID/FEWS-NET RFE2.0 forcing, from the LDT configuration file.
!
!EOP

  implicit none
  integer :: n, rc
  integer :: doy
  real    :: gmt
  integer :: status
  type(ESMF_TimeInterval) :: modelTimeStep 
  type(ESMF_Time)         :: LDTstartTime
  type(ESMF_Time)         :: LDTendTime
! _________________________________________________

  RFE2Daily_struc%hour_offset = 6

  write(LDT_logunit,*)" Using RFE2 Daily forcing"

! - Read config file entries:
! Read in RFE2 Daily PPT directory:
  call ESMF_ConfigFindLabel(LDT_config,"RFE2Daily forcing directory:",rc=rc)
  call LDT_verify(rc,"LDTconfig: RFE2Daily forcing directory not in ldt.config")
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,RFE2Daily_struc(n)%RFE2DailyDir,rc=rc)
     call LDT_verify(rc,"LDTconfig: RFE2Daily forcing dir value not correct/given")
  enddo
  write(LDT_logunit,*) " RFE2Daily forcing directory: ", &
        trim(RFE2Daily_struc(1)%RFE2DailyDir)

! Read in RFE2 Daily PPT time-based hour-offset (should be 6Z; but 0Z for WRSI model):
  call ESMF_ConfigFindLabel(LDT_config, "RFE2Daily time offset:",rc=rc) 
  call LDT_verify(rc,"LDTconfig: RFE2Daily time offset not in ldt.config")
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,RFE2Daily_struc(n)%hour_offset,rc=rc)
     write(LDT_logunit,*) "For nest,", n,", RFE2Daily hour offset:", &
           RFE2Daily_struc(n)%hour_offset
  enddo
! ---

!- Set time interval and check available dates of RFE2 data:
  do n=1,LDT_rc%nnest

   ! Data interval is daily
     call ESMF_TimeIntervalSet(RFE2Daily_struc(n)%timeStep, &
                               s=24*60*60,rc=status)

     call ESMF_TimeIntervalSet(modelTimeStep, &
                               s=nint(LDT_rc%nts(n)),rc=status)

     IF( modelTimeStep .GT. RFE2Daily_struc(n)%timeStep ) THEN
     ! Confirm if model run time step is really less than RFE2Daily forcing time step
       WRITE(LDT_logunit,*) "Model run time step should be sub-daily!!"
       WRITE(LDT_logunit,*) "Program stopping ... "
       CALL LDT_endrun()
     ENDIF 

     ! Start date of RFE2Daily availability: October 30, 2000 ! SY: (Changed to 
     ! reflect data representative period of 0600 - 0600 GMT as per Ronald W Lietzow & 
     ! James Rowland, USGS, E-mail Communication 01/07/2011)
     call ESMF_TimeSet(RFE2Daily_struc(n)%startTime, yy=2000, &
           mm=10, dd=30, h=RFE2Daily_struc(n)%hour_offset, &
           m=0, s=0, calendar=LDT_calendar, rc=status)

!     call ESMF_TimePrint(RFE2Daily_struc(n)%startTime)
!     call ESMF_CalendarPrint(LDT_calendar)

     call LDT_verify(status, &
             "readcrd_RFE2Daily: ESMF_TimeSet RFE2Daily_struc%startTime")

   ! SY: Changed here to reflect 0600 - 0600 GMT as for RFE2Daily_struc(n)%startTime above
     call LDT_date2time(RFE2Daily_struc(n)%st_real, doy, gmt, &
          2000, 10, 30, RFE2Daily_struc(n)%hour_offset, 0, 0)
     call ESMF_TimeSet(LDTstartTime, yy=LDT_rc%syr, &
           mm = LDT_rc%smo, dd = LDT_rc%sda, h=LDT_rc%shr, &
           m = LDT_rc%smn, s=LDT_rc%sss, calendar=LDT_calendar, &
           rc=status)
     call LDT_verify(status,"readcrd_RFE2Daily: ESMF_TimeSet LDTstartTime")

     call ESMF_TimeSet(LDTendTime, yy=LDT_rc%eyr1, &
           mm = LDT_rc%emo1, dd = LDT_rc%eda1, h=LDT_rc%ehr1, &
           m = LDT_rc%emn1, s=LDT_rc%ess1, calendar=LDT_calendar, &
           rc=status)
     call LDT_verify(status, 'readcrd_RFE2Daily: ESMF_TimeSet LDTendTime')

     IF( RFE2Daily_struc(n)%startTime .GE. LDTendTime ) THEN
     ! Confirm if at least some part of the run time duration uses RFE2Daily forcing
       WRITE(LDT_logunit,*) &
                  'No part of run time duration has RFE2 Daily availability!!'
       WRITE(LDT_logunit,*) 'RFE2Daily availability is October 30, 2000, 0600 GMT onwards'
       WRITE(LDT_logunit,*) 'Please modify ldt.config'
       WRITE(LDT_logunit,*) 'Program stopping ... '
       CALL LDT_endrun()
     ENDIF 

     IF( LDTstartTime .LT. RFE2Daily_struc(n)%startTime ) THEN
     ! Give warning that some part of the run time duration does not use RFE2Daily forcing 
       WRITE(LDT_logunit,*) &
            "* NOTE: RFE2Daily data not available till October 30, 2000 0600 GMT "
     ENDIF 

!------------------------------------------------------------------------
! Setting global observed precip times to LDT_rc%udef to ensure
! data is read in during first time step
!------------------------------------------------------------------------
     RFE2Daily_struc(n)%RFE2DailyEndTime = LDT_rc%udef

  enddo

end subroutine readcrd_RFE2Daily

