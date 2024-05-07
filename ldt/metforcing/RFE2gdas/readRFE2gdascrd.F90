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
! !ROUTINE: readRFE2gdascrd
! \label{readRFE2gdascrd}
!
! !REVISION HISTORY:
! 26 MAY 2010; Soni Yatheendradas, Initial LDT version for FEWSNET
! 26 JAN 2011; Clement Alo, Modified for 6-hourly RFE2gdas data
!
! !INTERFACE:
subroutine readRFE2gdascrd()
! !USES:
  use ESMF
  use LDT_coreMod,   only : LDT_rc, LDT_config
  use LDT_logMod,    only : LDT_logunit, LDT_endrun, LDT_verify
  use LDT_timeMgrMod, only : LDT_calendar, LDT_date2time 
  use RFE2gdas_forcingMod, only: RFE2gdas_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to the CPC-hosted 
!  USAID/FEWS-NET RFE2gdas.0 forcing, from the LDT configuration file.
!
!EOP

  implicit none
  integer  :: n, rc
  integer  :: doy
  real     :: gmt
  integer  :: status
  type(ESMF_TimeInterval)  :: modelTimeStep 
  type(ESMF_Time)          :: LDTstartTime
  type(ESMF_Time)          :: LDTendTime

  call ESMF_ConfigFindLabel(LDT_config,"RFE2gdas forcing directory:",rc=rc)
  call LDT_verify(rc,"readconfig: RFE2gdas forcing directory not in lis.config")
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,RFE2gdas_struc(n)%RFE2gdasDir,rc=rc)
     call LDT_verify(rc,"readconfig: RFE2gdas forcing dir value not correct/given")
   write(LDT_logunit,*) 'For nest: ', n
   write(LDT_logunit,*)' RFE2gdas forcing directory: ', RFE2gdas_struc(n)%RFE2gdasDir
  enddo

  do n=1,LDT_rc%nnest

!    data interval is 6-hourly
     call ESMF_TimeIntervalSet(RFE2gdas_struc(n)%timeStep, &
          s=6*60*60,rc=status)
     call ESMF_TimeIntervalSet(modelTimeStep, &
          s=nint(LDT_rc%nts(n)),rc=status)

     IF (modelTimeStep .GT. RFE2gdas_struc(n)%timeStep) THEN 
     ! Confirm if model run time step is really less than RFE2gdas forcing time step
!     IF (LDT_rc%nts(n) .GT. 6*60*60) THEN ! Confirm if model run time step is really less than RFE2gdas forcing time step
       WRITE(LDT_logunit,*) '[ERR] Model run time step should be sub-daily!!'
       WRITE(LDT_logunit,*) ' Program stopping ... '
       CALL LDT_endrun()
     ENDIF 

     ! Start date of RFE2gdas availability: October 30, 2000
     call ESMF_TimeSet(RFE2gdas_struc(n)%startTime, yy=2000, &
           mm = 10, dd = 30, h=0, m = 0, s=0, calendar=LDT_calendar, &
           rc=status)
     call LDT_verify(status, &
                 'readRFE2gdascrd: ESMF_TimeSet RFE2gdas_struc%startTime')
     call LDT_date2time(RFE2gdas_struc(n)%st_real, doy, gmt, &
          2000, 10, 30, 0, 0, 0)
     call ESMF_TimeSet(LDTstartTime, yy=LDT_rc%syr, &
           mm = LDT_rc%smo, dd = LDT_rc%sda, h=LDT_rc%shr, &
           m = LDT_rc%smn, s=LDT_rc%sss, calendar=LDT_calendar, &
           rc=status)
     call LDT_verify(status, 'readRFE2gdascrd: ESMF_TimeSet LDTstartTime')
     call ESMF_TimeSet(LDTendTime, yy=LDT_rc%eyr1, &
           mm = LDT_rc%emo1, dd = LDT_rc%eda1, h=LDT_rc%ehr1, &
           m = LDT_rc%emn1, s=LDT_rc%ess1, calendar=LDT_calendar, &
           rc=status)
     call LDT_verify(status, 'readRFE2gdascrd: ESMF_TimeSet LDTendTime')
     IF (RFE2gdas_struc(n)%startTime .GE. LDTendTime) THEN 
     ! Confirm if atleast some part of the run time duration uses RFE2gdas forcing 
       WRITE(LDT_logunit,*) &
                  '[ERR] No part of run time duration has RFE2gdas availability!!'
       WRITE(LDT_logunit,*) '  RFE2gdas availability is October 30, 2000 onwards'
       WRITE(LDT_logunit,*) '  Please modify lis.config'
       WRITE(LDT_logunit,*) ' Program stopping ... '
       CALL LDT_endrun()
     ENDIF 

     IF (LDTstartTime .LT. RFE2gdas_struc(n)%startTime) THEN 
    ! Give warning/notification that some part of the run time duration does not use RFE2gdas forcing 
       WRITE(LDT_logunit,*) &
            '[WARN] RFE2gdas not availabile for some part of run time duration'
     ENDIF 

!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure
! data is read in during first time step
!------------------------------------------------------------------------
     RFE2gdas_struc(n)%RFE2gdasEndTime = 0.0

  enddo

end subroutine readRFE2gdascrd

