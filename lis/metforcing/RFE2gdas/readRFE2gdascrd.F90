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
! 26 MAY 2010; Soni Yatheendradas, Initial LIS version for FEWSNET
! 26 JAN 2011; Clement Alo, Modified for 6-hourly RFE2gdas data
!
! !INTERFACE:
subroutine readRFE2gdascrd()
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_config
  use LIS_logMod,    only : LIS_logunit, LIS_endrun, LIS_verify
  use LIS_timeMgrMod, only : LIS_calendar, LIS_date2time 
  use RFE2gdas_forcingMod, only: RFE2gdas_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to the CPC-hosted 
!  USAID/FEWS-NET RFE2gdas.0 forcing, from the LIS configuration file.
!
!EOP

  implicit none
  integer    :: n, rc
  integer    :: doy
  real       :: gmt
  integer    :: status
  type(ESMF_TimeInterval)  :: modelTimeStep 
  type(ESMF_Time)          :: LISstartTime
  type(ESMF_Time)          :: LISendTime
! ___________________________________________

  call ESMF_ConfigFindLabel(LIS_config,"RFE2gdas forcing directory:",rc=rc)
  call LIS_verify(rc,"readconfig: RFE2gdas forcing directory not in lis.config")
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,RFE2gdas_struc(n)%RFE2gdasDir,rc=rc)
    call LIS_verify(rc,"readconfig: RFE2gdas forcing dir value not correct/given")
    write(LIS_logunit,*)'[INFO] For nest: ', n
    write(LIS_logunit,*)'[INFO] RFE2gdas forcing directory: ',trim(RFE2gdas_struc(n)%RFE2gdasDir)
  enddo

  do n=1,LIS_rc%nnest

     ! Data interval is 6-hourly
     call ESMF_TimeIntervalSet(RFE2gdas_struc(n)%timeStep, &
          s=6*60*60,rc=status)
     call ESMF_TimeIntervalSet(modelTimeStep, &
          s=nint(LIS_rc%nts(n)),rc=status)

     ! Confirm if model run time step is really less than RFE2gdas forcing time step
     IF (modelTimeStep .GT. RFE2gdas_struc(n)%timeStep) THEN 
!     IF (LIS_rc%nts(n) .GT. 6*60*60) THEN 
      ! Confirm if model run time step is really less than RFE2gdas forcing time step:
       WRITE(LIS_logunit,*)'[ERR] Model run time step should be sub-daily!'
       WRITE(LIS_logunit,*)' Program stopping ... '
       CALL LIS_endrun()
     ENDIF 

     ! Start date of RFE2gdas availability: October 30, 2000
     call ESMF_TimeSet(RFE2gdas_struc(n)%startTime, yy=2000, &
           mm = 10, dd = 30, h=0, m = 0, s=0, calendar=LIS_calendar, &
           rc=status)
     call LIS_verify(status, &
                 'readRFE2gdascrd: ESMF_TimeSet RFE2gdas_struc%startTime')
     call LIS_date2time(RFE2gdas_struc(n)%st_real, doy, gmt, &
          2000, 10, 30, 0, 0, 0)
     call ESMF_TimeSet(LISstartTime, yy=LIS_rc%syr, &
           mm = LIS_rc%smo, dd = LIS_rc%sda, h=LIS_rc%shr, &
           m = LIS_rc%smn, s=LIS_rc%sss, calendar=LIS_calendar, &
           rc=status)
     call LIS_verify(status, 'readRFE2gdascrd: ESMF_TimeSet LISstartTime')
     call ESMF_TimeSet(LISendTime, yy=LIS_rc%eyr1, &
           mm = LIS_rc%emo1, dd = LIS_rc%eda1, h=LIS_rc%ehr1, &
           m = LIS_rc%emn1, s=LIS_rc%ess1, calendar=LIS_calendar, &
           rc=status)
     call LIS_verify(status, 'readRFE2gdascrd: ESMF_TimeSet LISendTime')

     ! Confirm if at least some part of the run time duration uses RFE2gdas forcing :
     IF( RFE2gdas_struc(n)%startTime .GE. LISendTime ) THEN 
       WRITE(LIS_logunit,*) &
             '[ERR] No part of run time duration has RFE2gdas availability!'
       WRITE(LIS_logunit,*) ' RFE2gdas availability is from October 30, 2000 onwards.'
       WRITE(LIS_logunit,*) ' Please update your lis.config file ...'
!       WRITE(LIS_logunit,*) 'Program stopping ... '
!       CALL LIS_endrun()
     ENDIF 

   ! Warning that some part of the run-time duration does not use RFE2gdas forcing: 
     IF( LISstartTime .LT. RFE2gdas_struc(n)%startTime ) THEN 
        WRITE(LIS_logunit,*) &
             '[WARN] RFE2gdas not availabile for some part of run time duration'
     ENDIF 

!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure
! data is read in during first time step
!------------------------------------------------------------------------
     RFE2gdas_struc(n)%RFE2gdasEndTime = 0.0

  enddo

end subroutine readRFE2gdascrd

