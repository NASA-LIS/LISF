!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: NASA_AMSREsm_obsMod
! 
! !DESCRIPTION: 
!  This module handles the observation plugin for the standard NASA
!  AMSR-E soil moisture retrieval product
!  
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
!
module NASA_AMSREsm_obsMod
! !USES: 
  use ESMF
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: NASA_AMSREsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: NASA_AMSREsmobs
!EOP
  type, public :: nasaamsresmobsdec

     character(len=LDT_CONST_PATH_LEN)        :: odir

     integer         :: mo
     integer,allocatable :: n112(:)
     real            :: gridDesci(20)

     real,    allocatable    :: smobs(:,:)
     real,    allocatable    :: smtime(:,:)
     integer*2, allocatable  :: smqc(:,:)     
     logical             :: startflag     

  end type nasaamsresmobsdec

  type(nasaamsresmobsdec),save:: NASA_AMSREsmobs

contains
  
!BOP
! 
! !ROUTINE: NASA_AMSREsm_obsInit
! \label{NASA_AMSREsm_obsInit}
! 
! !INTERFACE: 
  subroutine NASA_AMSREsm_obsinit()
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading NASA AMSRE soil moisture data. 
! 
!EOP
    integer         ezlh_convert
    integer            :: npts
    integer,parameter  :: ease_nr=586
    integer,parameter  :: ease_nc=1383
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    integer                 :: n 

    n = 1
    
    call ESMF_ConfigGetAttribute(LDT_Config, NASA_AMSREsmobs%odir, &
         label='NASA AMSRE soil moisture observation directory:', rc=status)
    call LDT_verify(status, 'NASA AMSRE soil moisture observation directory: not defined')

    NASA_AMSREsmobs%gridDesci=0.0 
    
    !filling the items needed by the interpolation library
    NASA_AMSREsmobs%gridDesci(1) = 9  !input is EASE grid
    !these  corner coordinates were calculated based on ezlh_convert
    NASA_AMSREsmobs%gridDesci(4) = -90.0  !lat
    NASA_AMSREsmobs%gridDesci(5) = -179.6096 !lon
    NASA_AMSREsmobs%gridDesci(7) = 83.33788  !lat
    NASA_AMSREsmobs%gridDesci(8) = 180.1301  !lon
    
    
    NASA_AMSREsmobs%gridDesci(2) = ease_nc  !nx
    NASA_AMSREsmobs%gridDesci(3) = ease_nr  !ny
    
    NASA_AMSREsmobs%gridDesci(9)  = 1  !Ml

    npts= LDT_rc%lnc(n)*LDT_rc%lnr(n)
    NASA_AMSREsmobs%mo=npts
    allocate(NASA_AMSREsmobs%n112(npts))
    
    NASA_AMSREsmobs%n112=0.0
    call neighbor_interp_input(n, NASA_AMSREsmobs%gridDesci, &
         NASA_AMSREsmobs%n112)

    NASA_AMSREsmobs%startflag = .true. 
#if 0 
    call ESMF_TimeIntervalSet(alarmInterval, &
         d = 1, h = 0, m = 0, s = 0, &
         calendar = LDT_calendar, rc=status)

    call ESMF_TimeSet(alarmTime, yy = LDT_rc%yr, &
         mm = LDT_rc%mo, &
         dd = LDT_rc%da, &
         h  = 0, &
         m  = 0, & 
         s  = 0, &
         calendar = LDT_calendar,&
         rc = status)
    call LDT_verify(status,'alarmset in NASA_AMSREsmsnow_obsMod')

    if(LDT_rc%hr.gt.0) then 
       alarmTime = alarmTime + alarmInterval
    endif
    print*, LDT_rc%yr, LDT_rc%mo, LDT_rc%da, LDT_rc%hr, LDT_rc%mn, LDT_rc%ss
    
    NASA_AMSREsmobs%readAlarm = ESMF_AlarmCreate(clock=LDT_clock,&
         ringTime = alarmTime, ringInterval=alarmInterval, &
         enabled = .true., rc=status)

    call LDT_verify(status, "NASA_AMSREsm_setup : AlarmCreate")
#endif    
    allocate(NASA_AMSREsmobs%smobs(LDT_rc%lnc(n)*LDT_rc%lnr(n),2))
    allocate(NASA_AMSREsmobs%smtime(LDT_rc%lnc(n)*LDT_rc%lnr(n),2))
    allocate(NASA_AMSREsmobs%smqc(LDT_rc%lnc(n)*LDT_rc%lnr(n),2))

!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%soilmoist_obs, "m3/m3",1,1)
    LDT_DAobsData(1)%soilmoist_obs%selectStats = 1

  end subroutine NASA_AMSREsm_obsinit


end module NASA_AMSREsm_obsMod
