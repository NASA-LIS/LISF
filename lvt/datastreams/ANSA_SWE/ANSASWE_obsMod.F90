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
! !MODULE: ANSASWE_obsMod
!  \label(ANSASWE_obsMod)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the standard VU
!  AMSR-E soil moisture retrieval product
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
!EOP
! 
!  
! 
!
module ANSASWE_obsMod
! !USES: 
  use ESMF

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ANSASWE_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ANSASWEobs
!EOP
  type, public :: ansasweobsdec

     character*100        :: odir
     real                 :: gridDesc(6)
     integer              :: mi
     integer              :: nc, nr
     type(ESMF_Alarm)     :: readAlarm
     logical              :: startflag
     real, allocatable        :: rlat(:)
     real, allocatable        :: rlon(:)
     integer, allocatable     :: n11(:)
     integer, allocatable     :: n12(:)
     integer, allocatable     :: n21(:)
     integer, allocatable     :: n22(:)
     real, allocatable        :: w11(:)
     real, allocatable        :: w12(:)
     real, allocatable        :: w21(:)
     real, allocatable        :: w22(:)
     real,    allocatable     :: swe(:)
     integer              :: offset1, offset2

  end type ansasweobsdec

  type(ansasweobsdec), allocatable:: ANSASWEobs(:)

contains
  
!BOP
! 
! !ROUTINE: ANSASWE_obsInit
! \label{ANSASWE_obsInit}
! 
! !INTERFACE: 
  subroutine ANSASWE_obsinit(i)
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_config, LVT_domain
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod,     only : LVT_verify, LVT_logunit

    implicit none
! !ARGUMENTS: 
    integer,   intent(IN) :: i 
! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading VU AMSRE soil moisture data. 
! 
!EOP

    real                   :: gridDesci(50)
    real                   :: cornerlat1, cornerlat2
    real                   :: cornerlon1, cornerlon2
    real                   :: minlat, minlon, maxlon, maxlat, dx, dy
    integer                :: status
    type(ESMF_Time)         :: alarmTime
    type(ESMF_TimeInterval) :: alarmInterval

    if(.not.allocated(ANSASWEobs)) then 
       allocate(ANSASWEobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, ANSASWEobs(i)%odir, &
         label='ANSA SWE observation directory:', rc=status)
    call LVT_verify(status, 'ANSA SWE observation directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    allocate(ANSASWEobs(i)%swe(LVT_rc%lnc*LVT_rc%lnr))

    call ESMF_ConfigFindLabel(LVT_config,"ANSA SWE lower left lat:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASWEobs(i)%gridDesc(1),&
         rc=status)
    call LVT_verify(status,'ANSA SWE lower left lat: not defined')
    
    call ESMF_ConfigFindLabel(LVT_config,"ANSA SWE lower left lon:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASWEobs(i)%gridDesc(2),&
         rc=status)
    call LVT_verify(status,'ANSA SWE lower left lon: not defined')

    
    call ESMF_ConfigFindLabel(LVT_config,"ANSA SWE upper right lat:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASWEobs(i)%gridDesc(3),&
         rc=status)
    call LVT_verify(status,'ANSA SWE upper right lat: not defined')

    call ESMF_ConfigFindLabel(LVT_config,"ANSA SWE upper right lat:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASWEobs(i)%gridDesc(3),&
         rc=status)
    call LVT_verify(status,'ANSA SWE upper right lat: not defined')

    call ESMF_ConfigFindLabel(LVT_config,"ANSA SWE upper right lon:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASWEobs(i)%gridDesc(4),&
         rc=status)
    call LVT_verify(status,'ANSA SWE upper right lon: not defined')

    call ESMF_ConfigFindLabel(LVT_config,"ANSA SWE resolution (dx):",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASWEobs(i)%gridDesc(5),&
         rc=status)
    call LVT_verify(status,'ANSA SWE resolution (dx): not defined')

    call ESMF_ConfigFindLabel(LVT_config,"ANSA SWE resolution (dy):",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASWEobs(i)%gridDesc(6),&
         rc=status)
    call LVT_verify(status,'ANSA SWE resolution (dy): not defined')

    minlat = ANSASWEobs(i)%gridDesc(1)
    minlon = ANSASWEobs(i)%gridDesc(2)
    maxlat = ANSASWEobs(i)%gridDesc(3)
    maxlon = ANSASWEobs(i)%gridDesc(4)
    dx = ANSASWEobs(i)%gridDesc(5)
    dy = ANSASWEobs(i)%gridDesc(6)

!sets the local domain corner points with additional buffer 
    cornerlat1 = max(minlat, nint((LVT_domain%minlat-minlat)/dx)*dx+minlat-2*dx)
    cornerlon1 = max(minlon, nint((LVT_domain%minlon-minlon)/dy)*dy+minlon-2*dy)
    cornerlat2 = min(maxlat, nint((LVT_domain%maxlat-minlat)/dx)*dx+minlat+2*dx)
    cornerlon2 = min(maxlon, nint((LVT_domain%maxlon-minlon)/dy)*dy+minlon+2*dy)

    ANSASWEobs(i)%offset1 = nint((cornerlon1-minlon)/dy)
    ANSASWEobs(i)%offset2 = nint((cornerlat1-minlat)/dx)
    
    ANSASWEobs(i)%nr = nint((cornerlat2-cornerlat1)/dx)+1
    ANSASWEobs(i)%nc = nint((cornerlon2-cornerlon1)/dy)+1
    
    gridDesci(1) = 0 
    gridDesci(2) = ANSASWEobs(i)%nc
    gridDesci(3) = ANSASWEobs(i)%nr
    gridDesci(4) = cornerlat1
    gridDesci(5) = cornerlon1
    gridDesci(6) = 128
    gridDesci(7) = cornerlat2
    gridDesci(8) = cornerlon2
    gridDesci(9) = ANSASWEobs(i)%gridDesc(5)
    gridDesci(10) = ANSASWEobs(i)%gridDesc(6)
    gridDesci(20) = 64.0
    
    ANSASWEobs(i)%mi = ANSASWEobs(i)%nc*ANSASWEobs(i)%nr
    allocate(ANSASWEobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASWEobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(ANSASWEobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASWEobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASWEobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASWEobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(ANSASWEobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASWEobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASWEobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASWEobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(ANSASWEobs(i)%swe(LVT_rc%lnc*LVT_rc%lnr))
    
    ANSASWEobs(i)%swe = LVT_rc%udef
    
    call bilinear_interp_input(gridDesci, LVT_rc%gridDesc, & 
         LVT_rc%lnc*LVT_rc%lnr, & 
         ANSASWEobs(i)%rlat,ANSASWEobs(i)%rlon,&
         ANSASWEobs(i)%n11,ANSASWEobs(i)%n12, &
         ANSASWEobs(i)%n21,ANSASWEobs(i)%n22, &
         ANSASWEobs(i)%w11,ANSASWEobs(i)%w12, &
         ANSASWEobs(i)%w21,ANSASWEobs(i)%w22)

!--------------------------------------------------------------------------------
! The assimilation interval will be every 24 hours, at 10.30 localtime
! The data will be read and kept in memory at 0z
!--------------------------------------------------------------------------------
    call ESMF_TimeIntervalSet(alarmInterval, &
         d = 1, h = 0, m = 0, s = 0, &
         calendar = LVT_calendar, rc=status)
    
    call ESMF_TimeSet(alarmTime, yy = LVT_rc%yr, &
         mm = LVT_rc%mo, &
         dd = LVT_rc%da, &
         h  = 0, &
         m  = 0, & 
         s  = 0, &
         calendar = LVT_calendar,&
         rc = status)
    call LVT_verify(status,'alarmset in ANSASWEsnow_obsMod')
    
    if(LVT_rc%hr.gt.0) then 
       alarmTime = alarmTime + alarmInterval
    endif
    
    ANSASWEobs(i)%readAlarm = ESMF_AlarmCreate(clock=LVT_clock,&
         ringTime = alarmTime, ringInterval=alarmInterval, &
         rc=status)
    call LVT_verify(status, "ANSASWEsnow_setup : AlarmCreate")

    ANSASWEobs(i)%startflag = .true. 
!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------

  end subroutine ANSASWE_obsinit


end module ANSASWE_obsMod
