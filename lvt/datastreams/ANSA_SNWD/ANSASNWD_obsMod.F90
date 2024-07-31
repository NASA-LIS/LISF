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
! !MODULE: ANSASNWD_obsMod
!  \label(ANSASNWD_obsMod)
!
! !INTERFACE:
module ANSASNWD_obsMod
! 
! !USES: 
  use ESMF

  implicit none
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

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ANSASNWD_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ANSASNWDobs
!EOP
  type, public :: ansasnwdobsdec

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
     real,    allocatable     :: snwd(:)
     integer              :: offset1, offset2

  end type ansasnwdobsdec

  type(ansasnwdobsdec), allocatable:: ANSASNWDobs(:)

contains
  
!BOP
! 
! !ROUTINE: ANSASNWD_obsInit
! \label{ANSASNWD_obsInit}
!
! !INTERFACE: 
  subroutine ANSASNWD_obsinit(i)
! 
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_config, LVT_domain
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod,     only : LVT_verify, LVT_logunit

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading VU AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer,   intent(IN) :: i 
!EOP

    real                   :: gridDesci(50)
    real                   :: cornerlat1, cornerlat2
    real                   :: cornerlon1, cornerlon2
    real                   :: minlat, minlon, maxlon, maxlat, dx, dy
    integer                :: status
    type(ESMF_Time)         :: alarmTime
    type(ESMF_TimeInterval) :: alarmInterval

    if(.not.allocated(ANSASNWDobs)) then 
       allocate(ANSASNWDobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, ANSASNWDobs(i)%odir, &
         label='ANSA snow depth observation directory:', rc=status)
    call LVT_verify(status, 'ANSA snow depth observation directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    allocate(ANSASNWDobs(i)%snwd(LVT_rc%lnc*LVT_rc%lnr))

    call ESMF_ConfigFindLabel(LVT_config,"ANSA snow depth lower left lat:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASNWDobs(i)%gridDesc(1),&
         rc=status)
    call LVT_verify(status,'ANSA snow depth lower left lat: not defined')
    
    call ESMF_ConfigFindLabel(LVT_config,"ANSA snow depth lower left lon:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASNWDobs(i)%gridDesc(2),&
         rc=status)
    call LVT_verify(status,'ANSA snow depth lower left lon: not defined')

    
    call ESMF_ConfigFindLabel(LVT_config,"ANSA snow depth upper right lat:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASNWDobs(i)%gridDesc(3),&
         rc=status)
    call LVT_verify(status,'ANSA snow depth upper right lat: not defined')

    call ESMF_ConfigFindLabel(LVT_config,"ANSA snow depth upper right lat:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASNWDobs(i)%gridDesc(3),&
         rc=status)
    call LVT_verify(status,'ANSA snow depth upper right lat: not defined')

    call ESMF_ConfigFindLabel(LVT_config,"ANSA snow depth upper right lon:",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASNWDobs(i)%gridDesc(4),&
         rc=status)
    call LVT_verify(status,'ANSA snow depth upper right lon: not defined')

    call ESMF_ConfigFindLabel(LVT_config,"ANSA snow depth resolution (dx):",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASNWDobs(i)%gridDesc(5),&
         rc=status)
    call LVT_verify(status,'ANSA snow depth resolution (dx): not defined')

    call ESMF_ConfigFindLabel(LVT_config,"ANSA snow depth resolution (dy):",&
         rc=status)
    call ESMF_ConfigGetattribute(LVT_config,ANSASNWDobs(i)%gridDesc(6),&
         rc=status)
    call LVT_verify(status,'ANSA snow depth resolution (dy): not defined')

    minlat = ANSASNWDobs(i)%gridDesc(1)
    minlon = ANSASNWDobs(i)%gridDesc(2)
    maxlat = ANSASNWDobs(i)%gridDesc(3)
    maxlon = ANSASNWDobs(i)%gridDesc(4)
    dx = ANSASNWDobs(i)%gridDesc(5)
    dy = ANSASNWDobs(i)%gridDesc(6)

!sets the local domain corner points with additional buffer 
    cornerlat1 = max(minlat, nint((LVT_domain%minlat-minlat)/dx)*dx+minlat-2*dx)
    cornerlon1 = max(minlon, nint((LVT_domain%minlon-minlon)/dy)*dy+minlon-2*dy)
    cornerlat2 = min(maxlat, nint((LVT_domain%maxlat-minlat)/dx)*dx+minlat+2*dx)
    cornerlon2 = min(maxlon, nint((LVT_domain%maxlon-minlon)/dy)*dy+minlon+2*dy)

    ANSASNWDobs(i)%offset1 = nint((cornerlon1-minlon)/dy)
    ANSASNWDobs(i)%offset2 = nint((cornerlat1-minlat)/dx)
    
    ANSASNWDobs(i)%nr = nint((cornerlat2-cornerlat1)/dx)+1
    ANSASNWDobs(i)%nc = nint((cornerlon2-cornerlon1)/dy)+1
    
    gridDesci(1) = 0 
    gridDesci(2) = ANSASNWDobs(i)%nc
    gridDesci(3) = ANSASNWDobs(i)%nr
    gridDesci(4) = cornerlat1
    gridDesci(5) = cornerlon1
    gridDesci(6) = 128
    gridDesci(7) = cornerlat2
    gridDesci(8) = cornerlon2
    gridDesci(9) = ANSASNWDobs(i)%gridDesc(5)
    gridDesci(10) = ANSASNWDobs(i)%gridDesc(6)
    gridDesci(20) = 64.0
    
    ANSASNWDobs(i)%mi = ANSASNWDobs(i)%nc*ANSASNWDobs(i)%nr
    allocate(ANSASNWDobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASNWDobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(ANSASNWDobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASNWDobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASNWDobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASNWDobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(ANSASNWDobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASNWDobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASNWDobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ANSASNWDobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    ANSASNWDobs(i)%snwd = LVT_rc%udef
    
    call bilinear_interp_input(gridDesci, LVT_rc%gridDesc, & 
         LVT_rc%lnc*LVT_rc%lnr, & 
         ANSASNWDobs(i)%rlat,ANSASNWDobs(i)%rlon,&
         ANSASNWDobs(i)%n11,ANSASNWDobs(i)%n12, &
         ANSASNWDobs(i)%n21,ANSASNWDobs(i)%n22, &
         ANSASNWDobs(i)%w11,ANSASNWDobs(i)%w12, &
         ANSASNWDobs(i)%w21,ANSASNWDobs(i)%w22)

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
    call LVT_verify(status,'alarmset in ANSASNWDsnow_obsMod')
    
    if(LVT_rc%hr.gt.0) then 
       alarmTime = alarmTime + alarmInterval
    endif
    
    ANSASNWDobs(i)%readAlarm = ESMF_AlarmCreate(clock=LVT_clock,&
         ringTime = alarmTime, ringInterval=alarmInterval, &
         rc=status)
    call LVT_verify(status, "ANSASNWDsnow_setup : AlarmCreate")

    ANSASNWDobs(i)%startflag = .true. 
!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------
  end subroutine ANSASNWD_obsinit


end module ANSASNWD_obsMod
