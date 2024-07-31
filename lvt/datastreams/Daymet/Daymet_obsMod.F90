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
! !MODULE: Daymet_obsMod
!  \label(Daymet_obsMod)
!
! !INTERFACE:
module Daymet_obsMod
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
!  This module handles the observation plugin for the Daymet data
!  
! Daymet is a collection of algorithms and computer software designed
!  to interpolate and extrapolate from daily meteorological observations
!  to produce gridded estimates of daily weather parameters.
!  Weather parameters generated include daily surfaces of minimum 
!  and maximum temperature, precipitation, humidity, and radiation 
!  produced on a 1 km x 1 km gridded surface.
! 
! Website: https://daymet.ornl.gov/overview.html
! 
! !NOTES: Currently the plugin only handles the processing of Daymet SWE data
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  29 Sept 2016   Sujay Kumar  Initial Specification
! 
!EOP
!
! 
!
! 
!

  PRIVATE 

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: Daymet_obsinit !Initializes structures for reading Daymet data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: Daymetobs !Object to hold Daymet observation attributes
!EOP
  type, public :: daymetobsdec
     character*100               :: odir
     integer                     :: nc
     integer                     :: nr
     real                        :: udef
     logical                     :: startMode
     type(ESMF_Clock)            :: clock
     type(ESMF_Time)             :: startTime, startMonth
     type(ESMF_TimeInterval)     :: timestep
     real,  allocatable          :: swe(:,:)
     integer, allocatable        :: n11(:)
  end type daymetobsdec

  type(daymetobsdec), allocatable :: daymetobs(:)

contains
  
!BOP
! 
! !ROUTINE: Daymet_obsInit
! \label{Daymet_obsInit}
!
! !INTERFACE: 
  subroutine Daymet_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading Daymet data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer            :: status, rc
    integer            :: ftn, c,r,k
    real*8             :: tdur
    integer            :: syr, smo, sda, shr, smn, sss
    integer            :: ts
    integer            :: t, ios

    if(.not.allocated(daymetobs)) then 
       allocate(daymetobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, daymetobs(i)%odir, &
         label='Daymet observation directory:',rc=status)
    call LVT_verify(status, 'Daymet observation directory: not defined')

    daymetobs(i)%udef = -1
    daymetobs(i)%startMode = .true. 
    syr = 2002
    smo = 1
    sda = 1
    shr = 0
    smn = 0 
    call ESMF_TimeSet(daymetobs(i)%startTime,yy=syr, &
         mm=smo, dd=sda, h=shr, m=smn, &
         s = 0, calendar=LVT_calendar, rc=status)
    call LVT_verify(status, 'daymet: starttime set failed')
    
    ts = 86400
    call ESMF_TimeIntervalSet(daymetobs(i)%timestep, s=ts, rc=status)
    call LVT_verify(status, 'error in setting timestep (daymetobs)')

    call LVT_update_timestep(LVT_rc, 86400)

    daymetobs(i)%nc = 8075
    daymetobs(i)%nr = 7814
    
    allocate(daymetobs(i)%swe(LVT_rc%lnc,LVT_rc%lnr))
    daymetobs(i)%swe = LVT_rc%udef

  end subroutine Daymet_obsinit


end module Daymet_obsMod
