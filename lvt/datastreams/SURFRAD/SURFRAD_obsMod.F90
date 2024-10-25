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
! !MODULE: SURFRAD_obsMod
! \label(SURFRAD_obsMod)
!
! !INTERFACE:
module SURFRAD_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the NOAA Earth 
!  System Research Laboratory Global Monitoring Division's 
!  Surface Radiation (SURFRAD) network data. The dataset 
!  includes observations from seven stations operating in 
!  climatologically diverse regions: Montana, Colorado, 
!  Illinois, Mississippi, Pennsylvania, Nevada, and South 
!  Dakota. The data is available at a minute interval and
!  prior to 2009 at 3 minute intervals. 
!  
!  http://www.srrb.noaa.gov/surfrad/
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Apr 2009   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SURFRAD_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SURFRADobs
!EOP

  type, public :: SURFRADobsdec
     character*100 :: odir
     real, allocatable :: dw_psp(:,:)
     real, allocatable :: windspd(:,:)
     real, allocatable :: dw_pir(:,:)
     real, allocatable :: pres(:,:)
     
     type(ESMF_Time)	:: starttime
     type(ESMF_TimeInterval)	:: timestep
     real, allocatable :: stnlat(:), stnlon(:)
     integer :: day
     logical :: startflag
     integer :: num_stations = 7
     integer :: num_lines = 1441
  end type SURFRADobsdec

  type(SURFRADobsdec), allocatable :: SURFRADobs(:)

contains
  
!BOP
! 
! !ROUTINE: SURFRAD_obsInit
! \label{SURFRAD_obsInit}
!
! !INTERFACE: 
  subroutine SURFRAD_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_logMod
    use LVT_timeMgrMod
    use LVT_histDataMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading the SURFRAD data.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer :: status

    if(.not.allocated(SURFRADobs)) then 
       allocate(SURFRADobs(LVT_rc%nDataStreams))
    endif

!------------------------------------------------------------------------------
! Read any runtime specifications from the lvt.config file. 
!------------------------------------------------------------------------------

    call ESMF_ConfigGetAttribute(LVT_config, SURFRADobs(i)%odir, &
         label='SURFRAD observation directory:',rc=status)
    call LVT_verify(status, 'SURFRAD observation directory: not defined')
  
! 3 minute intervals - switches to 1minute in 2009. 
    call ESMF_TimeIntervalSet(SURFRADobs(i)%timestep, s=3*60, rc=status)
    call LVT_verify(status, 'error in setting timestep (SURFRADobs)')

    if(LVT_rc%yr.ge.2009.and.LVT_rc%mo.ge.1.and.LVT_rc%da.ge.1) then 
       call ESMF_TimeIntervalSet(SURFRADobs(i)%timestep, s=1*60, rc=status)
       call LVT_verify(status, 'error in setting timestep (SURFRADobs)')
    endif
!------------------------------------------------------------------------------
! Initialize any other variables.   
!------------------------------------------------------------------------------

    call LVT_update_timestep(LVT_rc, 60)

    allocate (SURFRADobs(i)%dw_psp(SURFRADobs(i)%num_stations,SURFRADobs(i)%num_lines))
    SURFRADobs(i)%dw_psp = LVT_rc%udef
    allocate (SURFRADobs(i)%dw_pir(SURFRADobs(i)%num_stations,SURFRADobs(i)%num_lines))
    SURFRADobs(i)%dw_pir = LVT_rc%udef
    allocate (SURFRADobs(i)%windspd(SURFRADobs(i)%num_stations,&
         SURFRADobs(i)%num_lines))
    SURFRADobs(i)%windspd = LVT_rc%udef
    allocate (SURFRADobs(i)%pres(SURFRADobs(i)%num_stations,&
         SURFRADobs(i)%num_lines))
    SURFRADobs(i)%pres = LVT_rc%udef

    allocate(SURFRADobs(i)%stnlat(SURFRADobs(i)%num_stations))
    allocate(SURFRADobs(i)%stnlon(SURFRADobs(i)%num_stations))
    
    SURFRADobs(i)%day = -1
    
    SURFRADobs(i)%stnlat(1) = 40.05
    SURFRADobs(i)%stnlon(1) = -88.37
    
    SURFRADobs(i)%stnlat(2) = 40.125
    SURFRADobs(i)%stnlon(2) =  -105.237
    
    SURFRADobs(i)%stnlat(3) = 36.624
    SURFRADobs(i)%stnlon(3) = -116.019
    
    SURFRADobs(i)%stnlat(4) = 48.31
    SURFRADobs(i)%stnlon(4) = -105.1
    
    SURFRADobs(i)%stnlat(5) = 34.25
    SURFRADobs(i)%stnlon(5) = -89.87
    
    SURFRADobs(i)%stnlat(6) = 40.72
    SURFRADobs(i)%stnlon(6) = -77.93
    
    SURFRADobs(i)%stnlat(7) = 43.73
    SURFRADobs(i)%stnlon(7) = -96.62

    SURFRADobs(i)%startFlag = .true. 

  end subroutine SURFRAD_obsinit


end module SURFRAD_obsMod
