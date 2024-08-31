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
! !MODULE: JULES_obsMod
! \label(JULES_obsMod)
!
! !INTERFACE:
module JULES_obsMod
! 
! !USES: 
  use ESMF

  implicit none

  PRIVATE 
!
! !DESCRIPTION: 
!   This module handles the observation plugin for the JULES model
!   output data. The code assumes that the output is provided in a
!   single file in NetCDF format. The time dimension is assumed
!   to be part of the file structure. 
!   
! !FILES USED:
!
! !REVISION HISTORY: 
!  8 July 2015   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: JULES_obsinit !Initializes structures for reading JULES data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: JULESobs !Object to hold JULES observation attributes
!EOP

  type, public :: julesdec
     character*500               :: odir
     type(ESMF_Time)             :: refTime
     integer                     :: ntimes
     logical                     :: startMode

     real,  allocatable      :: rainf_jules(:,:,:)
     real,  allocatable      :: snowf_jules(:,:,:)
     real,  allocatable      :: qle_jules(:,:,:)
     real,  allocatable      :: qh_jules(:,:,:)
     real,  allocatable      :: qtau_jules(:,:,:)
     real,  allocatable      :: evap_jules(:,:,:)
     real,  allocatable      :: smc_jules(:,:,:,:)
     real,  allocatable      :: stc_jules(:,:,:,:)
     real,  allocatable      :: pstar_jules(:,:,:)
     real,  allocatable      :: tstar_jules(:,:,:)
     real,  allocatable      :: gpp_jules(:,:,:)
     real,  allocatable      :: time_val(:)
     real,  allocatable      :: lat(:,:), lon(:,:)
     
     integer                 :: nx, ny, nsoil


  end type julesdec
     
  type(julesdec), allocatable :: JULESObs(:)

contains
  
!BOP
! 
! !ROUTINE: JULES_obsInit
! \label{JULES_obsInit}
!
! !INTERFACE: 
  subroutine JULES_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the JULES data. The user is required to provide the 
!   following details through the LVT config file: 
!    * Name of the JULES data file
!    * Data reference time for the JULES data. It should be specified as:
!       <yr> <mo> <da> <hr> <mn> <ss> where yr is a 4 digit year, 
!       mo is a 2 digit month, da is a 2 digit day, hr is a 2 digit hour
!       mn is a 2 digit minute and ss is a 2 digit second. 
!            
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: i
    integer               :: k
    integer               :: time(6)
    integer               :: status
    integer               :: tint

    if(.not.allocated(JULESobs)) then 
       allocate(JULESobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, JULESobs(i)%odir, &
         label='JULES data file: ',rc=status)
    call LVT_verify(status, 'JULES data file: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, tint, &
         label='JULES timestep: ',rc=status)
    call LVT_verify(status, 'JULES timestep: not defined')

    call ESMF_ConfigFindLabel(LVT_config,"JULES data reference time:",rc=status)
    do k=1,6
       call ESMF_ConfigGetAttribute(LVT_config,time(k),rc=status)
       if(status.ne.0) then 
          write(LVT_logunit,*) '[ERR] please input the reference time as '
          write(LVT_logunit,*) '[ERR] JULES data reference time: <yr> <mo> <da> <hr> <mn> <ss>'
          call LVT_endrun()
       endif
    enddo

     call ESMF_TimeSet(JULESobs(i)%refTime,  yy=time(1), &
          mm = time(2), &
          dd = time(3),&
          h = time(4),&
          m = time(5),&
          s = time(6),&
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status, 'ESMF_TimeSet error in JULES_obsInit')

    call LVT_update_timestep(LVT_rc, 1800)
    
    JULESobs(i)%startMode  = .true. 
  end subroutine JULES_obsinit

end module JULES_obsMod
