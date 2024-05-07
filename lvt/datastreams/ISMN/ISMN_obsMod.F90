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
! !MODULE: ISMN_obsMod
! \label(ISMN_obsMod)
!
! !INTERFACE:
module ISMN_obsMod
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
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  13 May 2011   Sujay Kumar  Initial Specification
! 
!EOP


  PUBLIC :: ISMN_obsInit
  PUBLIC :: ISMNobs
  
  type, public :: ismnstn
     integer                    :: vlevels
     real                       :: lat,lon
     character*500, allocatable :: fname(:)
     real,         allocatable  :: sm(:,:)
     real,         allocatable  :: sfsm(:)
     real,         allocatable  :: rzsm(:)

  end type ismnstn

  type, public :: ismnobsdec

     character*500               :: odir
     integer                     :: yr
     integer                     :: n_stns
     integer                     :: nts 
     type(ESMF_Time)             :: startTime
     type(ESMF_TimeInterval)     :: timestep
     type(ismnstn), allocatable  :: stn(:)

  end type ismnobsdec

  type(ismnobsdec), allocatable :: ISMNobs(:)

contains

!BOP
! 
! !ROUTINE: ISMN_obsInit
! \label(ISMN_obsInit)
!
! !INTERFACE:
  subroutine ISMN_obsInit(i)
! 
! !USES:   
    use LVT_coreMod,    only : LVT_rc, LVT_config
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod,     only : LVT_verify, LVT_logunit, &
         LVT_getNextUnitNumber, LVT_releaseUnitNumber
    
    implicit none
!
! !INPUT PARAMETERS: 
    integer,     intent(IN) :: i   ! index of the observation type
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer                 :: k 
    integer                 :: ftn
    integer                 :: status
    character*500           :: stnlist_file

    if(.not.allocated(ISMNobs)) then 
       allocate(ISMNobs(LVT_rc%nDataStreams))
    endif

    write(LVT_logunit,*) '[INFO] Initializing ISMN data reader....'
    call ESMF_ConfigGetAttribute(LVT_config, ISMNobs(i)%odir, &
         label='ISMN observation directory:',rc=status)
    call LVT_verify(status, 'ISMN observation directory: not defined')
      

    ISMNobs(i)%yr = -1

    call ESMF_TimeIntervalSet(ISMNobs(i)%timestep,s=3600,rc=status)
    call LVT_verify(status,"ESMF_TimeIntervalSet failed in ISMN_obsInit")

    ISMNobs(i)%nts = 8784 !24*366
    
    write(LVT_logunit,*) '[INFO] Finished initializing ISMN data reader....'

    allocate(ISMNobs(i)%stn(ISMNobs(i)%n_stns))

    call LVT_update_timestep(LVT_rc, 3600)

  end subroutine ISMN_obsInit


end module ISMN_obsMod
