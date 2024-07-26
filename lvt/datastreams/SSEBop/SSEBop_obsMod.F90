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
! !MODULE: SSEBop_obsMod
! \label(SSEBop_obsMod)
!
! !INTERFACE:
module SSEBop_obsMod
! 
! !USES:   
  use ESMF
  use map_utils

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
!  18 May 2011   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SSEBop_obsinit !Initializes structures for reading SSEBop data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SSEBopobs !Object to hold SSEBop observation attributes
!EOP

  type, public :: ssebopdec
     character*100           :: odir
     integer                 :: use_anomaly
     integer                 :: nc, nr
     integer, allocatable    :: n11(:)
     real,    allocatable    :: et_var(:)
     integer                 :: yr
     integer                 :: mo
     real                    :: gridDesc(50)
     logical                 :: startFlag
     type(proj_info)         :: map_proj
  end type ssebopdec
     
  type(ssebopdec), allocatable :: SSEBopObs(:)

contains
  
!BOP
! 
! !ROUTINE: SSEBop_obsInit
! \label{SSEBop_obsInit}
!
! !INTERFACE: 
  subroutine SSEBop_obsinit(i)
! 
! !USES: 
    use LVT_coreMod,   only : LVT_rc, LVT_Config
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the SSEBop data, including the computation of spatial 
!   interpolation weights. The SSEBop data is provides in the 
!   EASE grid projection. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status
    real                  :: cornerlat1, cornerlat2
    real                  :: cornerlon1, cornerlon2

    if(.not.allocated(SSEBopObs)) then 
       allocate(SSEBopObs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SSEBopObs(i)%odir, &
         label='SSEBop data directory: ',rc=status)
    call LVT_verify(status, 'SSEBop data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SSEBopObs(i)%use_anomaly, &
         label='SSEBop process anomaly data: ',rc=status)
    call LVT_verify(status, 'SSEBop process anomaly data: not defined')

    call LVT_update_timestep(LVT_rc, 2592000)

    allocate(SSEBopobs(i)%et_var(LVT_rc%lnc*LVT_rc%lnr))

    ssebopobs(i)%gridDesc = 0
    
    SSEBopobs(i)%mo = LVT_rc%mo
    SSEBopobs(i)%yr = -1
    ssebopobs(i)%startFlag = .true.

    if(LVT_rc%tavgInterval.lt.2592000) then 
       write(LVT_logunit,*) '[ERR] The time averaging interval must be greater than'
       write(LVT_logunit,*) '[ERR] equal to a month since the SSEBop data is monthly'
       call LVT_endrun()
    endif

  end subroutine SSEBop_obsinit


end module SSEBop_obsMod
