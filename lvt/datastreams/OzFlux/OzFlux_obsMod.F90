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
! !MODULE: OzFlux_obsMod
!  \label(OzFlux_obsMod)
!
! !INTERFACE:
module OzFlux_obsMod
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
!  This module handles the observation plugin for the OzFlux
!  in-situ data. 
!  
!  web link: http://www.ozflux.org.au
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 Apr 2020;   Sujay Kumar  Initial Specification
! 
!EOP
! 
! 
!
  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: OzFlux_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: OzFluxobs
!EOP
  type, public :: ozfluxobsdec
     character*100           :: odir
     integer                 :: n_stns 

     real,           allocatable :: stnlat(:)
     real,           allocatable :: stnlon(:)
     integer,        allocatable :: stncol(:)
     integer,        allocatable :: stnrow(:)
     character*100,  allocatable :: stn_name(:)
     integer, allocatable        :: tindex(:,:)

     logical                 :: startFlag
     integer                 :: nts
     real                    :: udef
     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: ts

     integer                 :: yr

     real, allocatable           :: qle(:,:)
     real, allocatable           :: qh(:,:)
     real, allocatable           :: qg(:,:)

  end type ozfluxobsdec

  type(ozfluxobsdec), allocatable:: ozfluxobs(:)

contains
  
!BOP
! 
! !ROUTINE: OzFlux_obsInit
! \label{OzFlux_obsInit}
!
! !INTERFACE: 
  subroutine OzFlux_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod
    use map_utils

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading OzFlux data. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!! !ARGUMENTS: 
    integer,   intent(IN) :: i 
!!EOP
    integer                 :: status
    integer                 :: ftn
    character*100           :: stnlist_file
    character*100           :: currentLine
    integer                 :: k, iloc
    real                    :: col,row

    if(.not.allocated(ozfluxobs)) then 
       allocate(ozfluxobs(LVT_rc%nDataStreams))
    endif
!
    ozfluxobs(i)%startFlag = .true. 

    call ESMF_ConfigGetAttribute(LVT_Config, ozfluxobs(i)%odir, &
         label='OzFlux observation directory:', rc=status)
    call LVT_verify(status, 'OzFlux observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, stnlist_file, &
         label='OzFlux station list file:', rc=status)
    call LVT_verify(status, 'OzFlux station list file: not defined')

    write(LVT_logunit,*) '[INFO] Processing OzFlux data locations '

    ftn = LVT_getNextUnitNumber()
    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) ozfluxobs(i)%n_stns
    read(ftn,*) 

    allocate(ozfluxobs(i)%stn_name(ozfluxobs(i)%n_stns))
    allocate(ozfluxobs(i)%stnlat(ozfluxobs(i)%n_stns))
    allocate(ozfluxobs(i)%stnlon(ozfluxobs(i)%n_stns))

    allocate(ozfluxobs(i)%stncol(ozfluxobs(i)%n_stns))
    allocate(ozfluxobs(i)%stnrow(ozfluxobs(i)%n_stns))

    do k=1,ozfluxobs(i)%n_stns
       read(ftn,'(a)') currentLine
       iloc = Index(currentLine, ";")
       read(currentLine(1:iloc -1), *) ozfluxobs(i)%stn_name(k)
       currentLine = currentLine(iloc+1:Len(currentLine))
       
       iloc = Index(currentLine, ";")
       read(currentLine(1:iloc -1), *) ozfluxobs(i)%stnlat(k)
       currentLine = currentLine(iloc+1:Len(currentLine))

       read(currentLine, *) ozfluxobs(i)%stnlon(k)
       
       write(LVT_logunit,*) '[INFO] ',ozfluxobs(i)%stn_name(k), &
            ozfluxobs(i)%stnlat(k), ozfluxobs(i)%stnlon(k)

       call latlon_to_ij(LVT_domain%lvtproj, ozfluxobs(i)%stnlat(k), &
            ozfluxobs(i)%stnlon(k),&
            col, row)
       
       ozfluxobs(i)%stncol(k) = nint(col)
       ozfluxobs(i)%stnrow(k) = nint(row)
    enddo
    
    call LVT_releaseUnitNumber(ftn)

    ozfluxobs(i)%udef   = 9999.9


    call ESMF_TimeIntervalSet(ozfluxobs(i)%ts, s=1800,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: OzFlux_obsMod ')

    allocate(ozfluxobs(i)%qle(ozfluxobs(i)%n_stns, 17569))
    allocate(ozfluxobs(i)%qh(ozfluxobs(i)%n_stns, 17569))
    allocate(ozfluxobs(i)%qg(ozfluxobs(i)%n_stns, 17569))
    allocate(ozfluxobs(i)%tindex(ozfluxobs(i)%n_stns, 17569))        
    ozfluxobs(i)%qle = LVT_rc%udef
    ozfluxobs(i)%qh  = LVT_rc%udef
    ozfluxobs(i)%qg  = LVT_rc%udef

    ozfluxobs(i)%yr = -1

    call LVT_update_timestep(LVT_rc, 1800)
  end subroutine OzFlux_obsinit

end module OzFlux_obsMod
