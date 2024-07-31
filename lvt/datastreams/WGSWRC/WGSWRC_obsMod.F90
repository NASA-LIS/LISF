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
! !MODULE: WGSWRC_obsMod
! \label(WGSWRC_obsMod)
!
! !INTERFACE:
module WGSWRC_obsMod
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


  PUBLIC :: WGSWRC_obsInit
  PUBLIC :: WGSWRCobs
  
  type, public :: wgswrcobsdec

     character*100         :: odir
     integer               :: yr
     integer               :: n_stns
     character*50, allocatable :: stn_name(:)
     type(ESMF_Time)       :: startTime
     type(ESMF_TimeInterval) :: timestep
     real,         allocatable :: sm(:,:,:)
     real,         allocatable :: rootsm(:,:,:)
     real,         allocatable :: st(:,:,:)
  end type wgswrcobsdec

  type(wgswrcobsdec), allocatable :: WGSWRCobs(:)

contains

!BOP
! 
! !ROUTINE: WGSWRC_obsInit
! \label(WGSWRC_obsInit)
!
! !INTERFACE:
  subroutine WGSWRC_obsInit(i)
! 
! !USES:   
    use LVT_coreMod
    use LVT_obsDataMod
    use LVT_logMod
    
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
    character*100           :: stnlist_file

    if(.not.allocated(WGSWRCobs)) then 
       allocate(WGSWRCobs(LVT_rc%nDataStreams))
    endif

    write(LVT_logunit,*) '[INFO] Initializing WGSWRC data reader....'
    call ESMF_ConfigGetAttribute(LVT_config, WGSWRCobs(i)%odir, &
         label='WG SWRC observation directory:',rc=status)
    call LVT_verify(status, 'WG SWRC observation directory: not defined')
  
    call ESMF_ConfigGetAttribute(LVT_config, stnlist_file, &
         label='WG SWRC station list file:',rc=status)
    call LVT_verify(status, 'WG SWRC station list file: not defined')

    ftn = LVT_getNextUnitNumber()

    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) WGSWRCobs(i)%n_stns
    read(ftn,*) 

    allocate(WGSWRCobs(i)%stn_name(WGSWRCobs(i)%n_stns))
    
    write(LVT_logunit,*) '[INFO] WGSWRC station list .. '
    do k=1,WGSWRCobs(i)%n_stns
       read(ftn,*) WGSWRCobs(i)%stn_name(k)
       write(LVT_logunit,*) '[INFO] ',trim(WGSWRCobs(i)%stn_name(k))
    enddo
    call LVT_releaseUnitNumber(ftn)

    WGSWRCobs(i)%yr = -1

    call ESMF_TimeIntervalSet(WGSWRCobs(i)%timestep,s=1800,rc=status)
    call LVT_verify(status,"ESMF_TimeIntervalSet failed in WGSWRC_obsInit")

    allocate(WGSWRCobs(i)%sm(LVT_rc%lnc,LVT_rc%lnr,366*48))
    allocate(WGSWRCobs(i)%rootsm(LVT_rc%lnc,LVT_rc%lnr,366*48))
    allocate(WGSWRCobs(i)%st(LVT_rc%lnc,LVT_rc%lnr,366*48))

    WGSWRCobs(i)%sm = -9999.0
    WGSWRCobs(i)%rootsm = -9999.0
    WGSWRCobs(i)%st = -9999.0

    write(LVT_logunit,*) '[INFO] Finished initializing WGSWRC data reader....'

  end subroutine WGSWRC_obsInit


end module WGSWRC_obsMod
