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
! !MODULE: PBOH2O_obsMod
!  \label(PBOH2O_obsMod)
!
! !INTERFACE:
module PBOH2O_obsMod
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
!  This module handles the observation plugin for the plate boundary
!  observatory (PBO) in-situ data. The datasets handled
!  in this plugin include measurements of soil moisture and snow. 
!  
!  web link: http://xenon.colorado.edu/portal/index.php
!
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 May 2014   Sujay Kumar  Initial Specification
! 
!EOP
! 
! 
!
  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: PBOH2O_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: PBOH2Oobs
!EOP
  type, public :: pboh2oobsdec
     character*100           :: odir
     integer                 :: n_stns 
     integer                 :: yr
     character*30            :: site_id
     type(ESMF_Time)         :: startTime
     type(ESMF_TimeInterval) :: timestep
     real,           allocatable :: stnlat(:)
     real,           allocatable :: stnlon(:)
     real,           allocatable :: snod(:,:)
     real,           allocatable :: swe(:,:)
     real,           allocatable :: sm(:,:)

     character*100,  allocatable :: stn_name(:)

     logical                 :: startFlag

  end type pboh2oobsdec

  type(pboh2oobsdec), allocatable:: pboh2oobs(:)

contains
  
!BOP
! 
! !ROUTINE: PBOH2O_obsInit
! \label{PBOH2O_obsInit}
!
! !INTERFACE: 
  subroutine PBOH2O_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod


    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading PBOH2O data. 
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

    if(.not.allocated(pboh2oobs)) then 
       allocate(pboh2oobs(LVT_rc%nDataStreams))
    endif

    pboh2oobs(i)%startFlag = .true. 

    call ESMF_ConfigGetAttribute(LVT_Config, pboh2oobs(i)%odir, &
         label='PBOH2O observation directory:', rc=status)
    call LVT_verify(status, 'PBOH2O observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, stnlist_file, &
         label='PBOH2O station list file:', rc=status)
    call LVT_verify(status, 'PBOH2O station list file: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    ftn = LVT_getNextUnitNumber()
    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) pboh2oobs(i)%n_stns
    read(ftn,*) 

    allocate(pboh2oobs(i)%stn_name(pboh2oobs(i)%n_stns))
    allocate(pboh2oobs(i)%stnlat(pboh2oobs(i)%n_stns))
    allocate(pboh2oobs(i)%stnlon(pboh2oobs(i)%n_stns))

    allocate(pboh2oobs(i)%snod(pboh2oobs(i)%n_stns,366))
    allocate(pboh2oobs(i)%swe(pboh2oobs(i)%n_stns,366))
    allocate(pboh2oobs(i)%sm(pboh2oobs(i)%n_stns,366))

    do k=1,pboh2oobs(i)%n_stns
       read(ftn,'(a)') pboh2oobs(i)%stn_name(k)
       write(LVT_logunit,*) '[INFO] ',pboh2oobs(i)%stn_name(k)
!            pboh2oobs(i)%stnlat(k), pboh2oobs(i)%stnlon(k),pboh2oobs(i)%stnbd(k)
    enddo
    
    call LVT_releaseUnitNumber(ftn)

!-------------------------------------------------------------------------
!  PBOH2O data contains soil moisture, soil temperature (5 layers) 
!-------------------------------------------------------------------------

    call ESMF_TimeIntervalSet(pboh2oobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (pboh2oobs')
    pboh2oobs(i)%yr = -1

  end subroutine PBOH2O_obsinit


end module PBOH2O_obsMod
