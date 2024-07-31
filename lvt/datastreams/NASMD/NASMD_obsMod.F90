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
! !MODULE: NASMD_obsMod
! \label(NASMD_obsMod)
!
! !INTERFACE:
module NASMD_obsMod
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
!  This module handles the observation plugin for the North American
!  Soil Moisture Database (NASMD).
! 
!   http://soilmoisturemaps.tamu.edu/
! 
!   The North American Soil Moisture Database (NASMD) is a harmonized and 
!   quality-controlled soil moisture dataset that helps investigate 
!   land-atmosphere interactions, validates the accuracy of soil moisture 
!   simulations in global land-surface models and from satellite platforms, 
!   and describes how soil moisture influences climate on seasonal to 
!   inter-annual timescales.
!
!   Download the data through the interactive map:
!    http://soilmoisturemaps.tamu.edu/Data/Map/
!
!   The NASMD was developed and constructed at the Department of Geography's
!   Climate Science Lab at Texas A&M University.
!
!  References:
!
!   Ford, T. W. and S. M. Quiring (submitted) Comparison and application of
!   multiple methods for the interpolation of soil moisture observations. 
!   International Journal of Climatology.
!
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  07 Mar 2014  Sujay Kumar  Initial Specification
! 
!EOP
! !USES: 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: NASMD_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: NASMDobs
!EOP
  type, public :: nasmdobsdec
     character*100             :: odir
     integer                   :: nstns
     character*10,     allocatable :: stnid(:)
     real,             allocatable :: stnlat(:)
     real,             allocatable :: stnlon(:)
     integer,          allocatable :: stnnlayers(:)
     real,             allocatable :: stndepths(:,:)
     logical                 :: startFlag
     integer                 :: nts
     real                    :: udef
     type(ESMF_Time)         :: startTime, stopTime
     type(ESMF_TimeInterval) :: timeStep
     
     integer                 :: yr 

     real, allocatable           :: sfsm(:,:) 
     real, allocatable           :: rzsm(:,:) 
  end type nasmdobsdec

  type(nasmdobsdec), allocatable:: nasmdobs(:)

contains
  
!BOP
! 
! !ROUTINE:
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
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ROUTINE: NASMD_obsInit
! \label{NASMD_obsInit}
! 
! !INTERFACE: 
  subroutine NASMD_obsinit(i)
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
! !ARGUMENTS: 
    integer,   intent(IN) :: i 
! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading NASMD data. 
! 
!EOP
    integer                 :: status, rc
    integer                 :: ftn, k,j
    real*8                  :: tdur
    integer                 :: syr, smo, sda, shr, smn, sss
    integer                 :: eyr, emo, eda, ehr, emn, ess
    integer                 :: ts
    character*2             :: cdum
    character*100           :: coordfile

    if(.not.allocated(nasmdobs)) then 
       allocate(nasmdobs(LVT_rc%nDataStreams))
    endif

    nasmdobs(i)%startFlag = .true. 

    call ESMF_ConfigGetAttribute(LVT_Config, nasmdobs(i)%odir, &
         label='NASMD observation directory:', rc=status)
    call LVT_verify(status, 'NASMD observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, nasmdobs(i)%nstns, &
         label='NASMD number of stations:', default = 1289,rc=status)
    call LVT_verify(status, 'NASMD number of stations: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_Config, coordfile, &
         label='NASMD coord file:',rc=status)
    call LVT_verify(status, 'NASMD coord file: not defined')

    call ESMF_TimeIntervalSet(nasmdobs(i)%timestep, s=ts, rc=status)
    call LVT_verify(status, 'error in setting timestep (nasmdobs(i))')

    call ESMF_TimeSet(nasmdobs(i)%startTime,  yy=LVT_rc%yr, &
         mm = 1, &
         dd = 1, &
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status, 'error in setting nasmd start time')

!    call ESMF_TimeSet(nasmdobs(i)%stopTime, yy=eyr, &
!         mm = emo, &
!         dd = eda, &
!         h = ehr, &
!         m = emn, &
!         calendar = LVT_calendar, &
!         rc=status)
!    call LVT_verify(status, 'error in setting nasmd stop time')

    call ESMF_TimeIntervalSet(nasmdobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (nasmdobs(i))')
        
    call LVT_update_timestep(LVT_rc, 86400)
!yearly data = 24*366
    nasmdobs(i)%nts = 366
!    nasmdobs(i)%nstns = 181

    allocate(nasmdobs(i)%stnid(nasmdobs(i)%nstns))
    allocate(nasmdobs(i)%stnlat(nasmdobs(i)%nstns))
    allocate(nasmdobs(i)%stnlon(nasmdobs(i)%nstns))
    allocate(nasmdobs(i)%stnnlayers(nasmdobs(i)%nstns))
    allocate(nasmdobs(i)%stndepths(nasmdobs(i)%nstns,9))

    nasmdobs(i)%stndepths = 0 

    ftn = LVT_getNextUnitNumber()
    open(ftn,file=trim(coordfile),form='formatted')
    write(LVT_logunit,*) '[INFO] opening coord file ',trim(coordfile)
    read(ftn,*) 

    do k=1,nasmdobs(i)%nstns
!       read(ftn,fmt='(A10,2F10.3,I2)') nasmdobs(i)%stnid(k), &
       read(ftn,fmt=*) nasmdobs(i)%stnid(k), &
            nasmdobs(i)%stnlat(k), nasmdobs(i)%stnlon(k),&
            nasmdobs(i)%stnnlayers(k),&
            (nasmdobs(i)%stndepths(k,j),j=1,nasmdobs(i)%stnnlayers(k))
       write (LVT_logunit,*)  '[INFO] ',nasmdobs(i)%stnid(k), &
            nasmdobs(i)%stnlat(k), nasmdobs(i)%stnlon(k),&
            nasmdobs(i)%stnnlayers(k),&
            (nasmdobs(i)%stndepths(k,j),j=1,nasmdobs(i)%stnnlayers(k))
    enddo
    call LVT_releaseUnitNumber(ftn)

!-------------------------------------------------------------------------
!  NASMD data contains soil moisture, soil temperature (5 layers) 
!-------------------------------------------------------------------------

    allocate(nasmdobs(i)%sfsm(nasmdobs(i)%nstns, nasmdobs(i)%nts))

    nasmdobs(i)%sfsm = LVT_rc%udef

    allocate(nasmdobs(i)%rzsm(nasmdobs(i)%nstns, nasmdobs(i)%nts))

    nasmdobs(i)%rzsm = LVT_rc%udef

    nasmdobs(i)%yr = -1
  end subroutine NASMD_obsinit


end module NASMD_obsMod
