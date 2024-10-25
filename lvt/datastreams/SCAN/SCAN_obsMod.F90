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
! !MODULE: SCAN_obsMod
! \label(SCAN_obsMod)
!
! !INTERFACE:
module SCAN_obsMod
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
!  This module handles the observation plugin for the NRCS's Soil 
!  Climate Analysis Network (SCAN) data. The network is operated
!  by the USDA NRCS National Water and Climate Center in Portland,
!  OR, with assistance from the USDA NRCS National Soil Survey 
!  Center in Lincon, NE. The system focuses on agricultural areas
!  of the U.S., monitoring soil temperature and soil moisture 
!  content at several depths. 
!  
! http://www.wcc.nrcs.usda.gov/scan/
! 
!  This implementation employs the quality controlled data processed
!  at the NASA Global Modeling and Assimilation Office (GMAO) by 
!  Rolf Reichle and Gabrielle de Lannoy. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Apr 2009   Sujay Kumar  Initial Specification
! 
!EOP
! !USES: 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SCAN_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SCANobs
!EOP
  type, public :: scanobsdec
     character*100        :: odir
     integer              :: nstns
     integer,     allocatable :: stnid(:)
     real,        allocatable :: stnlat(:)
     real,        allocatable :: stnlon(:)
     logical                 :: startFlag
     integer                 :: nts
     real                    :: udef
     type(ESMF_Time)         :: startTime, stopTime
     type(ESMF_TimeInterval) :: timeStep
     
     integer                 :: yr 

!     real                    :: st_wts(5)
!     real                    :: sm_wts(5)

     real, allocatable           :: soilt(:,:)  
     real, allocatable           :: sm(:,:) 
     real, allocatable           :: rootsm(:,:)  
     real, allocatable           :: roott(:,:)  

  end type scanobsdec

  type(scanobsdec), allocatable:: scanobs(:)

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
! !ROUTINE: SCAN_obsInit
! \label{SCAN_obsInit}
! 
! !INTERFACE: 
  subroutine SCAN_obsinit(i)
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
!  for reading SCAN data. 
! 
!EOP
    integer                 :: status, rc
    integer                 :: ftn, k
    real*8                  :: tdur
    integer                 :: syr, smo, sda, shr, smn, sss
    integer                 :: eyr, emo, eda, ehr, emn, ess
    integer                 :: ts
    character*100           :: coordfile

    if(.not.allocated(scanobs)) then 
       allocate(scanobs(LVT_rc%nDataStreams))
    endif

    scanobs(i)%startFlag = .true. 

    call ESMF_ConfigGetAttribute(LVT_Config, scanobs(i)%odir, &
         label='SCAN observation directory:', rc=status)
    call LVT_verify(status, 'SCAN observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, scanobs(i)%nstns, &
         label='SCAN number of stations:', rc=status)
    call LVT_verify(status, 'SCAN number of stations: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_Config, coordfile, &
         label='SCAN coord file:',rc=status)
    call LVT_verify(status, 'SCAN coord file: not defined')

    call ESMF_TimeIntervalSet(scanobs(i)%timestep, s=ts, rc=status)
    call LVT_verify(status, 'error in setting timestep (scanobs)')

    call ESMF_TimeSet(scanobs(i)%startTime,  yy=LVT_rc%yr, &
         mm = 1, &
         dd = 1, &
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status, 'error in setting scan start time')

!    call ESMF_TimeSet(scanobs(i)%stopTime, yy=eyr, &
!         mm = emo, &
!         dd = eda, &
!         h = ehr, &
!         m = emn, &
!         calendar = LVT_calendar, &
!         rc=status)
!    call LVT_verify(status, 'error in setting scan stop time')

    call ESMF_TimeIntervalSet(scanobs(i)%timestep, s=3600, rc=status)
    call LVT_verify(status, 'error in setting timestep (scanobs)')
        
    call LVT_update_timestep(LVT_rc, 3600)
!yearly data = 24*366
    scanobs(i)%nts = 8784
!    scanobs(i)%nstns = 181

    allocate(scanobs(i)%stnid(scanobs(i)%nstns))
    allocate(scanobs(i)%stnlat(scanobs(i)%nstns))
    allocate(scanobs(i)%stnlon(scanobs(i)%nstns))

    ftn = LVT_getNextUnitNumber()
    open(ftn,file=trim(coordfile),form='formatted')
    
    write(LVT_logunit,*) '[INFO] opening coord file ',trim(coordfile)
    do k=1,scanobs(i)%nstns
       read(ftn,*) scanobs(i)%stnid(k), &
            scanobs(i)%stnlat(k), scanobs(i)%stnlon(k)
       write (LVT_logunit,*) '[INFO] ',scanobs(i)%stnid(k), &
            scanobs(i)%stnlat(k), scanobs(i)%stnlon(k)
    enddo

    call LVT_releaseUnitNumber(ftn)

!-------------------------------------------------------------------------
!  SCAN data contains soil moisture, soil temperature (5 layers) 
!-------------------------------------------------------------------------
    allocate(scanobs(i)%soilt(scanobs(i)%nstns, scanobs(i)%nts))
    allocate(scanobs(i)%sm(scanobs(i)%nstns, scanobs(i)%nts))

    scanobs(i)%soilt = LVT_rc%udef
    scanobs(i)%sm = LVT_rc%udef

    allocate(scanobs(i)%roott(scanobs(i)%nstns, scanobs(i)%nts))
    allocate(scanobs(i)%rootsm(scanobs(i)%nstns, scanobs(i)%nts))

    scanobs(i)%roott = LVT_rc%udef
    scanobs(i)%rootsm = LVT_rc%udef


    scanobs(i)%yr = -1
  end subroutine SCAN_obsinit


end module SCAN_obsMod
