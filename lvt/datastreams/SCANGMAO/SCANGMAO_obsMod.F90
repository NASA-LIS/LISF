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
! !MODULE: SCANGMAO_obsMod
! \label(SCANGMAO_obsMod)
!
! !INTERFACE:
module SCANGMAO_obsMod
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
!  Climate Analysis Network (SCAN) data, quality controlled and 
!  reprocessed at NASA GMAO. The network is operated
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
  PUBLIC :: SCANGMAO_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SCANGMAOobs
!EOP
  type, public :: scangmaoobsdec
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

  end type scangmaoobsdec

  type(scangmaoobsdec), allocatable:: scangmaoobs(:)

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
! !ROUTINE: SCANGMAO_obsInit
! \label{SCANGMAO_obsInit}
! 
! !INTERFACE: 
  subroutine SCANGMAO_obsinit(i)
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
!  for reading SCANGMAO data. 
! 
!EOP
    integer                 :: status, rc
    integer                 :: ftn, k
    real*8                  :: tdur
    integer                 :: syr, smo, sda, shr, smn, sss
    integer                 :: eyr, emo, eda, ehr, emn, ess
    integer                 :: ts
    character*2             :: cdum
    character*100           :: coordfile

    if(.not.allocated(scangmaoobs)) then 
       allocate(scangmaoobs(LVT_rc%nDataStreams))
    endif

    scangmaoobs(i)%startFlag = .true. 

    call ESMF_ConfigGetAttribute(LVT_Config, scangmaoobs(i)%odir, &
         label='SCAN (GMAO) observation directory:', rc=status)
    call LVT_verify(status, 'SCANGMAO observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, scangmaoobs(i)%nstns, &
         label='SCAN (GMAO) number of stations:', rc=status)
    call LVT_verify(status, 'SCANGMAO number of stations: not defined')
    
    call ESMF_ConfigGetAttribute(LVT_Config, coordfile, &
         label='SCAN (GMAO) coord file:',rc=status)
    call LVT_verify(status, 'SCANGMAO coord file: not defined')

    call ESMF_TimeIntervalSet(scangmaoobs(i)%timestep, s=ts, rc=status)
    call LVT_verify(status, 'error in setting timestep (scangmaoobs)')

    call ESMF_TimeSet(scangmaoobs(i)%startTime,  yy=LVT_rc%yr, &
         mm = 1, &
         dd = 1, &
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status, 'error in setting scan start time')

!    call ESMF_TimeSet(scangmaoobs(i)%stopTime, yy=eyr, &
!         mm = emo, &
!         dd = eda, &
!         h = ehr, &
!         m = emn, &
!         calendar = LVT_calendar, &
!         rc=status)
!    call LVT_verify(status, 'error in setting scan stop time')

    call ESMF_TimeIntervalSet(scangmaoobs(i)%timestep, s=3600, rc=status)
    call LVT_verify(status, 'error in setting timestep (scangmaoobs)')
        
    call LVT_update_timestep(LVT_rc, 3600)
!yearly data = 24*366
    scangmaoobs(i)%nts = 8784
!    scangmaoobs(i)%nstns = 181

    allocate(scangmaoobs(i)%stnid(scangmaoobs(i)%nstns))
    allocate(scangmaoobs(i)%stnlat(scangmaoobs(i)%nstns))
    allocate(scangmaoobs(i)%stnlon(scangmaoobs(i)%nstns))

    ftn = LVT_getNextUnitNumber()
    open(ftn,file=trim(coordfile),form='formatted')
    
    write(LVT_logunit,*) '[INFO] opening coord file ',trim(coordfile)
    do k=1,scangmaoobs(i)%nstns
       read(ftn,*) cdum, scangmaoobs(i)%stnid(k), &
            scangmaoobs(i)%stnlat(k), scangmaoobs(i)%stnlon(k)
       write (LVT_logunit,*) '[INFO] ',scangmaoobs(i)%stnid(k), &
            scangmaoobs(i)%stnlat(k), scangmaoobs(i)%stnlon(k)
    enddo

    call LVT_releaseUnitNumber(ftn)

!-------------------------------------------------------------------------
!  SCANGMAO data contains soil moisture, soil temperature (5 layers) 
!-------------------------------------------------------------------------
    allocate(scangmaoobs(i)%soilt(scangmaoobs(i)%nstns, scangmaoobs(i)%nts))
    allocate(scangmaoobs(i)%sm(scangmaoobs(i)%nstns, scangmaoobs(i)%nts))

    scangmaoobs(i)%soilt = LVT_rc%udef
    scangmaoobs(i)%sm = LVT_rc%udef

    allocate(scangmaoobs(i)%roott(scangmaoobs(i)%nstns, scangmaoobs(i)%nts))
    allocate(scangmaoobs(i)%rootsm(scangmaoobs(i)%nstns, scangmaoobs(i)%nts))

    scangmaoobs(i)%roott = LVT_rc%udef
    scangmaoobs(i)%rootsm = LVT_rc%udef


    scangmaoobs(i)%yr = -1
  end subroutine SCANGMAO_obsinit


end module SCANGMAO_obsMod
