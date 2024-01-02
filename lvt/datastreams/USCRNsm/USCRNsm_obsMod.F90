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
! !MODULE: USCRNsm_obsMod
! \label(USCRNsm_obsMod)
!
! !INTERFACE:
module USCRNsm_obsMod
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
!  This module handles the observation plugin for soil moisture measurements
!  from the U.S. Climate Reference Network (USCRN), NOAA's premiere 
!  surface reference network. 
!
!  https://www.ncdc.noaa.gov/data-access/land-based-station-data/land-based-datasets/us-climate-reference-network-uscrn
! 
! !REVISION HISTORY: 
!  18 Jan 2017   Sujay Kumar  Initial Specification
! 
!EOP
! !USES: 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: USCRNsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: USCRNsmobs
!EOP
  type, public :: uscrnsmobsdec
     character*100        :: odir
     integer              :: nstns
     character*40, allocatable :: stnid(:)
     real,         allocatable :: stnlat(:)
     real,         allocatable :: stnlon(:)
     logical                 :: startFlag
     integer                 :: nts
     real                    :: udef
     type(ESMF_Time)         :: startTime, stopTime
     type(ESMF_TimeInterval) :: timeStep
     
     integer                 :: yr 

!     real                    :: st_wts(5)
!     real                    :: sm_wts(5)

     real, allocatable       :: soilt(:,:)  
     real, allocatable       :: sm(:,:) 
     real, allocatable       :: rootsm(:,:)  
     real, allocatable       :: roott(:,:)  

  end type uscrnsmobsdec

  type(uscrnsmobsdec), allocatable :: uscrnsmobs(:)

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
! !ROUTINE: USCRNsm_obsInit
! \label{USCRNsm_obsInit}
! 
! !INTERFACE: 
  subroutine USCRNsm_obsinit(i)
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
!  for reading USCRNsm data. 
! 
!EOP
    integer                 :: status, rc
    integer                 :: ftn, k
    real*8                  :: tdur
    integer                 :: syr, smo, sda, shr, smn, sss
    integer                 :: eyr, emo, eda, ehr, emn, ess
    integer                 :: ts, iloc
    character*100           :: coordfile, cline
    integer                 :: filecheck

    if(.not.allocated(uscrnsmobs)) then 
       allocate(uscrnsmobs(LVT_rc%nDataStreams))
    endif

    uscrnsmobs(i)%startFlag = .true. 

    call ESMF_ConfigGetAttribute(LVT_Config, uscrnsmobs(i)%odir, &
         label='USCRN soil moisture observation directory:', rc=status)
    call LVT_verify(status, 'USCRN soil moisture observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, coordfile, &
         label='USCRN soil moisture station file:', rc=status)
    call LVT_verify(status, 'USCRN soil moisture station file: not defined')

    ! This implementation uses daily files 
    call ESMF_TimeIntervalSet(uscrnsmobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (uscrnsmobs)')

    call ESMF_TimeSet(uscrnsmobs(i)%startTime,  yy=LVT_rc%yr, &
         mm = 1, &
         dd = 1, &
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status, 'error in setting uscrnsm start time')
        
    call LVT_update_timestep(LVT_rc, 86400)
    uscrnsmobs(i)%nts = 366

    write(LVT_logunit,*) '[INFO] Opening lat/lon coordinate file, ',trim(coordfile)

    ! Determine number of station points in lat/lon coordinate file:
    uscrnsmobs(i)%nstns = 0
    ftn = LVT_getNextUnitNumber()
    open(ftn,file=trim(coordfile),form='formatted')
    do; 
       read(ftn,*,iostat=filecheck) cLine
       if( filecheck < 0 ) exit
       uscrnsmobs(i)%nstns = uscrnsmobs(i)%nstns + 1
    end do;
    rewind(ftn)
    write(LVT_logunit,*) "[INFO] Number of USCRN Stations : ",uscrnsmobs(i)%nstns
!    uscrnsmobs(i)%nstns = 242  ! Total number of US CRN stations
!    uscrnsmobs(i)%nstns = 86   ! QC'd number based on 5-SM layers available

    allocate(uscrnsmobs(i)%stnid(uscrnsmobs(i)%nstns))
    allocate(uscrnsmobs(i)%stnlat(uscrnsmobs(i)%nstns))
    allocate(uscrnsmobs(i)%stnlon(uscrnsmobs(i)%nstns))

    ! Read in station coordinate information:
!    open(ftn,file=trim(coordfile),form='formatted')
    do k=1,uscrnsmobs(i)%nstns
       read(ftn,'(a)') cLine
       iloc = index(cline,",")
       read(cline(1:iloc-1),*) uscrnsmobs(i)%stnid(k)
       cline = cline(iloc+1:len(cline))

       iloc = index(cline,",")
       read(cline(1:iloc-1),*) uscrnsmobs(i)%stnlat(k)
       cline = cline(iloc+1:len(cline))

       iloc = index(cline,",")
       read(cline(1:iloc-1),*) uscrnsmobs(i)%stnlon(k)
       cline = cline(iloc+1:len(cline))

       write (LVT_logunit,*) '[INFO] ',uscrnsmobs(i)%stnid(k), &
            uscrnsmobs(i)%stnlat(k), uscrnsmobs(i)%stnlon(k)
    enddo

    call LVT_releaseUnitNumber(ftn)

!-------------------------------------------------------------------------
!  USCRNsm data contains soil moisture, soil temperature (5 layers) 
!-------------------------------------------------------------------------
    allocate(uscrnsmobs(i)%soilt(uscrnsmobs(i)%nstns, uscrnsmobs(i)%nts))
    allocate(uscrnsmobs(i)%sm(uscrnsmobs(i)%nstns, uscrnsmobs(i)%nts))

    uscrnsmobs(i)%soilt = LVT_rc%udef
    uscrnsmobs(i)%sm = LVT_rc%udef

    allocate(uscrnsmobs(i)%roott(uscrnsmobs(i)%nstns, uscrnsmobs(i)%nts))
    allocate(uscrnsmobs(i)%rootsm(uscrnsmobs(i)%nstns, uscrnsmobs(i)%nts))

    uscrnsmobs(i)%roott = LVT_rc%udef
    uscrnsmobs(i)%rootsm = LVT_rc%udef


    uscrnsmobs(i)%yr = -1

  end subroutine USCRNsm_obsinit


end module USCRNsm_obsMod
