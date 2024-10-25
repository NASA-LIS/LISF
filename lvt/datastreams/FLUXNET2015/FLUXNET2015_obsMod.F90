!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! 
! !MODULE: FLUXNET2015_obsMod
!  \label(FLUXNET2015_obsMod)
!
! !INTERFACE:
! 
! !USES:   
!
! !DESCRIPTION: 
!  
!   This module provides the plugin for handling the processing of 
!   FLUXNET2015 dataset. The FLUXNET2015 Dataset includes data 
!   collected at sites from multiple regional flux networks. 
!   The previous versions of FLUXNET Dataset releases are the 
!   FLUXNET Marconi Dataset (2000) and the FLUXNET LaThuile Dataset (2007). 
!   The FLUXNET2015 Dataset includes several improvements to 
!   the data quality control protocols and the data processing pipeline.
!  
!   Website link: http://fluxnet.fluxdata.org/data/fluxnet2015-dataset/
!  
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 Feb 2017;   Sujay Kumar  Initial Specification
! 
!EOP
! 
! 
!
module FLUXNET2015_obsMod

  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: FLUXNET2015_obsinit
  PUBLIC :: FLUXNET2015obs

  type, public :: FLUXNET2015obsdec
     character*500 :: odir
     integer                 :: n_stns
     character*100, allocatable  :: stn_name(:)
     real,          allocatable  :: stnlat(:)
     real,          allocatable  :: stnlon(:)
     real                    :: change
     logical                 :: startflag
     integer                 :: yr
     integer                 :: nsmlayers
     integer                 :: nstlayers

     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: timestep

     real,          allocatable  :: Qle(:,:)        ! Latent heat flux
     real,          allocatable  :: Qh(:,:)         ! Sensible heat flux

  end type FLUXNET2015obsdec

  type(FLUXNET2015obsdec), allocatable :: FLUXNET2015obs(:)

contains
  
!BOP
! 
! !ROUTINE: FLUXNET2015_obsInit
! \label{FLUXNET2015_obsInit}
!
! !INTERFACE: 
 subroutine FLUXNET2015_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine should contain code to initialize and set up the data 
!  structures required for reading FLUXNET2015 data.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! 
! !ARGUMENTS: 
    integer,   intent(IN) :: i 
! 
!
!EOP
    character*100 :: stnlist_file
    integer       :: ftn 
    integer       :: k,iloc
    integer       :: arrayLen
    character*200 :: currentLine
    real          :: swcd(2), tsd(2)
    integer       :: status

    if(.not.allocated(FLUXNET2015obs)) then 
       allocate(FLUXNET2015obs(LVT_rc%nDataStreams))
    endif

!------------------------------------------------------------------------------
! Read any runtime specifications from the lvt.config file. 
!------------------------------------------------------------------------------

    call ESMF_ConfigGetAttribute(LVT_config, FLUXNET2015obs(i)%odir, &
         label='FLUXNET2015 observation directory:',rc=status)
    call LVT_verify(status, 'FLUXNET2015 observation directory: not defined')
  
    call ESMF_ConfigGetAttribute(LVT_config, stnlist_file, &
         label='FLUXNET2015 station list file:',rc=status)
    call LVT_verify(status, 'FLUXNET2015 station list file: not defined')

    ftn = LVT_getNextUnitNumber()
    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) FLUXNET2015obs(i)%n_stns
    read(ftn,*) 

    allocate(FLUXNET2015obs(i)%stn_name(FLUXNET2015obs(i)%n_stns))
    allocate(FLUXNET2015obs(i)%stnlat(FLUXNET2015obs(i)%n_stns))
    allocate(FLUXNET2015obs(i)%stnlon(FLUXNET2015obs(i)%n_stns))

!---------------------------------------------------------------------------
! Leap year test, following rules by Microsoft. The depth of this test may not
! be entirely necessary due to how recent the data is, but it is correct.
!---------------------------------------------------------------------------   
    if (Mod(LVT_rc%yr, 4) == 0) Then            !Step1
        if (Mod(LVT_rc%yr, 100) == 0) Then      !Step2
           if (Mod(LVT_rc%yr, 400) == 0) Then   !Step3
              arrayLen = 17568                  !Step4, leap LVT_rc%yr
           else                           
              arrayLen = 17520                  !Step5, not leap LVT_rc%yr
           end if
        else
           arrayLen = 17568                     !Step4, leap LVT_rc%yr
        end if
     else
        arrayLen = 17520                        !Step5, not leap LVT_rc%yr
     End if
!svk : to be safe:
     arrayLen = 17570

    allocate(FLUXNET2015obs(i)%Qle(FLUXNET2015obs(i)%n_stns, arrayLen))
    allocate(FLUXNET2015obs(i)%Qh(FLUXNET2015obs(i)%n_stns, arrayLen))

    FLUXNET2015obs(i)%Qle = LVT_rc%udef
    FLUXNET2015obs(i)%Qh = LVT_rc%udef

    write(LVT_logunit,*) '[INFO] Processing FLUXNET2015 stations'
!------------------------------------------------------------------------------
! For each station, this reads the site name, station name, and station
! position in lat, lon coordinates.
!------------------------------------------------------------------------------
    do k=1,FLUXNET2015obs(i)%n_stns
       read(ftn,'(a)') currentLine

       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) FLUXNET2015obs(i)%stn_name(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) FLUXNET2015obs(i)%stnlat(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) FLUXNET2015obs(i)%stnlon(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       write(LVT_logunit,*) '[INFO] ',k, &
            FLUXNET2015obs(i)%stn_name(k), &
            FLUXNET2015obs(i)%stnlat(k), FLUXNET2015obs(i)%stnlon(k)
    end do
    
    call LVT_releaseUnitNumber(ftn)
!------------------------------------------------------------------------------
! Sets the start flag to true since this is the start of data. The flag is then
! set to false after the initial file is read so that a switch in year will not
! be mistaken for a new data set.
!------------------------------------------------------------------------------

    FLUXNET2015obs(i)%startflag = .true. 
    FLUXNET2015obs(i)%yr = -1

    call ESMF_TimeIntervalSet(FLUXNET2015obs(i)%timestep, s=1800, rc=status)
    call LVT_verify(status, 'error in setting timestep (FLUXNET2015obs)')

    call LVT_update_timestep(LVT_rc, 1800)

  end subroutine FLUXNET2015_obsinit

end module FLUXNET2015_obsMod
