!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: SMOSL2sm_obsMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module SMOSL2sm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSL2sm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSL2smobs
!EOP
  type, public :: smossmobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     integer                :: mo
     real,    allocatable       :: smobs(:,:)
     logical                :: startmode 
  end type smossmobsdec

  type(smossmobsdec), allocatable:: SMOSL2smobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMOSL2sm_obsInit
! \label{SMOSL2sm_obsInit}
! 
! !INTERFACE: 
  subroutine SMOSL2sm_obsinit()
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading SMOSL2 soil moisture data. 
! 
!EOP
    integer                 :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(SMOSL2smobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOS L2 soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOSL2smobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'SMOS L2 soil moisture observation directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       SMOSL2smobs(n)%startmode = .true. 

       allocate(SMOSL2smobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       SMOSL2smobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
    enddo
  end subroutine SMOSL2sm_obsinit
     
end module SMOSL2sm_obsMod
