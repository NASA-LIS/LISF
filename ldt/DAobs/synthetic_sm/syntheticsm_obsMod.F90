!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: syntheticsm_obsMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module syntheticsm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: syntheticsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: syntheticsmobs
!EOP
  type, public :: syntheticsmobsdec

     character(len=LDT_CONST_PATH_LEN) :: odir
     integer                    :: mo
     real,    allocatable       :: smobs(:,:)
     logical                    :: startmode 
  end type syntheticsmobsdec

  type(syntheticsmobsdec), allocatable:: syntheticsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: syntheticsm_obsInit
! \label{syntheticsm_obsInit}
! 
! !INTERFACE: 
  subroutine syntheticsm_obsinit()
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
!  for reading SYNTHETIC soil moisture data. 
! 
!EOP
    integer                 :: npts
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(syntheticsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'Synthetic soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, syntheticsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'Synthetic soil moisture observation directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       syntheticsmobs(n)%startmode = .true. 

       allocate(syntheticsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       syntheticsmobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
    enddo
  end subroutine syntheticsm_obsinit
     
end module syntheticsm_obsMod
