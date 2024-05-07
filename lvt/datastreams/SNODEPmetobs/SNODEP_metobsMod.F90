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
! !MODULE: SNODEP_metobsMod
! \label(SNODEP_metobsMod)
!
! !INTERFACE:
module SNODEP_metobsMod
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
!  This module handles the observation plugin for the MET obs used in 
!  generating the SNODEP product at the Air Force Weather Agency (AFWA). 
!  The plugin handles the snow depth measurements of SNODEP. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  13 Dec 2010   Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SNODEP_metobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SNODEPmetobs
!EOP
  type, public :: SNODEPmetobsdec

     character*100           :: odir
     integer                 :: da
     type(ESMF_Time)         :: startTime
     real, allocatable       :: snod(:,:,:)
     type(ESMF_TimeInterval) :: ts
     logical                 :: startflag
  end type SNODEPmetobsdec

  type(SNODEPmetobsdec), allocatable :: SNODEPmetobs(:)

contains
  
!BOP
! 
! !ROUTINE: SNODEP_metobsInit
! \label{SNODEP_metobsInit}
!
! !INTERFACE: 
  subroutine SNODEP_metobsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine should contain code to initialize and 
!  set up the data structures required for reading the specific data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer :: status

    if(.not.allocated(SNODEPmetobs)) then 
       allocate(SNODEPmetobs(LVT_rc%nDataStreams))
    endif
!---------------------------------------------------------------------------
! Read runtime specifications from the lvt.config file. 
!---------------------------------------------------------------------------

    call ESMF_ConfigGetAttribute(LVT_config, SNODEPmetobs(i)%odir, &
         label='SNODEP metobs directory:',rc=status)
    call LVT_verify(status, 'SNODEP metobs directory: not defined')
  
!---------------------------------------------------------------------------
! Initialize variables.   
!---------------------------------------------------------------------------
    allocate(SNODEPmetobs(i)%snod(LVT_rc%lnc,LVT_rc%lnr,24))

    SNODEPmetobs(i)%snod = LVT_rc%udef

    call ESMF_TimeIntervalSet(SNODEPmetobs(i)%ts, s=3600,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: SNODEP_metobs')
    
    call LVT_update_timestep(LVT_rc, 3600)

    SNODEPmetobs(i)%da = -1
    SNODEPmetobs(i)%startflag = .true. 

  end subroutine SNODEP_metobsinit

  
end module SNODEP_metobsMod
