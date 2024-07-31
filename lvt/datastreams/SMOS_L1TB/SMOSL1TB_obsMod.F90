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
! !MODULE: SMOSL1TB_obsMod
! \label(SMOSL1TB_obsMod)
!
! !INTERFACE:
module SMOSL1TB_obsMod
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
!  This module handles the observation plugin for the Land Parameter
!  Retrieval Model (LPRM) AMSR-E soil moisture product
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  12 Dec 2014: Sujay Kumar, Initial Specification
! 
!EOP
! 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSL1TB_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSL1TBobs
!EOP
  type, public :: smosl2smobsdec

     character*100        :: odir
     logical             :: startmode     

  end type smosl2smobsdec

  type(smosl2smobsdec), allocatable:: SMOSL1TBobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMOSL1TB_obsInit
! \label{SMOSL1TB_obsInit}
!
! !INTERFACE: 
  subroutine SMOSL1TB_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading LPRM AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts, status
    real                  :: gridDesci(50)

    if(.not.allocated(SMOSL1TBobs)) then 
       allocate(SMOSL1TBobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMOSL1TBobs(i)%odir, &
         label='SMOS L1 TB observation directory:', rc=status)
    call LVT_verify(status, 'SMOS L1 TB observation directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)
        
    SMOSL1TBobs(i)%startmode = .true. 

  end subroutine SMOSL1TB_obsinit


end module SMOSL1TB_obsMod
