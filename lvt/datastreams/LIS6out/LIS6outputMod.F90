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
! !MODULE: LIS6outputMod
! \label(LIS6outputMod)
!
! !INTERFACE:
module LIS6outputMod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the use of a LIS model simulation output as 
!  "observations". 
!  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LIS6outputInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lis6output
!
!EOP
  
  type, public :: lis6obsdec
     character*100  :: odir
  end type lis6obsdec

  type(lis6obsdec), allocatable :: lis6output(:)

contains

!BOP
! 
! !ROUTINE: LIS6outputInit
! \label{LIS6outputInit}
!
! !INTERFACE: 
  subroutine LIS6outputInit(i)
! 
! !USES: 
    use ESMF
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,     intent(IN) :: i   ! index of the observation type
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of a 
! land surface model output (from a LIS simulation) as observations.  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real                    :: run_dd(8)
    integer                 :: t
    integer                 :: ts
    type(ESMF_Config)       :: modelSpecConfig
    character*20            :: domain
    character*10            :: time
    character*20            :: obstype
    integer                 :: rc

    if(.not.allocated(lis6output)) then 
       allocate(lis6output(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config,lis6output(i)%odir, &
         label="LIS6 output directory:",rc=rc)
    call LVT_verify(rc,'LIS6 output directory: not defined')

    ts = 86400 !daily
    call LVT_update_timestep(LVT_rc, ts)

  end subroutine LIS6outputInit
  
end module LIS6outputMod
