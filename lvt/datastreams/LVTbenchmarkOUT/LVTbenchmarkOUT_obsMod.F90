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
! !MODULE: LVTbenchmarkOUT_obsMod
! \label(LVTbenchmarkOUT_obsMod)
!
! !INTERFACE:
module LVTbenchmarkOUT_obsMod
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
  PUBLIC :: LVTbenchmarkOUT_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lvtbenchobs
!
!EOP
  
  type, public :: daobsdec
     character*100  :: odir
     character*50   :: vname
  end type daobsdec

  type(daobsdec), allocatable  :: lvtbenchobs(:)

contains

!BOP
! 
! !ROUTINE: LVTbenchmarkOUT_obsInit
! \label{LVTbenchmarkOUT_obsInit}
!
! !INTERFACE: 
  subroutine LVTbenchmarkOUT_obsInit(i)
! 
! !USES: 
    use ESMF
    use LVT_coreMod,    only : LVT_rc, LVT_config
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod,     only : LVT_verify

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

    integer                 :: t 
    type(ESMF_Config)       :: modelSpecConfig
    character*20            :: domain
    integer                 :: rc

    if(.not.allocated(lvtbenchobs)) then 
       allocate(lvtbenchobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config,lvtbenchobs(i)%odir, &
         label="LVT benchmark output directory:",rc=rc)
    call LVT_verify(rc,'LVT benchmark output directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,lvtbenchobs(i)%vname, &
         label="LVT benchmark variable:",rc=rc)
    call LVT_verify(rc,'LVT benchmark variable: not defined')

!Right now hardcoded
    call LVT_update_timestep(LVT_rc, 1800)
    
  end subroutine LVTbenchmarkOUT_obsInit
  
end module LVTbenchmarkOUT_obsMod
