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
! !MODULE: LISDAdiagOutputMod
! \label(LISDAdiagOutputMod)
!
! !INTERFACE:
module LISDAdiagOutputMod
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
  PUBLIC :: LISDAdiagOutputInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: lisdadiagoutput
!
!EOP
  
  type, public :: lisdadiagobsdec
     character*100  :: odir
     integer        :: computeInnovDist
     integer        :: computeGain
     integer        :: computeSpread
     integer        :: computeAnlIncr
     integer        :: nstvars
     integer        :: instance

  end type lisdadiagobsdec

  type(lisdadiagobsdec), allocatable  :: lisdadiagoutput(:)

contains

!BOP
! 
! !ROUTINE: LISDAdiagOutputInit
! \label{LISDAdiagOutputInit}
!
! !INTERFACE: 
  subroutine LISDAdiagOutputInit(i)
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

    if(.not.allocated(lisdadiagoutput)) then 
       allocate(lisdadiagoutput(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config,lisdadiagoutput(i)%odir, &
         label="LIS DA output directory:",rc=rc)
    call LVT_verify(rc,'LIS DA output directory: not defined')

     call ESMF_ConfigGetAttribute(LVT_config,&
          lisdadiagoutput(i)%computeInnovDist, &
          label="LIS DA process innovation distribution:", rc=rc)
     call LVT_verify(rc, 'LIS DA process innovation distribution: not defined')
     call ESMF_ConfigGetAttribute(LVT_config,&
          lisdadiagoutput(i)%computeGain, &
          label="LIS DA process analysis gain:", rc=rc)
     call LVT_verify(rc, 'LIS DA process analysis gain: not defined')
     call ESMF_ConfigGetAttribute(LVT_config,&
          lisdadiagoutput(i)%computeSpread, &
          label="LIS DA process ensemble spread:", rc=rc)
     call LVT_verify(rc, 'LIS DA process ensemble spread: not defined')

     call ESMF_ConfigGetAttribute(LVT_config,&
          lisdadiagoutput(i)%computeAnlIncr, &
          label="LIS DA process analysis increments:", rc=rc)
     call LVT_verify(rc, 'LIS DA process analysis increments: not defined')

     call ESMF_ConfigGetAttribute(LVT_config,&
          lisdadiagoutput(i)%nstvars, &
          label="LIS DA Number of state variables in the DA update:", rc=rc)
     call LVT_verify(rc, &
          'LIS DA Number of state variables in the DA update: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,lisdadiagoutput(i)%instance, &
         label="LIS DA instance index:",rc=rc)
    call LVT_verify(rc,'LIS DA instance index: not defined')

    ts = 86400 !daily
    call LVT_update_timestep(LVT_rc, ts)

  end subroutine LISDAdiagOutputInit
  
end module LISDAdiagOutputMod
