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
! !MODULE: FluxSat_obsMod
! \label(FluxSat_obsMod)
!
! !INTERFACE:
module FluxSat_obsMod
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
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Feb 2021   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: FluxSat_obsinit !Initializes structures for reading FluxSat data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: FluxSatobs !Object to hold FluxSat observation attributes
!EOP

  type, public :: fluxcomdec
     character*100           :: odir
     real, allocatable       :: rlat(:)
     real, allocatable       :: rlon(:)
     integer, allocatable    :: n11(:)
     logical                 :: startflag
  end type fluxcomdec
     
  type(fluxcomdec), allocatable :: FluxSatObs(:)

contains
  
!BOP
! 
! !ROUTINE: FluxSat_obsInit
! \label{FluxSat_obsInit}
!
! !INTERFACE: 
  subroutine FluxSat_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the FluxSat data, including the computation of spatial 
!   interpolation weights. The FluxSat data is provides in the 
!   EASE grid projection. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: i 
    integer               :: status
    real                  :: gridDesci(50)

    if(.not.allocated(FluxSatobs)) then 
       allocate(FluxSatobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, FluxSatobs(I)%odir, &
         label='FluxSat data directory: ',rc=status)
    call LVT_verify(status, 'FluxSat data directory: not defined')

    allocate(FluxSatobs(I)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(FluxSatobs(I)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(FluxSatobs(I)%n11(LVT_rc%lnc*LVT_rc%lnr))

    gridDesci = 0
    
    !filling the items needed by the interpolation library
    gridDesci(1) = 0  !input is EASE grid
    gridDesci(2) = 7200
    gridDesci(3) = 3600
    gridDesci(4) = -89.975
    gridDesci(5) = -179.975
    gridDesci(7) = 89.975
    gridDesci(8) = 179.975
    gridDesci(6) = 128
    gridDesci(9) = 0.05
    gridDesci(10) = 0.05
    gridDesci(20) = 64
 
    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         FluxSatobs(I)%rlat, FluxSatobs(I)%rlon,&
         FluxSatobs(I)%n11)

    FluxSatobs(I)%startflag = .true. 

    call LVT_update_timestep(LVT_rc, 86400)

  end subroutine FluxSat_obsinit


end module FluxSat_obsMod
