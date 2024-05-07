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
! !MODULE: FLUXCOM_obsMod
! \label(FLUXCOM_obsMod)
!
! !INTERFACE:
module FLUXCOM_obsMod
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
!  5 April 2018   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: FLUXCOM_obsinit !Initializes structures for reading FLUXCOM data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: FLUXCOMobs !Object to hold FLUXCOM observation attributes
!EOP

  type, public :: fluxcomdec
     character*100           :: odir
     character*10            :: method
     real, allocatable       :: rlat(:)
     real, allocatable       :: rlon(:)
     integer, allocatable    :: n11(:)
     logical                 :: startflag
  end type fluxcomdec
     
  type(fluxcomdec), allocatable :: FLUXCOMObs(:)

contains
  
!BOP
! 
! !ROUTINE: FLUXCOM_obsInit
! \label{FLUXCOM_obsInit}
!
! !INTERFACE: 
  subroutine FLUXCOM_obsinit(i)
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
!   for reading the FLUXCOM data, including the computation of spatial 
!   interpolation weights. The FLUXCOM data is provides in the 
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

    if(.not.allocated(FLUXCOMobs)) then 
       allocate(FLUXCOMobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, FLUXCOMobs(I)%odir, &
         label='FLUXCOM data directory: ',rc=status)
    call LVT_verify(status, 'FLUXCOM data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, FLUXCOMobs(I)%method, &
         label='FLUXCOM data product method: ',rc=status)

    !ANN - Artificial Neural Networks
    !RF -  Random Forest
    !MARS - Multivariate Adaptive Regression Splines

    if(status.ne.0) then 
       write(LVT_logunit,*) '[ERR] FLUXCOM data product method: not defined'
       write(LVT_logunit,*) '[ERR] options are..'
       write(LVT_logunit,*) "[ERR] 'ANN', 'RF' or 'MARS'"
       call LVT_endrun()
    endif

    allocate(FLUXCOMobs(I)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(FLUXCOMobs(I)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(FLUXCOMobs(I)%n11(LVT_rc%lnc*LVT_rc%lnr))

    gridDesci = 0
    
    !filling the items needed by the interpolation library
    gridDesci(1) = 0  !input is EASE grid
    gridDesci(2) = 720
    gridDesci(3) = 360
    gridDesci(4) = -89.75
    gridDesci(5) = -179.75
    gridDesci(7) = 89.75
    gridDesci(8) = 179.75
    gridDesci(6) = 128
    gridDesci(9) = 0.50
    gridDesci(10) = 0.50
    gridDesci(20) = 64
 
    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         FLUXCOMobs(I)%rlat, FLUXCOMobs(I)%rlon,&
         FLUXCOMobs(I)%n11)

    FLUXCOMobs(I)%startflag = .true. 

    call LVT_update_timestep(LVT_rc, 86400)

  end subroutine FLUXCOM_obsinit


end module FLUXCOM_obsMod
