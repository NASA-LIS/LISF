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
! !MODULE: FLUXNETmte_obsMod
! \label(FLUXNETmte_obsMod)
!
! !INTERFACE:
module FLUXNETmte_obsMod
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
!  18 May 2011   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: FLUXNETmte_obsinit !Initializes structures for reading FLUXNETmte data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: FLUXNETmteobs !Object to hold FLUXNETmte observation attributes
!EOP

  type, public :: fluxnetdec
     character*100           :: odir
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
     integer, allocatable        :: n11(:)
     real,    allocatable        :: qle(:,:)
     real,    allocatable        :: qh(:,:)
     integer                 :: yr
     integer                 :: mo
  end type fluxnetdec
     
  type(fluxnetdec), allocatable :: FLUXNETmteObs(:)

contains
  
!BOP
! 
! !ROUTINE: FLUXNETmte_obsInit
! \label{FLUXNETmte_obsInit}
!
! !INTERFACE: 
  subroutine FLUXNETmte_obsinit(i)
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
!   for reading the FLUXNETmte data, including the computation of spatial 
!   interpolation weights. The FLUXNETmte data is provides in the 
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

    if(.not.allocated(FLUXNETmteobs)) then 
       allocate(FLUXNETmteobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, FLUXNETmteobs(I)%odir, &
         label='FLUXNETmte data directory: ',rc=status)
    call LVT_verify(status, 'FLUXNETmte data directory: not defined')

    allocate(FLUXNETmteobs(I)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(FLUXNETmteobs(I)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(FLUXNETmteobs(I)%n11(LVT_rc%lnc*LVT_rc%lnr))

    allocate(FLUXNETmteobs(I)%qle(LVT_rc%lnc*LVT_rc%lnr,12))
    allocate(FLUXNETmteobs(I)%qh(LVT_rc%lnc*LVT_rc%lnr,12))

    gridDesci = 0
    
    !filling the items needed by the interpolation library
    gridDesci(1) = 0  !input is EASE grid
    gridDesci(2) = 720
    gridDesci(3) = 291
    gridDesci(4) = -55.25
    gridDesci(5) = -179.75
    gridDesci(7) = 89.75
    gridDesci(8) = 179.75
    gridDesci(6) = 128
    gridDesci(9) = 0.50
    gridDesci(10) = 0.50
    gridDesci(20) = 64
 
    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         FLUXNETmteobs(I)%rlat, FLUXNETmteobs(I)%rlon,&
         FLUXNETmteobs(I)%n11)

    FLUXNETmteobs(I)%yr = -1
    FLUXNETmteobs(I)%mo = LVT_rc%mo

    call LVT_update_timestep(LVT_rc, 2592000)

  end subroutine FLUXNETmte_obsinit


end module FLUXNETmte_obsMod
