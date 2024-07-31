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
! !MODULE: ALEXIesi_obsMod
! \label(ALEXIesi_obsMod)
!
! !INTERFACE:
module ALEXIesi_obsMod
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
!  This module handles the observation plugin for the 
!  Evaporative Stress Index (ESI) estimates from 
!  the Atmospheric Land Exchange Inverse (ALEXI) model
!  (Anderson et al. 1997, A two-source time-integrated
!  model for estimating surface fluxes using thermal
!  remote sensing, RSE)
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  08 Feb 2020:   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ALEXIesi_obsinit !Initializes structures for reading ALEXI ESI data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ALEXIesiobs !Object to hold ALEXI ESI observation attributes
!EOP

  type, public :: alexidec
     character*100           :: odir
     character*50            :: extent
     integer                 :: nc,nr
     integer                 :: res
     logical                 :: startFlag
     real                    :: datares
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
     integer, allocatable        :: n11(:)
  end type alexidec
     
  type(alexidec), allocatable :: ALEXIesiObs(:)

contains
  
!BOP
! 
! !ROUTINE: ALEXIesi_obsInit
! \label{ALEXIesi_obsInit}
!
! !INTERFACE: 
  subroutine ALEXIesi_obsinit(source)
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
!   for reading the ALEXI ESI data, including the computation of spatial 
!   interpolation weights. The current plugin supports the data over
!   CONUS (in lat/lon projection) at 4 km spatial resolution. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer               :: source
    integer               :: status
    real                  :: gridDesci(50)


    if(.not.allocated(ALEXIesiobs)) then 
       allocate(ALEXIesiobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, ALEXIesiobs(source)%odir, &
         label='ALEXI ESI data directory: ',rc=status)
    call LVT_verify(status, 'ALEXI ESI data directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    allocate(ALEXIesiobs(source)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ALEXIesiobs(source)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    gridDesci = 0
  
    gridDesci(1) = 0  
    gridDesci(2) = 1456
    gridDesci(3) = 625
    gridDesci(4) = 24.80
    gridDesci(5) = -125.0
    gridDesci(7) = 49.76
    gridDesci(8) = -66.8
    gridDesci(6) = 128
    gridDesci(9) = 0.04
    gridDesci(10) = 0.04
    gridDesci(20) = 64
    
    ALEXIesiobs(source)%datares = 0.04
    
    ALEXIesiobs(source)%nc = 1456
    ALEXIesiobs(source)%nr = 625
      
    if(LVT_isAtAfinerResolution(ALEXIesiobs(source)%datares)) then
       allocate(ALEXIesiobs(source)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
            LVT_rc%lnc*LVT_rc%lnr, &
            ALEXIesiobs(source)%rlat, ALEXIesiobs(source)%rlon,&
            ALEXIesiobs(source)%n11)
    else
       allocate(ALEXIesiobs(source)%n11(&
            ALEXIesiobs(source)%nc*ALEXIesiobs(source)%nr))
       call upscaleByAveraging_input(gridDesci,&
            LVT_rc%gridDesc,ALEXIesiobs(source)%nc*ALEXIesiobs(source)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,ALEXIesiobs(source)%n11)       
    endif

    ALEXIesiobs(source)%startFlag = .true. 

  end subroutine ALEXIesi_obsinit

end module ALEXIesi_obsMod
