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
! !MODULE: FMISWE_obsMod
! \label(FMISWE_obsMod)
!
! !INTERFACE:
module FMISWE_obsMod
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
!  This module handles the observation plugin for the Finnish Meteorological
!  Institute snow course data. The snow course measurements are made only 
!  two times a month. Here we use the daily snow water equivalent values 
!  interpolated on 10kmx10km grid based on original snow water equivalent 
!  measurements and weather data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Apr 2009   Sujay Kumar  Initial Specification
! 
!EOP
!
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: FMISWE_obsinit !Initializes structures for reading FMI SWE data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: FMISWEobs !Object to hold FMISWE observation attributes
!EOP
  type, public :: fmisweobsdec
     character*100        :: odir
     real                 :: udef
     logical              :: startflag
     integer              :: yr
     integer              :: mo
     integer              :: da
     integer              :: nc, nr
     integer, allocatable        :: n11(:)
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
  end type fmisweobsdec

  type(fmisweobsdec), allocatable :: fmisweobs(:)

contains
  
!BOP
! 
! !ROUTINE: FMISWE_obsInit
! \label{FMISWE_obsInit}
!
! !INTERFACE: 
  subroutine FMISWE_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod
    use map_utils

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading FMISWE data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer            :: status
    real               :: gridDesci(50)
    
    if(.not.allocated(fmisweobs)) then 
       allocate(fmisweobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, fmisweobs(i)%odir, &
         label='FMISWE observation directory:',rc=status)
    call LVT_verify(status, 'FMISWE observation directory: not defined')
    
    call LVT_update_timestep(LVT_rc, 86400)

    fmisweobs(i)%nc = 131
    fmisweobs(i)%nr = 107

    gridDesci = 0 
    gridDesci(1) = 0
    gridDesci(2) = fmisweobs(i)%nc
    gridDesci(3) = fmisweobs(i)%nr
    griddesci(4) = 59.45
    gridDesci(5) = 19.15
    gridDesci(6) = 128
    gridDesci(7) = 70.05
    gridDesci(8) = 32.15
    gridDesci(9) = 0.10
    gridDesci(10) =0.10
    gridDesci(20) = 64

    allocate(fmisweobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(fmisweobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    allocate(fmisweobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))

    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr,fmisweobs(i)%rlat,fmisweobs(i)%rlon,&
         fmisweobs(i)%n11)

    fmisweobs(i)%startflag = .true. 
   
  end subroutine FMISWE_obsinit


end module FMISWE_obsMod
