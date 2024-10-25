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
! !MODULE: GLEAM_obsMod
! \label(GLEAM_obsMod)
!
! !INTERFACE:
module GLEAM_obsMod
! 
! !USES: 
  use ESMF
  use LVT_logMod

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
!  14 Feb 2017   Sujay Kumar  Initial Specification
! 
!EOP
!

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GLEAM_obsinit !Initializes structures for reading GLEAM data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GLEAMobs !Object to hold GLEAM observation attributes
!EOP

  type, public :: gleamdec
     character*100           :: odir
     character*10            :: version
     integer                 :: nc,nr
     logical                 :: startFlag
     real, allocatable       :: rlat(:)
     real, allocatable       :: rlon(:)
     integer, allocatable    :: n11(:)
  end type gleamdec
     
  type(gleamdec), allocatable :: GLEAMObs(:)

contains
  
!BOP
! 
! !ROUTINE: GLEAM_obsInit
! \label{GLEAM_obsInit}
!
! !INTERFACE: 
  subroutine GLEAM_obsinit(source)
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
!   for reading the GLEAM data, including the computation of spatial 
!   interpolation weights. The GLEAM data is provides in the 
!   EASE grid projection. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer               :: source
    integer               :: status
    real                  :: gridDesci(50)

    if(.not.allocated(GLEAMobs)) then 
       allocate(GLEAMobs(LVT_rc%nDataStreams))
    endif

    write(LVT_logunit,*) "[INFO] DATA STREAM:  GLEAM "

    call ESMF_ConfigGetAttribute(LVT_Config, GLEAMobs(source)%odir, &
         label='GLEAM data directory: ',rc=status)
    call LVT_verify(status, 'GLEAM data directory: not defined')

! Enter 3.0a/3.0b/3.0c

    call ESMF_ConfigGetAttribute(LVT_Config, GLEAMobs(source)%version, &
         label='GLEAM data version: ',rc=status)
    call LVT_verify(status, 'GLEAM data version: not defined')


    call LVT_update_timestep(LVT_rc, 86400)

    allocate(GLEAMobs(source)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(GLEAMobs(source)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(GLEAMobs(source)%n11(LVT_rc%lnc*LVT_rc%lnr))

    gridDesci = 0
    
    gridDesci(1) = 0  
    gridDesci(2) = 1440
    gridDesci(3) = 720
    gridDesci(4) = -89.875
    gridDesci(5) = -179.875
    gridDesci(7) = 89.875
    gridDesci(8) = 179.875
    gridDesci(6) = 128
    gridDesci(9) = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64
    
    GLEAMobs(source)%nc = 1440
    GLEAMobs(source)%nr = 720      
 
    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         GLEAMobs(source)%rlat, GLEAMobs(source)%rlon,&
         GLEAMobs(source)%n11)

    GLEAMobs(source)%startFlag = .true. 

  end subroutine GLEAM_obsinit

end module GLEAM_obsMod
