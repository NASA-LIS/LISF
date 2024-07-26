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
! !MODULE: ALEXI_obsMod
! \label(ALEXI_obsMod)
!
! !INTERFACE:
module ALEXI_obsMod
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
  PUBLIC :: ALEXI_obsinit !Initializes structures for reading ALEXI data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ALEXIobs !Object to hold ALEXI observation attributes
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
     
  type(alexidec), allocatable :: ALEXIObs(:)

contains
  
!BOP
! 
! !ROUTINE: ALEXI_obsInit
! \label{ALEXI_obsInit}
!
! !INTERFACE: 
  subroutine ALEXI_obsinit(source)
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
!   for reading the ALEXI data, including the computation of spatial 
!   interpolation weights. The ALEXI data is provides in the 
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


    if(.not.allocated(ALEXIobs)) then 
       allocate(ALEXIobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, ALEXIobs(source)%odir, &
         label='ALEXI data directory: ',rc=status)
    call LVT_verify(status, 'ALEXI data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, ALEXIobs(source)%extent, &
         label='ALEXI data domain extent: ',rc=status)

    if(status.ne.0) then 
       write(LVT_logunit,*) '[ERR] ALEXI data domain extent: not defined'
       write(LVT_logunit,*) '[ERR] supported options are: "CONUS" or "GLOBAL"'
       call LVT_endrun()
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, ALEXIobs(source)%res, &
         label='ALEXI data resolution (in km): ',rc=status)
    call LVT_verify(status, 'ALEXI data resolution (in km): not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    allocate(ALEXIobs(source)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ALEXIobs(source)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    gridDesci = 0
    
    if(ALEXIobs(source)%extent.eq."CONUS") then 
       if(ALEXIobs(source)%res.eq.10) then 
          !filling the items needed by the interpolation library
          gridDesci(1) = 0  
          gridDesci(2) = 667
          gridDesci(3) = 291
          gridDesci(4) = 23.9584
          gridDesci(5) = -124.9474
          gridDesci(7) = 50.038709
          gridDesci(8) = -65.0526214
          gridDesci(6) = 128
          gridDesci(9) = 0.0899321
          gridDesci(10) = 0.0899321
          gridDesci(20) = 64
          
          ALEXIobs(source)%nc = 667
          ALEXIobs(source)%nr = 291  

          ALEXIobs(source)%datares = 0.0899321
          
       elseif(ALEXIobs(source)%res.eq.4) then 
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

          ALEXIobs(source)%datares = 0.04

          ALEXIobs(source)%nc = 1456
          ALEXIobs(source)%nr = 625
       else
          write(LVT_logunit,*) "[ERR] The ALEXI data resolution must be either 10 or 4"
          call LVT_endrun()
       endif
    elseif(ALEXIobs(source)%extent.eq."GLOBAL") then 
       gridDesci(1) = 0  
       gridDesci(2) = 7200
       gridDesci(3) = 3000
       gridDesci(4) = -59.975
       gridDesci(5) = -179.975
       gridDesci(7) = 89.975
       gridDesci(8) = 179.975
       gridDesci(6) = 128
       gridDesci(9) = 0.05
       gridDesci(10) = 0.05
       gridDesci(20) = 64
       
       ALEXIobs(source)%nc = 7200
       ALEXIobs(source)%nr = 3000

       ALEXIobs(source)%datares = 0.05
    endif
      
    if(LVT_isAtAfinerResolution(ALEXIobs(source)%datares)) then
       allocate(ALEXIobs(source)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(gridDesci,LVT_rc%gridDesc,&
            LVT_rc%lnc*LVT_rc%lnr, &
            ALEXIobs(source)%rlat, ALEXIobs(source)%rlon,&
            ALEXIobs(source)%n11)
    else
       allocate(ALEXIobs(source)%n11(&
            ALEXIobs(source)%nc*ALEXIobs(source)%nr))
       call upscaleByAveraging_input(gridDesci,&
            LVT_rc%gridDesc,ALEXIobs(source)%nc*ALEXIobs(source)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,ALEXIobs(source)%n11)       
    endif

    ALEXIobs(source)%startFlag = .true. 

  end subroutine ALEXI_obsinit

end module ALEXI_obsMod
