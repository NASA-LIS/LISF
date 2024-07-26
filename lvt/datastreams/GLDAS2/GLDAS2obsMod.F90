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
! !MODULE: GLDAS2obsMod
! \label(GLDAS2obsMod)
!
! !INTERFACE:
module GLDAS2obsMod
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
!  7 Mar 2015   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GLDAS2obsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GLDAS2obs !Object to hold GLDAS2 observation attributes
!EOP

  type, public :: gldas2dec
     character*100           :: odir
     character*20            :: model_name
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: yr
     integer                 :: mo
     logical                 :: startFlag
     real                    :: datares
  end type gldas2dec
     
  type(gldas2dec), allocatable :: GLDAS2Obs(:)

contains
  
!BOP
! 
! !ROUTINE: GLDAS2obsInit
! \label{GLDAS2obsInit}
!
! !INTERFACE: 
  subroutine GLDAS2obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the GIMMS NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    if(.not.allocated(GLDAS2obs)) then 
       allocate(GLDAS2obs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, GLDAS2obs(i)%odir, &
         label='GLDAS2 data directory:',rc=status)
    call LVT_verify(status, 'GLDAS2 data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, GLDAS2obs(i)%model_name, &
         label='GLDAS2 data model name:',rc=status)
    call LVT_verify(status, 'GLDAS2 data model name: not defined')

    call LVT_update_timestep(LVT_rc, 2592000)

    gldas2obs(i)%gridDesc = 0
        
    gldas2obs(i)%nc = 1440
    gldas2obs(i)%nr = 600

    !filling the items needed by the interpolation library
    gldas2obs(i)%gridDesc(1) = 0  
    gldas2obs(i)%gridDesc(2) = gldas2obs(i)%nc
    gldas2obs(i)%gridDesc(3) = gldas2obs(i)%nr
    gldas2obs(i)%gridDesc(4) = -59.875
    gldas2obs(i)%gridDesc(5) = -179.875
    gldas2obs(i)%gridDesc(7) = 89.875
    gldas2obs(i)%gridDesc(8) = 179.875
    gldas2obs(i)%gridDesc(6) = 128
    gldas2obs(i)%gridDesc(9) = 0.25
    gldas2obs(i)%gridDesc(10) = 0.25
    gldas2obs(i)%gridDesc(20) = 64

    gldas2obs(i)%datares  = 0.25

    if(LVT_isAtAfinerResolution(gldas2obs(i)%datares)) then
       
       allocate(gldas2obs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gldas2obs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(gldas2obs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(gldas2obs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            gldas2obs(i)%rlat, &
            gldas2obs(i)%rlon, &
            gldas2obs(i)%n11)
    else
       allocate(gldas2obs(i)%n11(gldas2obs(i)%nc*gldas2obs(i)%nr))
       call upscaleByAveraging_input(gldas2obs(i)%gridDesc,&
            LVT_rc%gridDesc,gldas2obs(i)%nc*gldas2obs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,gldas2obs(i)%n11)
    endif

    gldas2obs(i)%yr = -1
    gldas2obs(i)%mo = LVT_rc%mo
    gldas2obs(i)%startFlag = .false. 

  end subroutine GLDAS2obsinit


end module GLDAS2obsMod
