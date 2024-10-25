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
! !MODULE: NASA_AMSREsm_obsMod
! \label(NASA_AMSREsm_obsMod)
!
! !INTERFACE:
module NASA_AMSREsm_obsMod
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
!  This module handles the observation plugin for the standard NASA
!  AMSR-E soil moisture retrieval product
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: NASA_AMSREsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: NASA_AMSREsmobs
!EOP
  type, public :: nasaamsresmobsdec

     character*100        :: odir

     integer         :: mo
     integer,allocatable :: n112(:)
     real,allocatable    :: rlat2(:)
     real,allocatable    :: rlon2(:)
     real            :: gridDesci(50)

     real,    allocatable    :: smobs(:,:)
     real,    allocatable    :: smtime(:,:)
     integer*2, allocatable  :: smqc(:,:)     
     logical             :: startflag     

  end type nasaamsresmobsdec

  type(nasaamsresmobsdec), allocatable:: NASA_AMSREsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: NASA_AMSREsm_obsInit
! \label{NASA_AMSREsm_obsInit}
!
! !INTERFACE: 
  subroutine NASA_AMSREsm_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading NASA AMSR-E soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer         ezlh_convert
    integer            :: npts
    integer,parameter  :: ease_nr=586
    integer,parameter  :: ease_nc=1383
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc

    if(.not.allocated(NASA_AMSREsmobs)) then 
       allocate(NASA_AMSREsmobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, NASA_AMSREsmobs(i)%odir, &
         label='NASA AMSR-E soil moisture observation directory:', rc=status)
    call LVT_verify(status, 'NASA AMSR-E soil moisture observation directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    NASA_AMSREsmobs(i)%gridDesci=0.0 
    
    !filling the items needed by the interpolation library
    NASA_AMSREsmobs(i)%gridDesci(1) = 9  !input is EASE grid
    !these  corner coordinates were calculated based on ezlh_convert
    NASA_AMSREsmobs(i)%gridDesci(4) = -90.0  !lat
    NASA_AMSREsmobs(i)%gridDesci(5) = -179.6096 !lon
    NASA_AMSREsmobs(i)%gridDesci(7) = 83.33788  !lat
    NASA_AMSREsmobs(i)%gridDesci(8) = 180.1301  !lon
    
    
    NASA_AMSREsmobs(i)%gridDesci(2) = ease_nc  !nx
    NASA_AMSREsmobs(i)%gridDesci(3) = ease_nr  !ny
    
    npts= LVT_rc%lnc*LVT_rc%lnr
    NASA_AMSREsmobs(i)%mo=npts
    allocate(NASA_AMSREsmobs(i)%rlat2(npts))
    allocate(NASA_AMSREsmobs(i)%rlon2(npts))
    allocate(NASA_AMSREsmobs(i)%n112(npts))
    
    NASA_AMSREsmobs(i)%rlat2=0.0
    NASA_AMSREsmobs(i)%rlon2=0.0
    NASA_AMSREsmobs(i)%n112=0.0
    call neighbor_interp_input(NASA_AMSREsmobs(i)%gridDesci,LVT_rc%gridDesc,&
         npts,NASA_AMSREsmobs(i)%rlat2,&
         NASA_AMSREsmobs(i)%rlon2,NASA_AMSREsmobs(i)%n112)
    

    NASA_AMSREsmobs(i)%startflag = .true. 

    allocate(NASA_AMSREsmobs(i)%smobs(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(NASA_AMSREsmobs(i)%smtime(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(NASA_AMSREsmobs(i)%smqc(LVT_rc%lnc*LVT_rc%lnr,2))

!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------

  end subroutine NASA_AMSREsm_obsinit


end module NASA_AMSREsm_obsMod
