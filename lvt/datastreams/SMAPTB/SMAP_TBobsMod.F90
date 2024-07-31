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
! !MODULE: SMAP_TBobsMod
! \label(SMAP_TBobsMod)
!
! !INTERFACE:
module SMAP_TBobsMod
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
  PUBLIC :: SMAP_TBobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_TBobs
!EOP
  type, public :: smapobsdec

     character*100        :: odir
     character*20         :: data_designation

     integer              :: nc
     integer              :: nr
     integer              :: mo
     integer,allocatable  :: n112(:)
     real,allocatable     :: rlat2(:)
     real,allocatable     :: rlon2(:)
     real                 :: gridDesci(50)

     real,    allocatable    :: TBobs(:,:)
     real,    allocatable    :: smtime(:,:)
     integer*2, allocatable  :: smqc(:,:)     
     logical                 :: startflag     

  end type smapobsdec

  type(smapobsdec), allocatable:: SMAP_TBobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMAP_TBobsInit
! \label{SMAP_TBobsInit}
!
! !INTERFACE: 
  subroutine SMAP_TBobsinit(i)
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
    integer                 :: ease_nc,ease_nr
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc


    if(.not.allocated(SMAP_TBobs)) then 
       allocate(SMAP_TBobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_TBobs(i)%odir, &
         label='SMAP TB observation directory:', rc=status)
    call LVT_verify(status, &
         'SMAP TB observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_TBobs(i)%data_designation, &
         label='SMAP TB data designation:', rc=status)
    call LVT_verify(status, &
         'SMAP TB data designation: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    if(SMAP_TBobs(i)%data_designation.eq."SPL1CTB") then 
       !SMAP L3 radiometer daily 36km 
    !filling the items needed by the interpolation library
       SMAP_TBobs(i)%gridDesci(1) = 9  !input is EASE grid
       ease_nr=406
       ease_nc=964       
       
       SMAP_TBobs(i)%gridDesci(9) = 4 !M36 grid
       SMAP_TBobs(i)%gridDesci(20) = 64
       
       SMAP_TBobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_TBobs(i)%gridDesci(3) = ease_nr  !ny

       SMAP_TBobs(i)%nc=ease_nc
       SMAP_TBobs(i)%nr=ease_nr       
    endif

    npts= LVT_rc%lnc*LVT_rc%lnr
    SMAP_TBobs(i)%mo=npts
    
    
    allocate(SMAP_TBobs(i)%rlat2(npts))
    allocate(SMAP_TBobs(i)%rlon2(npts))
    allocate(SMAP_TBobs(i)%n112(npts))
    
    SMAP_TBobs(i)%rlat2=0.0
    SMAP_TBobs(i)%rlon2=0.0
    SMAP_TBobs(i)%n112=0.0
    call neighbor_interp_input(SMAP_TBobs(i)%gridDesci,LVT_rc%gridDesc,&
         npts,SMAP_TBobs(i)%rlat2,&
         SMAP_TBobs(i)%rlon2,SMAP_TBobs(i)%n112)
    

    SMAP_TBobs(i)%startflag = .true. 

    allocate(SMAP_TBobs(i)%TBobs(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMAP_TBobs(i)%smtime(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMAP_TBobs(i)%smqc(LVT_rc%lnc*LVT_rc%lnr,2))

!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------

  end subroutine SMAP_TBobsinit


end module SMAP_TBobsMod
