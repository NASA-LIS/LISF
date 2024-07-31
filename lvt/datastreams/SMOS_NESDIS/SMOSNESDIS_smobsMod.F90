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
! !MODULE: SMOSNESDIS_smobsMod
! \label(SMOSNESDIS_smobsMod)
!
! !INTERFACE:
module SMOSNESDIS_smobsMod
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
!  retrieval products from the Soil Moisture Active Passive
!  (SMOSNESDIS) mission
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  21 July 2016: Sujay Kumar, Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSNESDIS_smobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSNESDIS_smobs
!EOP
  type, public :: smosnesdisobsdec

     character*100        :: odir

     integer              :: nc
     integer              :: nr
     integer              :: mo
     integer,allocatable  :: n112(:)
     real,allocatable     :: rlat2(:)
     real,allocatable     :: rlon2(:)
     real                 :: gridDesci(50)

     real,    allocatable    :: smobs(:,:)
     logical                 :: startflag     

  end type smosnesdisobsdec

  type(smosnesdisobsdec), allocatable:: SMOSNESDIS_smobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMOSNESDIS_smobsInit
! \label{SMOSNESDIS_smobsInit}
!
! !INTERFACE: 
  subroutine SMOSNESDIS_smobsinit(i)
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

    if(.not.allocated(SMOSNESDIS_smobs)) then 
       allocate(SMOSNESDIS_smobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMOSNESDIS_smobs(i)%odir, &
         label='SMOS (NESDIS) soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'SMOS (NESDIS) soil moisture observation directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    SMOSNESDIS_smobs(i)%gridDesci=0.0 

    SMOSNESDIS_smobs(i)%nc = 1440
    SMOSNESDIS_smobs(i)%nr = 720

    SMOSNESDIS_smobs(i)%gridDesci(1) = 0
    SMOSNESDIS_smobs(i)%gridDesci(2) = 1440
    SMOSNESDIS_smobs(i)%gridDesci(3) = 720
    SMOSNESDIS_smobs(i)%gridDesci(4) = -89.875
    SMOSNESDIS_smobs(i)%gridDesci(5) = -179.875
    SMOSNESDIS_smobs(i)%gridDesci(6) = 128
    SMOSNESDIS_smobs(i)%gridDesci(7) = 89.875
    SMOSNESDIS_smobs(i)%gridDesci(8) = 179.875
    SMOSNESDIS_smobs(i)%gridDesci(9) =  0.25
    SMOSNESDIS_smobs(i)%gridDesci(10) = 0.25
    SMOSNESDIS_smobs(i)%gridDesci(20) = 64

    npts= LVT_rc%lnc*LVT_rc%lnr
    SMOSNESDIS_smobs(i)%mo=npts
    
    allocate(SMOSNESDIS_smobs(i)%rlat2(npts))
    allocate(SMOSNESDIS_smobs(i)%rlon2(npts))
    allocate(SMOSNESDIS_smobs(i)%n112(npts))
    
    SMOSNESDIS_smobs(i)%rlat2=0.0
    SMOSNESDIS_smobs(i)%rlon2=0.0
    SMOSNESDIS_smobs(i)%n112=0.0
    call neighbor_interp_input(SMOSNESDIS_smobs(i)%gridDesci,LVT_rc%gridDesc,&
         npts,SMOSNESDIS_smobs(i)%rlat2,&
         SMOSNESDIS_smobs(i)%rlon2,SMOSNESDIS_smobs(i)%n112)
    

    SMOSNESDIS_smobs(i)%startflag = .true. 

    allocate(SMOSNESDIS_smobs(i)%smobs(LVT_rc%lnc*LVT_rc%lnr,2))

!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------

  end subroutine SMOSNESDIS_smobsinit


end module SMOSNESDIS_smobsMod
