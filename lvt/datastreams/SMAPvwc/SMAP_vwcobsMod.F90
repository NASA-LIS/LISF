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
! !MODULE: SMAP_vwcobsMod
! \label(SMAP_vwcobsMod)
!
! !INTERFACE:
module SMAP_vwcobsMod
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
!  (SMAP) mission
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  28 Aug 2018: Mahdi Navari,  Sujay Kumar Initial Specification
!  31 July 2019 Mahdi Navari : SMAP Composite Release ID was added (this option asks a user to 
!         enter the part of Composite Release ID a three-character string like R16 ) 
!
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_vwcobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_vwcobs
!EOP
  type, public :: smapobsdec

     character*100        :: odir
     character*20         :: data_designation
     character*3           :: release_number
     integer              :: nc
     integer              :: nr
     integer              :: mo
     integer,allocatable  :: n112(:)
     real,allocatable     :: rlat2(:)
     real,allocatable     :: rlon2(:)
     real                 :: gridDesci(50)

     real,    allocatable    :: vwcobs(:,:)
     real,    allocatable    :: vwctime(:,:)
     integer*2, allocatable  :: vwcqc(:,:)     
     logical                 :: startflag     

  end type smapobsdec

  type(smapobsdec), allocatable:: SMAP_vwcobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMAP_vwcobsInit
! \label{SMAP_vwcobsInit}
!
! !INTERFACE: 
  subroutine SMAP_vwcobsinit(i)
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
!  for reading NASA vegetation water content data. 
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

    if(.not.allocated(SMAP_vwcobs)) then 
       allocate(SMAP_vwcobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_vwcobs(i)%odir, &
         label='SMAP soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'SMAP soil moisture observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_vwcobs(i)%data_designation, &
         label='SMAP soil moisture data designation:', rc=status)
    call LVT_verify(status, &
         'SMAP soil moisture data designation: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_vwcobs(i)%release_number, &
         label='SMAP(NASA) soil moisture Composite Release ID (e.g., R16):', rc=status)
    call LVT_verify(status, &
         'SMAP(NASA) soil moisture Composite Release ID (e.g., R16): not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    if(SMAP_vwcobs(i)%data_designation.eq."SPL3SMAP") then 
       !SMAP L3 radar/radiometer daily 9km
       SMAP_vwcobs(i)%gridDesci=0.0 

       ease_nr=1624
       ease_nc=3856       

    !filling the items needed by the interpolation library
       SMAP_vwcobs(i)%gridDesci(1) = 9  !input is EASE grid
    !these  corner coordinates were calculated based on ezlh_convert
       SMAP_vwcobs(i)%gridDesci(4) = -84.6564  !lat
       SMAP_vwcobs(i)%gridDesci(5) = -179.953 !lon
       SMAP_vwcobs(i)%gridDesci(7) = 84.6564  !lat
       SMAP_vwcobs(i)%gridDesci(8) = 179.953  !lon
       SMAP_vwcobs(i)%gridDesci(9) = 5 !M09 grid
       SMAP_vwcobs(i)%gridDesci(20) = 64
       
       SMAP_vwcobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_vwcobs(i)%gridDesci(3) = ease_nr  !ny

       SMAP_vwcobs(i)%nc=ease_nc
       SMAP_vwcobs(i)%nr=ease_nr       
    elseif(SMAP_vwcobs(i)%data_designation.eq."SPL3SMP") then 
       !SMAP L3 radiometer daily 36km 
    !filling the items needed by the interpolation library
       ease_nc=964       
       ease_nr=406
       SMAP_vwcobs(i)%gridDesci = 0
       SMAP_vwcobs(i)%gridDesci(1) = 9  !input is EASE grid
       SMAP_vwcobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_vwcobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_vwcobs(i)%gridDesci(9) = 4 !M36 grid
       SMAP_vwcobs(i)%gridDesci(20) = 64
       SMAP_vwcobs(i)%gridDesci(10) = 0.36
       SMAP_vwcobs(i)%gridDesci(11) = 1

       SMAP_vwcobs(i)%nc=ease_nc
       SMAP_vwcobs(i)%nr=ease_nr       
    elseif(SMAP_vwcobs(i)%data_designation.eq."SPL3SMP_E") then 
       ease_nc=3856       
       ease_nr=1624
       SMAP_vwcobs(i)%gridDesci = 0
       SMAP_vwcobs(i)%gridDesci(1) = 9  
       SMAP_vwcobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_vwcobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_vwcobs(i)%gridDesci(9) = 5 !M09 grid
       SMAP_vwcobs(i)%gridDesci(20) = 64
       SMAP_vwcobs(i)%gridDesci(10) = 0.09
       SMAP_vwcobs(i)%gridDesci(11) = 1

       SMAP_vwcobs(i)%nc=ease_nc
       SMAP_vwcobs(i)%nr=ease_nr       
    endif

    npts= LVT_rc%lnc*LVT_rc%lnr
    SMAP_vwcobs(i)%mo=npts
    
    
    allocate(SMAP_vwcobs(i)%rlat2(npts))
    allocate(SMAP_vwcobs(i)%rlon2(npts))
    allocate(SMAP_vwcobs(i)%n112(npts))
    
    SMAP_vwcobs(i)%rlat2=0.0
    SMAP_vwcobs(i)%rlon2=0.0
    SMAP_vwcobs(i)%n112=0.0
    call neighbor_interp_input(SMAP_vwcobs(i)%gridDesci,LVT_rc%gridDesc,&
         npts,SMAP_vwcobs(i)%rlat2,&
         SMAP_vwcobs(i)%rlon2,SMAP_vwcobs(i)%n112)
    

    SMAP_vwcobs(i)%startflag = .true. 

    allocate(SMAP_vwcobs(i)%vwcobs(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMAP_vwcobs(i)%vwctime(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMAP_vwcobs(i)%vwcqc(LVT_rc%lnc*LVT_rc%lnr,2))

!-------------------------------------------------------------------------
!  AMSRE data contains the a top soil soil moisture data
!-------------------------------------------------------------------------

  end subroutine SMAP_vwcobsinit


end module SMAP_vwcobsMod
