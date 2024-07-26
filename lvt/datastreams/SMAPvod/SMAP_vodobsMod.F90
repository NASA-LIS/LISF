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
! !MODULE: SMAP_vodobsMod
! \label(SMAP_vodobsMod)
!
! !INTERFACE:
module SMAP_vodobsMod
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
!  vegetation optical depth 
!  retrieval products from the Soil Moisture Active Passive
!  (SMAP) mission
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Apr 2019: Sujay Kumar, Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_vodobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAP_vodobs
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

     real,    allocatable    :: vodobs(:,:)
     real,    allocatable    :: vodtime(:,:)
     integer*2, allocatable  :: vodqc(:,:)     
     logical                 :: startflag     

  end type smapobsdec

  type(smapobsdec), allocatable:: SMAP_vodobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMAP_vodobsInit
! \label{SMAP_vodobsInit}
!
! !INTERFACE: 
  subroutine SMAP_vodobsinit(i)
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

    if(.not.allocated(SMAP_vodobs)) then 
       allocate(SMAP_vodobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_vodobs(i)%odir, &
         label='SMAP vegetation optical depth observation directory:', rc=status)
    call LVT_verify(status, &
         'SMAP vegetation optical depth observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMAP_vodobs(i)%data_designation, &
         label='SMAP vegetation optical depth data designation:', rc=status)
    call LVT_verify(status, &
         'SMAP vegetation optical depth data designation: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    if(SMAP_vodobs(i)%data_designation.eq."SPL3SMAP") then 
       !SMAP L3 radar/radiometer daily 9km
       SMAP_vodobs(i)%gridDesci=0.0 

       ease_nr=1624
       ease_nc=3856       

    !filling the items needed by the interpolation library
       SMAP_vodobs(i)%gridDesci(1) = 9  !input is EASE grid
    !these  corner coordinates were calculated based on ezlh_convert
       SMAP_vodobs(i)%gridDesci(4) = -84.6564  !lat
       SMAP_vodobs(i)%gridDesci(5) = -179.953 !lon
       SMAP_vodobs(i)%gridDesci(7) = 84.6564  !lat
       SMAP_vodobs(i)%gridDesci(8) = 179.953  !lon
       SMAP_vodobs(i)%gridDesci(9) = 5 !M09 grid
       SMAP_vodobs(i)%gridDesci(20) = 64
       
       SMAP_vodobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_vodobs(i)%gridDesci(3) = ease_nr  !ny

       SMAP_vodobs(i)%nc=ease_nc
       SMAP_vodobs(i)%nr=ease_nr       
    elseif(SMAP_vodobs(i)%data_designation.eq."SPL3SMP") then 
       !SMAP L3 radiometer daily 36km 
    !filling the items needed by the interpolation library
       ease_nc=964       
       ease_nr=406
       SMAP_vodobs(i)%gridDesci = 0
       SMAP_vodobs(i)%gridDesci(1) = 9  !input is EASE grid
       SMAP_vodobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_vodobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_vodobs(i)%gridDesci(9) = 4 !M36 grid
       SMAP_vodobs(i)%gridDesci(20) = 64
       SMAP_vodobs(i)%gridDesci(10) = 0.36
       SMAP_vodobs(i)%gridDesci(11) = 1

       SMAP_vodobs(i)%nc=ease_nc
       SMAP_vodobs(i)%nr=ease_nr       
    elseif(SMAP_vodobs(i)%data_designation.eq."SPL3SMP_E") then 
       ease_nc=3856       
       ease_nr=1624
       SMAP_vodobs(i)%gridDesci = 0
       SMAP_vodobs(i)%gridDesci(1) = 9  
       SMAP_vodobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_vodobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_vodobs(i)%gridDesci(9) = 5 !M09 grid
       SMAP_vodobs(i)%gridDesci(20) = 64
       SMAP_vodobs(i)%gridDesci(10) = 0.09
       SMAP_vodobs(i)%gridDesci(11) = 1

       SMAP_vodobs(i)%nc=ease_nc
       SMAP_vodobs(i)%nr=ease_nr    
 elseif(SMAP_vodobs(i)%data_designation.eq."SPL2SMP") then 
       ease_nc=964       
       ease_nr=406
       SMAP_vodobs(i)%gridDesci = 0
       SMAP_vodobs(i)%gridDesci(1) = 9  
       SMAP_vodobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_vodobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_vodobs(i)%gridDesci(9) = 4 !M36 grid
       SMAP_vodobs(i)%gridDesci(20) = 64
       SMAP_vodobs(i)%gridDesci(10) = 0.36
       SMAP_vodobs(i)%gridDesci(11) = 1

       SMAP_vodobs(i)%nc=ease_nc
       SMAP_vodobs(i)%nr=ease_nr

       call LVT_update_timestep(LVT_rc, 3600)   
    elseif(SMAP_vodobs(i)%data_designation.eq."SPL2SMP_E") then 
       ease_nc=3856       
       ease_nr=1624
       SMAP_vodobs(i)%gridDesci = 0
       SMAP_vodobs(i)%gridDesci(1) = 9  
       SMAP_vodobs(i)%gridDesci(2) = ease_nc  !nx
       SMAP_vodobs(i)%gridDesci(3) = ease_nr  !ny
       SMAP_vodobs(i)%gridDesci(9) = 5 !M09 grid
       SMAP_vodobs(i)%gridDesci(20) = 64
       SMAP_vodobs(i)%gridDesci(10) = 0.09
       SMAP_vodobs(i)%gridDesci(11) = 1

       SMAP_vodobs(i)%nc=ease_nc
       SMAP_vodobs(i)%nr=ease_nr

       call LVT_update_timestep(LVT_rc, 3600)
       
    endif

    npts= LVT_rc%lnc*LVT_rc%lnr
    SMAP_vodobs(i)%mo=npts
    
!    if(LVT_isAtAfinerResolution(SMAP_vodobs(i)%gridDesci(10))) then 
       allocate(SMAP_vodobs(i)%rlat2(npts))
       allocate(SMAP_vodobs(i)%rlon2(npts))
       allocate(SMAP_vodobs(i)%n112(npts))
       
       SMAP_vodobs(i)%rlat2=0.0
       SMAP_vodobs(i)%rlon2=0.0
       SMAP_vodobs(i)%n112=0.0
       call neighbor_interp_input(SMAP_vodobs(i)%gridDesci,&
            LVT_rc%gridDesc,&
            npts,SMAP_vodobs(i)%rlat2,&
            SMAP_vodobs(i)%rlon2,SMAP_vodobs(i)%n112)
!    else
!       allocate(SMAP_vodobs(i)%n112(ease_nc*ease_nr))
!       call upscaleByAveraging_input(&
!            SMAP_vodobs(i)%gridDesci,&
!            LVT_rc%gridDesc,&
!            ease_nc*ease_nr,&
!            npts, &
!            SMAP_vodobs(i)%n112)       
!    endif

    SMAP_vodobs(i)%startflag = .true. 

    allocate(SMAP_vodobs(i)%vodobs(LVT_rc%lnc,LVT_rc%lnr))
    allocate(SMAP_vodobs(i)%vodtime(LVT_rc%lnc,LVT_rc%lnr))
    allocate(SMAP_vodobs(i)%vodqc(LVT_rc%lnc*LVT_rc%lnr,2))

    call system("mkdir -p "//trim('SMAPvod'))     

  end subroutine SMAP_vodobsinit


end module SMAP_vodobsMod
