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
! !MODULE: UASMAP_obsMod
! \label(UASMAP_obsMod)
!
! !INTERFACE:
module UASMAP_obsMod
! 
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the
!  UASMAP soil moisture product
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  5 Apr 2021: Sujay Kumar, Initial Specification
! 
!EOP
! 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: UASMAP_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: UASMAPobs
!EOP
  type, public :: thysmobsdec
     character*100              :: odir
     real                       :: gridDesci(50)
     integer                    :: nc, nr
     type(proj_info)            :: proj
     integer, allocatable       :: n11(:)
     integer, allocatable       :: n12(:)
     integer, allocatable       :: n21(:)
     integer, allocatable       :: n22(:)
     real,  allocatable         :: w11(:)
     real,  allocatable         :: w12(:)
     real,  allocatable         :: w21(:)
     real,  allocatable         :: w22(:)
     real,  allocatable         :: rlat(:)
     real,  allocatable         :: rlon(:)
     logical                    :: startmode     

     type(ESMF_TimeInterval)    :: dt

  end type thysmobsdec

  type(thysmobsdec), allocatable:: UASMAPobs(:)

contains
  
!BOP
! 
! !ROUTINE: UASMAP_obsInit
! \label{UASMAP_obsInit}
!
! !INTERFACE: 
  subroutine UASMAP_obsinit(i)
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
!  for reading UA SMAP soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts, status
    real                  :: gridDesci(50)

    real                  :: depth(4)

    if(.not.allocated(UASMAPobs)) then 
       allocate(UASMAPobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, UASMAPobs(i)%odir, &
         label='UA SMAP data directory:', rc=status)
    call LVT_verify(status, 'UA SMAP data directory: not defined')
    
    call LVT_update_timestep(LVT_rc, 86400)

    UASMAPobs(i)%startmode = .true. 

    UASMAPobs(i)%nc = 6946
    UASMAPobs(i)%nr = 2988
    
    call map_set(PROJ_LATLON,24.5243432,-124.76103,&
         0.0, 0.0083228, 0.0083228,0.0,&
         UASMAPobs(i)%nc,UASMAPobs(i)%nr,&
         UASMAPobs(i)%proj)

    gridDesci = 0 
    gridDesci(1) = 0 
    gridDesci(2) = 6946
    gridDesci(3) = 2988
    gridDesci(4) = 24.52434324
    gridDesci(5) = -124.7610304
    gridDesci(6) = 128
    gridDesci(7) = 49.38818468
    gridDesci(8) = -66.95072617
    gridDesci(9) = 0.008322282
    gridDesci(10) = 0.008322282
    gridDesci(20) = 64
    
    allocate(UASMAPobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(UASMAPobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(UASMAPobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(UASMAPobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(UASMAPobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(UASMAPobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(UASMAPobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(UASMAPobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(UASMAPobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(UASMAPobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,UASMAPobs(i)%rlat, &
         UASMAPobs(i)%rlon,UASMAPobs(i)%n11, &
         UASMAPobs(i)%n12, UASMAPobs(i)%n21, &
         UASMAPobs(i)%n22, UASMAPobs(i)%w11, &
         UASMAPobs(i)%w12, UASMAPobs(i)%w21, &
         UASMAPobs(i)%w22)

    
    call ESMF_TimeIntervalSet(UASMAPobs(i)%dt,s=86400,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: UASMAP_obsinit')
    

  end subroutine UASMAP_obsinit


end module UASMAP_obsMod
