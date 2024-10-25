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
! !MODULE: THySM_obsMod
! \label(THySM_obsMod)
!
! !INTERFACE:
module THySM_obsMod
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
!  THySM soil moisture product
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
  PUBLIC :: THySM_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: THySMobs
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

  type(thysmobsdec), allocatable:: THySMobs(:)

contains
  
!BOP
! 
! !ROUTINE: THySM_obsInit
! \label{THySM_obsInit}
!
! !INTERFACE: 
  subroutine THySM_obsinit(i)
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
!  for reading JASMIN AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts, status
    real                  :: gridDesci(50)

    real                  :: depth(4)

    if(.not.allocated(THySMobs)) then 
       allocate(THySMobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, THySMobs(i)%odir, &
         label='THySM data directory:', rc=status)
    call LVT_verify(status, 'THySM data directory: not defined')
    
    call LVT_update_timestep(LVT_rc, 86400)

    THySMobs(i)%startmode = .true. 

    THySMobs(i)%nc = 5800
    THySMobs(i)%nr = 2800
    
    call map_set(PROJ_LATLON, 25.005,-124.995,&
         0.0, 0.01,0.01, 0.0,&
         THySMobs(i)%nc,THySMobs(i)%nr,&
         THySMobs(i)%proj)
    
    gridDesci = 0 
    gridDesci(1) = 0 
    gridDesci(2) = 5800
    gridDesci(3) = 2800
    gridDesci(4) = 25.005
    gridDesci(5) = -124.995
    gridDesci(6) = 128
    gridDesci(7) = 52.995
    gridDesci(8) = -67.005
    gridDesci(9) = 0.01
    gridDesci(10) = 0.01
    gridDesci(20) = 64
    
    allocate(THySMobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(THySMobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(THySMobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(THySMobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(THySMobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(THySMobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(THySMobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(THySMobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(THySMobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(THySMobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,THySMobs(i)%rlat, &
         THySMobs(i)%rlon,THySMobs(i)%n11, &
         THySMobs(i)%n12, THySMobs(i)%n21, &
         THySMobs(i)%n22, THySMobs(i)%w11, &
         THySMobs(i)%w12, THySMobs(i)%w21, &
         THySMobs(i)%w22)

    
    call ESMF_TimeIntervalSet(THySMobs(i)%dt,s=86400,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: THySM_obsinit')
    

  end subroutine THySM_obsinit


end module THySM_obsMod
