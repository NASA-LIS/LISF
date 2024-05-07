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
! !MODULE: GCOMW_AMSR2L3sm_obsMod
! \label(GCOMW_AMSR2L3sm_obsMod)
!
! !INTERFACE:
module GCOMW_AMSR2L3sm_obsMod
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
!  This module handles the observation plugin for the Land Parameter
!  Retrieval Model (LPRM) AMSR-E soil moisture product
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  12 Dec 2014: Sujay Kumar, Initial Specification
! 
!EOP
! 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GCOMW_AMSR2L3sm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GCOMW_AMSR2L3smobs
!EOP
  type, public :: amsr2smobsdec

     character*100          :: odir
     integer                :: mo
     logical                :: startmode
     integer                :: amsr2nc, amsr2nr
     type(proj_info)        :: amsr2proj
     real                   :: datares
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

  end type amsr2smobsdec

  type(amsr2smobsdec), allocatable:: GCOMW_AMSR2L3smobs(:)

contains
  
!BOP
! 
! !ROUTINE: GCOMW_AMSR2L3sm_obsInit
! \label{GCOMW_AMSR2L3sm_obsInit}
!
! !INTERFACE: 
  subroutine GCOMW_AMSR2L3sm_obsinit(i)
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
!  for reading LPRM AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts, status
    real                  :: gridDesci(50)

    if(.not.allocated(GCOMW_AMSR2L3smobs)) then 
       allocate(GCOMW_AMSR2L3smobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, GCOMW_AMSR2L3smobs(i)%odir, &
         label='GCOMW AMSR2 L3 soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'GCOMW AMSR2 L3 soil moisture observation directory: not defined')
   
    call LVT_update_timestep(LVT_rc, 86400)

    GCOMW_AMSR2L3smobs(i)%startmode = .true. 

    GCOMW_AMSR2L3smobs(i)%amsr2nc = 3600
    GCOMW_AMSR2L3smobs(i)%amsr2nr = 1800
    
    call map_set(PROJ_LATLON, -89.95,-179.95,&
         0.0, 0.10,0.10, 0.0,&
         GCOMW_AMSR2L3smobs(i)%amsr2nc,GCOMW_AMSR2L3smobs(i)%amsr2nr,&
         GCOMW_AMSR2L3smobs(i)%amsr2proj)
    
    gridDesci = 0
    gridDesci(1) = 0
    gridDesci(2) = 3600
    gridDesci(3) = 1800
    gridDesci(4) = -89.95
    gridDesci(5) = -179.95
    gridDesci(6) = 128
    gridDesci(7) = 89.95
    gridDesci(8) = 179.95
    gridDesci(9) = 0.10
    gridDesci(10) = 0.10
    gridDesci(20) = 64
    
    GCOMW_AMSR2L3smobs(i)%datares = 0.10
    
    if(LVT_isatAfinerResolution(GCOMW_AMSR2L3smobs(i)%datares)) then
       allocate(GCOMW_AMSR2L3smobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GCOMW_AMSR2L3smobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       
       allocate(GCOMW_AMSR2L3smobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GCOMW_AMSR2L3smobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GCOMW_AMSR2L3smobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GCOMW_AMSR2L3smobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       
       allocate(GCOMW_AMSR2L3smobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GCOMW_AMSR2L3smobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GCOMW_AMSR2L3smobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(GCOMW_AMSR2L3smobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
       
       call bilinear_interp_input(gridDesci,LVT_rc%gridDesc(:),&
            LVT_rc%lnc*LVT_rc%lnr,GCOMW_AMSR2L3smobs(i)%rlat, &
            GCOMW_AMSR2L3smobs(i)%rlon,GCOMW_AMSR2L3smobs(i)%n11, &
            GCOMW_AMSR2L3smobs(i)%n12, GCOMW_AMSR2L3smobs(i)%n21, &
            GCOMW_AMSR2L3smobs(i)%n22, GCOMW_AMSR2L3smobs(i)%w11, &
            GCOMW_AMSR2L3smobs(i)%w12, GCOMW_AMSR2L3smobs(i)%w21, &
            GCOMW_AMSR2L3smobs(i)%w22)
    else
       allocate(GCOMW_AMSR2L3smobs(i)%n11(GCOMW_AMSR2L3smobs(i)%amsr2nc*&
            GCOMW_AMSR2L3smobs(i)%amsr2nr))
       call upscaleByAveraging_input(gridDesci,&
            LVT_rc%gridDesc(:),&
            GCOMW_AMSR2L3smobs(i)%amsr2nc*GCOMW_AMSR2L3smobs(i)%amsr2nr,&
            LVT_rc%lnc*LVT_rc%lnr,&
            GCOMW_AMSR2L3smobs(i)%n11)
    endif

  end subroutine GCOMW_AMSR2L3sm_obsinit


end module GCOMW_AMSR2L3sm_obsMod
