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
! !MODULE: JASMINsm_obsMod
! \label(JASMINsm_obsMod)
!
! !INTERFACE:
module JASMINsm_obsMod
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
!  21 July 2010: Sujay Kumar, Initial Specification
! 
!EOP
! 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: JASMINsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: JASMINsmobs
!EOP
  type, public :: jasminsmobsdec
     character*100              :: odir
     real                       :: gridDesci(50)
     integer                    :: jasminnc, jasminnr
     type(proj_info)            :: jasminproj
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

     type(ESMF_Time)            :: refTime
     type(ESMF_TimeInterval)    :: dt
     real                       :: sf_wt(4), rz_wt(4)
  end type jasminsmobsdec

  type(jasminsmobsdec), allocatable:: JASMINsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: JASMINsm_obsInit
! \label{JASMINsm_obsInit}
!
! !INTERFACE: 
  subroutine JASMINsm_obsinit(i)
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

    if(.not.allocated(JASMINsmobs)) then 
       allocate(JASMINsmobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, JASMINsmobs(i)%odir, &
         label='JASMIN soil moisture data directory:', rc=status)
    call LVT_verify(status, 'JASMIN soil moisture data directory: not defined')
    
    call LVT_update_timestep(LVT_rc, 86400)

    JASMINsmobs(i)%startmode = .true. 

    JASMINsmobs(i)%jasminnc = 811
    JASMINsmobs(i)%jasminnr = 669
    
    call map_set(PROJ_LATLON, -43.95,113.15,&
         0.0, 0.05,0.05, 0.0,&
         JASMINsmobs(i)%jasminnc,JASMINsmobs(i)%jasminnr,&
         JASMINsmobs(i)%jasminproj)
    
    gridDesci = 0 
    gridDesci(1) = 0 
    gridDesci(2) = 811
    gridDesci(3) = 669
    gridDesci(4) = -43.95
    gridDesci(5) = 113.15
    gridDesci(6) = 128
    gridDesci(7) = -10.55
    gridDesci(8) = 153.65
    gridDesci(9) = 0.05
    gridDesci(10) = 0.05
    gridDesci(20) = 64
    
    allocate(JASMINsmobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(JASMINsmobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(JASMINsmobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(JASMINsmobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(JASMINsmobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(JASMINsmobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(JASMINsmobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(JASMINsmobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(JASMINsmobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(JASMINsmobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,JASMINsmobs(i)%rlat, &
         JASMINsmobs(i)%rlon,JASMINsmobs(i)%n11, &
         JASMINsmobs(i)%n12, JASMINsmobs(i)%n21, &
         JASMINsmobs(i)%n22, JASMINsmobs(i)%w11, &
         JASMINsmobs(i)%w12, JASMINsmobs(i)%w21, &
         JASMINsmobs(i)%w22)

    
    call ESMF_TimeSet(JASMINsmobs(i)%reftime, yy=2010,&
         mm = 1,&
         dd = 1,&
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status,'error in timeset: JASMINsm_obsinit')
    
    call ESMF_TimeIntervalSet(JASMINsmobs(i)%dt,s=86400,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: JASMINsm_obsinit')
    
    depth(1) = 0.1
    depth(2) = 0.25
    depth(3) = 0.65
    depth(4) = 2.0

    call compute_vinterp_weights(&
         4,LVT_rc%lis_sf_d,&
         LVT_rc%lis_rz_d,&
         depth(1:4),&
         JASMINsmobs(i)%sf_wt,&
         JASMINsmobs(i)%rz_wt)

  end subroutine JASMINsm_obsinit


end module JASMINsm_obsMod
