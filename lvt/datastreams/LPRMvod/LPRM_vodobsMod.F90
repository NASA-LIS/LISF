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
! !MODULE: LPRM_vodobsMod
! \label(LPRM_vodobsMod)
!
! !INTERFACE:
module LPRM_vodobsMod
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
!  retrieval products from the Land Parameter Retrieval Model 
!  (LPRM) 
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
  PUBLIC :: LPRM_vodobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LPRM_vodobs
!EOP
  type, public :: lprmobsdec

     character*100        :: odir
     character*20         :: data_designation

     integer              :: nc
     integer              :: nr
     integer              :: mo
     integer,allocatable  :: n11(:)
     integer,allocatable  :: n12(:)
     integer,allocatable  :: n21(:)
     integer,allocatable  :: n22(:)

     real,allocatable     :: rlat(:)
     real,allocatable     :: rlon(:)

     real,allocatable     :: w11(:)
     real,allocatable     :: w12(:)
     real,allocatable     :: w21(:)
     real,allocatable     :: w22(:)

     logical                 :: startflag     

  end type lprmobsdec

  type(lprmobsdec), allocatable:: LPRM_vodobs(:)

contains
  
!BOP
! 
! !ROUTINE: LPRM_vodobsInit
! \label{LPRM_vodobsInit}
!
! !INTERFACE: 
  subroutine LPRM_vodobsinit(i)
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

    real                 :: gridDesci(50)
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc

    if(.not.allocated(LPRM_vodobs)) then 
       allocate(LPRM_vodobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, LPRM_vodobs(i)%odir, &
         label='LPRM vegetation optical depth observation directory:', rc=status)
    call LVT_verify(status, &
         'LPRM vegetation optical depth observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, LPRM_vodobs(i)%data_designation, &
         label='LPRM vegetation optical depth data designation:', rc=status)
    if(status.ne.0) then 
       write(LVT_logunit,*) "[ERR] 'LPRM vegetation optical depth data designation:' not defined"
       write(LVT_logunit,*) "[ERR] supported options are 'C-band' or 'X-band'"
       call LVT_endrun()
    endif

    call LVT_update_timestep(LVT_rc, 86400)

    gridDesci=0.0 
    LPRM_vodobs(i)%nc = 1440
    LPRM_vodobs(i)%nr = 720

    !filling the items needed by the interpolation library
    gridDesci(1) = 0
    gridDesci(2) = 1440
    gridDesci(3) = 720
    gridDesci(4) = -89.875
    gridDesci(5) = -179.875
    gridDesci(6) = 128 
    gridDesci(7) = 89.875
    gridDesci(8) = 179.875
    gridDesci(9) = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64
    
    npts= LVT_rc%lnc*LVT_rc%lnr
    LPRM_vodobs(i)%mo=npts
    
    allocate(LPRM_vodobs(i)%rlat(npts))
    allocate(LPRM_vodobs(i)%rlon(npts))
    allocate(LPRM_vodobs(i)%n11(npts))
    allocate(LPRM_vodobs(i)%n12(npts))
    allocate(LPRM_vodobs(i)%n21(npts))
    allocate(LPRM_vodobs(i)%n22(npts))

    allocate(LPRM_vodobs(i)%w11(npts))
    allocate(LPRM_vodobs(i)%w12(npts))
    allocate(LPRM_vodobs(i)%w21(npts))
    allocate(LPRM_vodobs(i)%w22(npts))
      
    call bilinear_interp_input(gridDesci,&
         LVT_rc%gridDesc,&
         npts,&
         LPRM_vodobs(i)%rlat,&
         LPRM_vodobs(i)%rlon,&
         LPRM_vodobs(i)%n11,&
         LPRM_vodobs(i)%n12,&
         LPRM_vodobs(i)%n21,&
         LPRM_vodobs(i)%n22,&
         LPRM_vodobs(i)%w11,&
         LPRM_vodobs(i)%w12,&
         LPRM_vodobs(i)%w21,&
         LPRM_vodobs(i)%w22)

    LPRM_vodobs(i)%startflag = .true. 


  end subroutine LPRM_vodobsinit

end module LPRM_vodobsMod
