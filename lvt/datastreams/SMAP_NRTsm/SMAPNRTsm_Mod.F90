!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.6
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: SMAPNRTsm_Mod
! \label(SMAPNRTsm_Mod)
!
! !INTERFACE:
module SMAPNRTsm_Mod
  !
  ! !USES:
  use ESMF
  use LVT_constantsMod, only: LVT_CONST_PATH_LEN

  implicit none

  PRIVATE
  !
  ! !DESCRIPTION:
  !  This module handles the observation plugin for the
  !  retrieval products from the Soil Moisture Active Passive
  !  (SMAP) mission
  !
  ! !FILES USED:
  !
  ! !REVISION HISTORY:
  !  2 July 2025: Mahdi Navari, Initial Specification
  !
  !EOP

  !------------------------------------------------------------------------
  ! !PUBLIC MEMBER FUNCTIONS:
  !------------------------------------------------------------------------

  PUBLIC :: SMAPNRT_smobsinit

  !------------------------------------------------------------------------
  ! !PUBLIC TYPES:
  !------------------------------------------------------------------------

  PUBLIC :: SMAPNRT_smobs

  !EOP

  type, public :: smapobsdec

     character(LVT_CONST_PATH_LEN) :: odir
     character*20         :: data_designation
     character*3          :: release_number
     integer              :: nc
     integer              :: nr
     integer              :: mo
     integer, allocatable :: n112(:)
     real, allocatable    :: rlat2(:)
     real, allocatable    :: rlon2(:)
     real                 :: gridDesci(50)

     real, allocatable    :: smobs(:,:)
     real, allocatable    :: smtime(:,:)
     integer*2, allocatable  :: smqc(:,:)
     logical                 :: startflag

  end type smapobsdec

  type(smapobsdec), allocatable:: SMAPNRT_smobs(:)

contains

  !BOP
  !
  ! !ROUTINE: SMAPNRT_smobsInit
  ! \label{SMAPNRT_smobsInit}
  !
  ! !INTERFACE:
  subroutine SMAPNRT_smobsinit(i)
    !
    ! !USES:
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
    !
    ! !INPUT PARAMETERS:
    integer,   intent(IN) :: i
    !
    ! !OUTPUT PARAMETERS:
    !
    ! !DESCRIPTION:
    !  This subroutine initializes and sets up the data structures
    !  required for reading SMAP L2 NRT soil moisture data.
    !
    ! !FILES USED:
    !
    ! !REVISION HISTORY:
    !
    !EOP

    integer :: npts
    integer :: ease_nc, ease_nr
    integer :: status

    external :: neighbor_interp_input
    external :: system

    if (.not.allocated(SMAPNRT_smobs)) then
       allocate(SMAPNRT_smobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMAPNRT_smobs(i)%odir, &
         label='SMAP L2 NRT soil moisture observation directory:', &
         rc=status)
    call LVT_verify(status, &
         'SMAP L2 NRT soil moisture observation directory: not defined')

    !SMAP L2 36km
    !filling the items needed by the interpolation library
    ease_nc = 964
    ease_nr = 406
    SMAPNRT_smobs(i)%gridDesci = 0
    SMAPNRT_smobs(i)%gridDesci(1) = 9  !input is EASE grid
    SMAPNRT_smobs(i)%gridDesci(2) = ease_nc  !nx
    SMAPNRT_smobs(i)%gridDesci(3) = ease_nr  !ny
    SMAPNRT_smobs(i)%gridDesci(9) = 4 !M36 grid
    SMAPNRT_smobs(i)%gridDesci(20) = 64
    SMAPNRT_smobs(i)%gridDesci(10) = 0.36
    SMAPNRT_smobs(i)%gridDesci(11) = 1

    SMAPNRT_smobs(i)%nc=ease_nc
    SMAPNRT_smobs(i)%nr=ease_nr

    call LVT_update_timestep(LVT_rc, 3600)

    npts= LVT_rc%lnc*LVT_rc%lnr
    SMAPNRT_smobs(i)%mo=npts

    allocate(SMAPNRT_smobs(i)%rlat2(npts))
    allocate(SMAPNRT_smobs(i)%rlon2(npts))
    allocate(SMAPNRT_smobs(i)%n112(npts))

    SMAPNRT_smobs(i)%rlat2 = 0.0
    SMAPNRT_smobs(i)%rlon2 = 0.0
    SMAPNRT_smobs(i)%n112 = 0.0
    call neighbor_interp_input(SMAPNRT_smobs(i)%gridDesci, &
         LVT_rc%gridDesc, &
         npts, SMAPNRT_smobs(i)%rlat2, &
         SMAPNRT_smobs(i)%rlon2, SMAPNRT_smobs(i)%n112)

    SMAPNRT_smobs(i)%startflag = .true.

    allocate(SMAPNRT_smobs(i)%smobs(LVT_rc%lnc*LVT_rc%lnr,2))
    allocate(SMAPNRT_smobs(i)%smtime(LVT_rc%lnc,LVT_rc%lnr))
    allocate(SMAPNRT_smobs(i)%smqc(LVT_rc%lnc*LVT_rc%lnr,2))

    call system("mkdir -p "//trim('SMAPsm'))

  end subroutine SMAPNRT_smobsinit
end module SMAPNRTsm_Mod
