!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: SMOSNRTNNL2sm_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! SMOS NRT NN L2 soil moisture retrievals
!
!
!   
! !REVISION HISTORY: 
!  31 Dec 2020: Yonghwan Kwon, Initial Specification
!  21 Feb. 2021: Mahdi Navari, code modified to write the DGG 
!                lookup table into a netCDF file
!
module SMOSNRTNNL2sm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSNRTNNL2sm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSNRTNNsmobs
!EOP
  type, public :: SMOS_in_lis_gridbox
     integer, allocatable :: dgg_indices(:)
     logical              :: dgg_assign
  end type

  type, public :: smosnrtnnl2smdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     real                   :: search_radius
     integer                :: mo
     real,    allocatable   :: smobs(:,:)
     integer                :: nc, nr
     type(proj_info)        :: proj
     !integer, allocatable   :: n11(:)
     integer                :: start_day, count_day
     integer, allocatable   :: dgg_lookup_1d(:)

     type(SMOS_in_lis_gridbox), pointer :: SMOS_lookup(:,:)

  end type smosnrtnnl2smdec

  type(smosnrtnnl2smdec), allocatable:: SMOSNRTNNsmobs(:)

contains

!BOP
! 
! !ROUTINE: SMOSNRTNNsm_obsInit
! \label{SMOSNRTNNsm_obsInit}
! 
! !INTERFACE: 
  subroutine SMOSNRTNNL2sm_obsinit()

! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading SMOS NRT NN L2 soil moisture data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n

    allocate(SMOSNRTNNsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOS NRT NN soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOSNRTNNsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'SMOS NRT NN soil moisture observation directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOS NRT NN search radius for openwater proximity detection:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOSNRTNNsmobs(n)%search_radius, &
            rc=status)
       call LDT_verify(status, &
            'SMOS NRT NN search radius for openwater proximity detection: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(SMOSNRTNNsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       allocate(SMOSNRTNNsmobs(n)%SMOS_lookup(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       SMOSNRTNNsmobs(n)%smobs = -9999.0
       SMOSNRTNNsmobs(n)%SMOS_lookup%dgg_assign = .false.
       SMOSNRTNNsmobs(n)%start_day = LDT_rc%da
       SMOSNRTNNsmobs(n)%count_day = 0
       SMOSNRTNNsmobs(n)%dgg_lookup_1d = 0
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1

    enddo
  end subroutine SMOSNRTNNL2sm_obsinit

end module SMOSNRTNNL2sm_obsMod
