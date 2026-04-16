!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: WSFsm_obsMod
!
! !DESCRIPTION:
!   This module contains interfaces and subroutines to
!   handle WSF soil moisture retrievals for DA preprocessing
!   (CDF generation) in LDT.
!
! !REVISION HISTORY:
!   2025: Ehsan Jalilvand; Initial Specification
!
module WSFsm_obsMod
! !USES:
  use ESMF
  use map_utils
  use LDT_constantsMod, only: LDT_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: WSFsmobsinit  !Initializes structures for reading WSF SM data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: WSFsmobs      !Object to hold WSF SM observation attributes
!EOP
  type, public :: wsfsmobs_dec

     character(len=LDT_CONST_PATH_LEN) :: odir
     integer                :: nc, nr
     real                   :: gridDesci(50)
     real,    allocatable   :: smobs(:,:)
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real,    allocatable   :: w11(:)
     real,    allocatable   :: w12(:)
     real,    allocatable   :: w21(:)
     real,    allocatable   :: w22(:)

  end type wsfsmobs_dec

  type(wsfsmobs_dec), allocatable :: WSFsmobs(:)

contains

!BOP
!
! !ROUTINE: WSFsmobsinit
! \label{WSFsmobsinit}
!
! !INTERFACE:
  subroutine WSFsmobsinit()
!
! !USES:
    use ESMF
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !ARGUMENTS:

!
! !DESCRIPTION:
!   This subroutine initializes and sets up the data structures required
!   for reading the WSF soil moisture data for CDF generation.
!
!EOP

    integer          :: status
    integer          :: n

    external :: bilinear_interp_input
    external :: upscaleByAveraging_input

    allocate(WSFsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'WSF soil moisture data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, WSFsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'WSF soil moisture data directory: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(WSFsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       WSFsmobs(n)%smobs = -9999.0

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
             "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1

       WSFsmobs(n)%nc = 2560
       WSFsmobs(n)%nr = 1920

       WSFsmobs(n)%gridDesci(1) = 0
       WSFsmobs(n)%gridDesci(2) = WSFsmobs(n)%nc
       WSFsmobs(n)%gridDesci(3) = WSFsmobs(n)%nr
       WSFsmobs(n)%gridDesci(4) = -89.953125
       WSFsmobs(n)%gridDesci(5) = -179.929688
       WSFsmobs(n)%gridDesci(6) = 128
       WSFsmobs(n)%gridDesci(7) = 89.953125
       WSFsmobs(n)%gridDesci(8) = 179.929688
       WSFsmobs(n)%gridDesci(9) = 0.140570   !dlon
       WSFsmobs(n)%gridDesci(10) = 0.093701  !dlat
       WSFsmobs(n)%gridDesci(20) = 64

       if(LDT_isLDTatAfinerResolution(n,0.093701)) then

          allocate(WSFsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WSFsmobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WSFsmobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WSFsmobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WSFsmobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WSFsmobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WSFsmobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(WSFsmobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n, &
               WSFsmobs(n)%gridDesci,&
               WSFsmobs(n)%n11,&
               WSFsmobs(n)%n12,&
               WSFsmobs(n)%n21,&
               WSFsmobs(n)%n22,&
               WSFsmobs(n)%w11,&
               WSFsmobs(n)%w12,&
               WSFsmobs(n)%w21,&
               WSFsmobs(n)%w22)

       else

          allocate(WSFsmobs(n)%n11(WSFsmobs(n)%nc*&
               WSFsmobs(n)%nr))

          call upscaleByAveraging_input(&
               WSFsmobs(n)%gridDesci,&
               LDT_rc%gridDesc(n,:),&
               WSFsmobs(n)%nc*&
               WSFsmobs(n)%nr,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               WSFsmobs(n)%n11)

       endif
    enddo
  end subroutine WSFsmobsinit

end module WSFsm_obsMod
