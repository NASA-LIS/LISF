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
!   Oct 2026: The grid of the WSF SM files is read from ldt.config
!             ("WSF soil moisture grid ..." entries). All six entries
!             are required; there is no default grid.
!
module WSFsm_obsMod
! !USES:
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
     real                   :: lat_ll, lon_ll   ! lower-left cell centre
     real                   :: dlat, dlon       ! grid spacing (deg)
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
    integer          :: n, i
    real, allocatable :: grid_dd(:,:)
    character(len=48), parameter :: grid_labels(6) = [character(len=48) :: &
         'WSF soil moisture grid lower left lat:', &
         'WSF soil moisture grid lower left lon:', &
         'WSF soil moisture grid upper right lat:', &
         'WSF soil moisture grid upper right lon:', &
         'WSF soil moisture grid resolution (dx):', &
         'WSF soil moisture grid resolution (dy):' ]

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

    ! Grid of the WSF soil moisture files (cell centres). All six
    ! entries are required; there is no default grid.
    allocate(grid_dd(LDT_rc%nnest,6))
    do i = 1, 6
       call ESMF_ConfigFindLabel(LDT_config, trim(grid_labels(i)), rc=status)
       if (status /= 0) then
          write(LDT_logunit,*) '[ERR] ', trim(grid_labels(i)), &
               ' not defined in ldt.config'
          write(LDT_logunit,*) '[ERR] Give the grid of the WSF soil ', &
               'moisture files with all six entries:'
          write(LDT_logunit,*) '[ERR]   ', trim(grid_labels(1)), ' ... ', &
               trim(grid_labels(6))
       endif
       call LDT_verify(status, trim(grid_labels(i))//' not defined')
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config, grid_dd(n,i), rc=status)
          call LDT_verify(status, trim(grid_labels(i))//' not defined')
       enddo
    enddo
    do n=1,LDT_rc%nnest
       if (grid_dd(n,5) <= 0.0 .or. grid_dd(n,6) <= 0.0 .or. &
            grid_dd(n,3) <= grid_dd(n,1) .or. &
            grid_dd(n,4) <= grid_dd(n,2)) then
          write(LDT_logunit,*) '[ERR] Invalid WSF soil moisture grid in ', &
               'ldt.config (lower left, upper right, dx, dy):'
          write(LDT_logunit,*) '[ERR] ', grid_dd(n,:)
          call LDT_endrun()
       endif
    enddo

    do n=1,LDT_rc%nnest

       allocate(WSFsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       WSFsmobs(n)%smobs = -9999.0

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
             "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1

       WSFsmobs(n)%lat_ll = grid_dd(n,1)
       WSFsmobs(n)%lon_ll = grid_dd(n,2)
       WSFsmobs(n)%dlon   = grid_dd(n,5)
       WSFsmobs(n)%dlat   = grid_dd(n,6)
       WSFsmobs(n)%nc = nint((grid_dd(n,4) - grid_dd(n,2))/grid_dd(n,5)) + 1
       WSFsmobs(n)%nr = nint((grid_dd(n,3) - grid_dd(n,1))/grid_dd(n,6)) + 1

       WSFsmobs(n)%gridDesci = 0
       WSFsmobs(n)%gridDesci(1) = 0
       WSFsmobs(n)%gridDesci(2) = WSFsmobs(n)%nc
       WSFsmobs(n)%gridDesci(3) = WSFsmobs(n)%nr
       WSFsmobs(n)%gridDesci(4) = grid_dd(n,1)
       WSFsmobs(n)%gridDesci(5) = grid_dd(n,2)
       WSFsmobs(n)%gridDesci(6) = 128
       WSFsmobs(n)%gridDesci(7) = grid_dd(n,3)
       WSFsmobs(n)%gridDesci(8) = grid_dd(n,4)
       WSFsmobs(n)%gridDesci(9) = grid_dd(n,5)   !dlon
       WSFsmobs(n)%gridDesci(10) = grid_dd(n,6)  !dlat
       WSFsmobs(n)%gridDesci(20) = 64

       write(LDT_logunit,'(a,i6,a,i6,a,2f12.6,a,2f10.6)') &
            ' [INFO] WSF SM grid: ', WSFsmobs(n)%nc, ' x ', WSFsmobs(n)%nr, &
            ', lower left ', grid_dd(n,1), grid_dd(n,2), &
            ', dx dy ', grid_dd(n,5), grid_dd(n,6)

       ! Bilinear if the LDT grid is finer than the WSF grid, else
       ! averaging (identity when the two grids are the same).
       if(LDT_isLDTatAfinerResolution(n,WSFsmobs(n)%dlat)) then

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
    deallocate(grid_dd)
  end subroutine WSFsmobsinit

end module WSFsm_obsMod
