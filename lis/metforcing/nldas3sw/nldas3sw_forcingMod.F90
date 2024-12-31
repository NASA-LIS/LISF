!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nldas3sw_forcingMod
!BOP
! !MODULE: nldas3sw_forcingMod
!
! !REVISION HISTORY:
! 27 Dec 2024: David Mocko, Initial Specification
!                           (derived from nldas20_forcingMod.F90)
!
! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_nldas3sw    !defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nldas3sw_struc
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used for
!  the reading of the CERES hourly 4-km SWdown data, for the downscaling
!  and slope-aspect correction down to the NLDAS-3 1-km grid.
!
!  The implementation in LIS has the derived data type {\tt nldas3sw\_struc}
!  that includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nldas3swtime1]
!    The nearest, previous hourly instance of the incoming
!    data (as a real time).
!  \item[nldas3swtime2]
!    The nearest, next hourly instance of the incoming
!    data (as a real time).
!  \item[nldas3swfordir]
!    Directory containing the CERES 4-km binary input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for nearest neighbor interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \end{description}
!
!EOP
  type, public :: nldas3sw_type_dec
     real         :: ts
     integer      :: ncold, nrold
     real*8       :: nldas3swtime1,nldas3swtime2
     character(len=LIS_CONST_PATH_LEN) :: nldas3swfordir

     real                   :: gridDesc(50)
     integer                :: mi
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)
     integer, allocatable   :: n113(:)
     integer                :: findtime1,findtime2
     real, allocatable      :: metdata1(:,:)
     real, allocatable      :: metdata2(:,:)

  end type nldas3sw_type_dec

  type(nldas3sw_type_dec), allocatable :: nldas3sw_struc(:)

contains
!BOP
!
! !ROUTINE: init_nldas3sw
! \label{init_nldas3sw}
!
! !INTERFACE:
  subroutine init_nldas3sw(findex)
! !USES:
    use LIS_coreMod,    only            : LIS_rc
    use LIS_timeMgrMod, only            : LIS_update_timestep
    use LIS_logMod,     only            : LIS_logunit,LIS_endrun
    use map_utils,      only            : proj_latlon

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for CERES
!  SWdown data.  The grid description arrays are based on the
!  binary data and followed in the LIS interpolation schemes
!  (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_nldas3sw](\ref{readcrd_nldas3sw}) \newline
!     reads the runtime options specified for NLDAS-3 SWdown data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!
!EOP
    integer :: n

    external :: readcrd_nldas3sw
    external :: bilinear_interp_input
    external :: conserv_interp_input
    external :: neighbor_interp_input

    allocate(nldas3sw_struc(LIS_rc%nnest))
    call readcrd_nldas3sw()

    do n = 1,LIS_rc%nnest
       nldas3sw_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc,n,nldas3sw_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 1

! Set CERES 4-km grid dimensions and extent information:
    nldas3sw_struc(:)%ncold = 2925
    nldas3sw_struc(:)%nrold = 1625

    do n = 1,LIS_rc%nnest

! Regular retrospective or non-forecast mode:
       allocate(nldas3sw_struc(n)%metdata1(LIS_rc%met_nf(findex),      &
                LIS_rc%ngrid(n)))
       allocate(nldas3sw_struc(n)%metdata2(LIS_rc%met_nf(findex),      &
                LIS_rc%ngrid(n)))

       nldas3sw_struc(n)%metdata1 = 0
       nldas3sw_struc(n)%metdata2 = 0
       nldas3sw_struc(n)%gridDesc = 0
       nldas3sw_struc(n)%findtime1 = 0
       nldas3sw_struc(n)%findtime2 = 0
       nldas3sw_struc(n)%gridDesc(1) = 0
       nldas3sw_struc(n)%gridDesc(2) = nldas3sw_struc(n)%ncold
       nldas3sw_struc(n)%gridDesc(3) = nldas3sw_struc(n)%nrold
       nldas3sw_struc(n)%gridDesc(4) = 7.02
       nldas3sw_struc(n)%gridDesc(5) = -168.98
       nldas3sw_struc(n)%gridDesc(6) = 128
       nldas3sw_struc(n)%gridDesc(7) = 71.98
       nldas3sw_struc(n)%gridDesc(8) = -52.02
       nldas3sw_struc(n)%gridDesc(9) = 0.04
       nldas3sw_struc(n)%gridDesc(10) = 0.04
       nldas3sw_struc(n)%gridDesc(20) = 64

! Check for grid and interp option selected:
       if ((nldas3sw_struc(n)%gridDesc(9).eq.LIS_rc%gridDesc(n,9)).and.&
         (nldas3sw_struc(n)%gridDesc(10).eq.LIS_rc%gridDesc(n,10)).and.&
            (LIS_rc%gridDesc(n,1).eq.proj_latlon).and.                 &
            (LIS_rc%met_interp(findex).ne."neighbor")) then
          write(LIS_logunit,*)                                         &
               "[ERR] The NLDAS grid was selected for the LIS domain;"
          write(LIS_logunit,*)                                         &
               "[ERR] however, 'bilinear', 'budget-bilinear', or some"
          write(LIS_logunit,*)                                         &
               "[ERR] other unknown option was selected to spatially"
          write(LIS_logunit,*)                                         &
               "[ERR] downscale the grid, which will cause errors"
          write(LIS_logunit,*)                                         &
               "[ERR] during runtime.  Please select 'neighbor'."
          write(LIS_logunit,*) "[ERR] Program stopping ..."
          call LIS_endrun()
       endif

       nldas3sw_struc(n)%mi = nldas3sw_struc(n)%ncold*nldas3sw_struc(n)%nrold

! Setting up weights for spatial interpolation:
       select case(LIS_rc%met_interp(findex))

       case ("bilinear")
          allocate(nldas3sw_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,nldas3sw_struc(n)%gridDesc(:),  &
               nldas3sw_struc(n)%n111,nldas3sw_struc(n)%n121,          &
               nldas3sw_struc(n)%n211,nldas3sw_struc(n)%n221,          &
               nldas3sw_struc(n)%w111,nldas3sw_struc(n)%w121,          &
               nldas3sw_struc(n)%w211,nldas3sw_struc(n)%w221)

       case ("budget-bilinear")
          allocate(nldas3sw_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas3sw_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,nldas3sw_struc(n)%gridDesc(:),  &
               nldas3sw_struc(n)%n111,nldas3sw_struc(n)%n121,          &
               nldas3sw_struc(n)%n211,nldas3sw_struc(n)%n221,          &
               nldas3sw_struc(n)%w111,nldas3sw_struc(n)%w121,          &
               nldas3sw_struc(n)%w211,nldas3sw_struc(n)%w221)

          allocate(nldas3sw_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas3sw_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas3sw_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas3sw_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas3sw_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas3sw_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas3sw_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas3sw_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,nldas3sw_struc(n)%gridDesc(:),   &
               nldas3sw_struc(n)%n112,nldas3sw_struc(n)%n122,          &
               nldas3sw_struc(n)%n212,nldas3sw_struc(n)%n222,          &
               nldas3sw_struc(n)%w112,nldas3sw_struc(n)%w122,          &
               nldas3sw_struc(n)%w212,nldas3sw_struc(n)%w222)

       case ("neighbor")
          allocate(nldas3sw_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n,nldas3sw_struc(n)%gridDesc(:),  &
               nldas3sw_struc(n)%n113)

       case default
          write(LIS_logunit,*)                                         &
               "[ERR] Interpolation option not specified for NLDAS-3 SW"
          write(LIS_logunit,*) "[ERR] Program stopping ..."
          call LIS_endrun()
       end select

    enddo

  end subroutine init_nldas3sw

end module nldas3sw_forcingMod

