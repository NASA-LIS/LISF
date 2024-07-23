!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module nldas20_forcingMod
!BOP
! !MODULE: nldas20_forcingMod
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from nldas2_forcingMod.F90)
!
! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_NLDAS20    !defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: nldas20_struc
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data used in the North American
!  Land Data Assimilation System Phase 2.  The variables are produced
!  at 0.125 degree spatial resolution, and at hourly intervals.
!  For more details please view the forcing files manual available
!  at the following URL:
!    https://ldas.gsfc.nasa.gov/nldas/v2/forcing
!
!  The implementation in LIS has the derived data type {\tt nldas20\_struc}
!  that includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nldas20time1]
!    The nearest, previous hourly instance of the incoming
!    data (as a real time).
!  \item[nldas20time2]
!    The nearest, next hourly instance of the incoming
!    data (as a real time).
!  \item[nldas20foradir]
!    Directory containing the FORA netCDF-4 input data
!  \item[nldas20foradir]
!    Directory containing the FORB netCDF-4 input data
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
  type, public :: nldas20_type_dec
     real         :: ts
     integer      :: ncold, nrold ! AWIPS 212 dimensions
     character(len=LIS_CONST_PATH_LEN) :: nldas20foradir,nldas20forbdir
     real*8       :: nldas20time1,nldas20time2
     integer      :: model_level_data
     integer      :: model_level_press
     integer      :: model_pcp_data
     integer      :: model_dswrf_data

     real, allocatable      :: orig_ediff(:)
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

     integer           :: findtime1, findtime2
     integer           :: nIter, st_iterid,en_iterid ! Forecast parameters

     real, allocatable :: metdata1(:,:,:)
     real, allocatable :: metdata2(:,:,:)

  end type nldas20_type_dec

  type(nldas20_type_dec), allocatable :: nldas20_struc(:)

contains
!BOP
!
! !ROUTINE: init_NLDAS20
! \label{init_NLDAS20}
!
! !INTERFACE:
  subroutine init_NLDAS20(findex)
! !USES:
    use LIS_coreMod,    only            : LIS_rc
    use LIS_timeMgrMod, only            : LIS_update_timestep
    use LIS_logMod,     only            : LIS_logunit,LIS_endrun
    use map_utils,      only            : proj_latlon
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use LIS_forecastMod

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for NLDAS-2
!  data. The grid description arrays are based on the netCDF-4
!  data and followed in the LIS interpolation schemes
!  (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_nldas20](\ref{readcrd_nldas20}) \newline
!     reads the runtime options specified for NLDAS-2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[read\_nldas20\_elev](\ref{read_nldas20_elev}) \newline
!    reads the native elevation of the NLDAS-2 data to be used
!    for topographic adjustments to the forcing
!  \end{description}
!
!EOP
    integer :: n

    external :: readcrd_nldas20
    external :: bilinear_interp_input
    external :: conserv_interp_input
    external :: neighbor_interp_input
    external :: read_orig_nldas20_elevdiff
    external :: read_nldas20_elev


    allocate(nldas20_struc(LIS_rc%nnest))
    call readcrd_nldas20()

    do n = 1,LIS_rc%nnest
       nldas20_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc,n,nldas20_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 13

! Set NLDAS-2 grid dimensions and extent information:
    nldas20_struc(:)%ncold = 464
    nldas20_struc(:)%nrold = 224

    do n = 1,LIS_rc%nnest

! Forecast mode:
       if (LIS_rc%forecastMode.eq.1) then
          if (mod(LIS_rc%nensem(n),                                  &
               LIS_forecast_struc(1)%niterations).ne.0) then
             write(LIS_logunit,*)                                    &
                  "[ERR] The number of ensembles must be a multiple"
             write(LIS_logunit,*)                                    &
                  "[ERR] of the number of iterations"
             write(LIS_logunit,*)                                    &
                  "[ERR] nensem = ",LIS_rc%nensem(n)
             write(LIS_logunit,*)                                    &
                  "[ERR] niter = ",LIS_forecast_struc(1)%niterations
             call LIS_endrun()
          endif

          allocate(nldas20_struc(n)%metdata1(                        &
               LIS_forecast_struc(1)%niterations,                    &
               LIS_rc%met_nf(findex),LIS_rc%ngrid(n)))
          allocate(nldas20_struc(n)%metdata2(                        &
               LIS_forecast_struc(1)%niterations,                    &
               LIS_rc%met_nf(findex),LIS_rc%ngrid(n)))

          nldas20_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          nldas20_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          nldas20_struc(n)%nIter = LIS_forecast_struc(1)%niterations

! Regular retrospective or non-forecast mode:
       else
          allocate(nldas20_struc(n)%metdata1(1,LIS_rc%met_nf(findex), &
               LIS_rc%ngrid(n)))
          allocate(nldas20_struc(n)%metdata2(1,LIS_rc%met_nf(findex), &
               LIS_rc%ngrid(n)))

          nldas20_struc(n)%st_iterid = 1
          nldas20_struc(n)%en_iterId = 1
          nldas20_struc(n)%nIter = 1
       endif

       nldas20_struc(n)%metdata1 = 0
       nldas20_struc(n)%metdata2 = 0
       nldas20_struc(n)%gridDesc = 0
       nldas20_struc(n)%findtime1 = 0
       nldas20_struc(n)%findtime2 = 0
       nldas20_struc(n)%gridDesc(1) = 0
       nldas20_struc(n)%gridDesc(2) = nldas20_struc(n)%ncold
       nldas20_struc(n)%gridDesc(3) = nldas20_struc(n)%nrold
       nldas20_struc(n)%gridDesc(4) = 25.0625
       nldas20_struc(n)%gridDesc(5) = -124.9375
       nldas20_struc(n)%gridDesc(6) = 128
       nldas20_struc(n)%gridDesc(7) = 52.9375
       nldas20_struc(n)%gridDesc(8) = -67.0625
       nldas20_struc(n)%gridDesc(9) = 0.125
       nldas20_struc(n)%gridDesc(10) = 0.125
       nldas20_struc(n)%gridDesc(20) = 64

! Check for grid and interp option selected:
       if ((nldas20_struc(n)%gridDesc(9).eq.LIS_rc%gridDesc(n,9)).and. &
            (nldas20_struc(n)%gridDesc(10).eq.LIS_rc%gridDesc(n,10)).and. &
            (LIS_rc%gridDesc(n,1).eq.proj_latlon).and.               &
            (LIS_rc%met_interp(findex).ne."neighbor")) then
          write(LIS_logunit,*)                                       &
               "[ERR] The NLDAS grid was selected for the LIS run domain;"
          write(LIS_logunit,*)                                       &
               "[ERR] however, 'bilinear', 'budget-bilinear', or some"
          write(LIS_logunit,*)                                       &
               "[ERR] other unknown option was selected to spatially"
          write(LIS_logunit,*)                                       &
               "[ERR] downscale the grid, which will cause errors"
          write(LIS_logunit,*)                                       &
               "[ERR] during runtime.  Please select 'neighbor'."
          write(LIS_logunit,*) "[ERR] Program stopping ..."
          call LIS_endrun()
       endif

       nldas20_struc(n)%mi = nldas20_struc(n)%ncold*nldas20_struc(n)%nrold

! Setting up weights for spatial interpolation:
       select case(LIS_rc%met_interp(findex))

       case ("bilinear")
          allocate(nldas20_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,nldas20_struc(n)%gridDesc(:), &
               nldas20_struc(n)%n111,nldas20_struc(n)%n121,          &
               nldas20_struc(n)%n211,nldas20_struc(n)%n221,          &
               nldas20_struc(n)%w111,nldas20_struc(n)%w121,          &
               nldas20_struc(n)%w211,nldas20_struc(n)%w221)

       case ("budget-bilinear")
          allocate(nldas20_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(nldas20_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,nldas20_struc(n)%gridDesc(:), &
               nldas20_struc(n)%n111,nldas20_struc(n)%n121,          &
               nldas20_struc(n)%n211,nldas20_struc(n)%n221,          &
               nldas20_struc(n)%w111,nldas20_struc(n)%w121,          &
               nldas20_struc(n)%w211,nldas20_struc(n)%w221)

          allocate(nldas20_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas20_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas20_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas20_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas20_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas20_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas20_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(nldas20_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,nldas20_struc(n)%gridDesc(:),  &
               nldas20_struc(n)%n112,nldas20_struc(n)%n122,          &
               nldas20_struc(n)%n212,nldas20_struc(n)%n222,          &
               nldas20_struc(n)%w112,nldas20_struc(n)%w122,          &
               nldas20_struc(n)%w212,nldas20_struc(n)%w222)

       case ("neighbor")
          allocate(nldas20_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n,nldas20_struc(n)%gridDesc(:), &
               nldas20_struc(n)%n113)

       case default
          write(LIS_logunit,*)                                       &
               "[ERR] Interpolation option not specified for NLDAS-2"
          write(LIS_logunit,*) "[ERR] Program stopping ..."
          call LIS_endrun()
       end select

! Read in elevation difference and NLDAS-2 elevation maps:
       if (LIS_rc%met_ecor(findex).ne."none") then
          allocate(nldas20_struc(n)%orig_ediff(                      &
               nldas20_struc(n)%ncold*nldas20_struc(n)%nrold))
          call read_orig_nldas20_elevdiff(n)
          call read_nldas20_elev(n,findex)
       endif

! Set up precipitation climate downscaling:
       if (LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,                    &
               nint(nldas20_struc(n)%gridDesc(2)),                   &
               nint(nldas20_struc(n)%gridDesc(3)))
       endif
    enddo

  end subroutine init_NLDAS20

end module nldas20_forcingMod

