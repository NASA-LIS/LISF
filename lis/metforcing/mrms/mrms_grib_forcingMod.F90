!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module mrms_grib_forcingMod

!BOP
! !MODULE: mrms_grib_forcingMod
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  National Center for Environmental Prediction (NCEP) MRMS  
!  (MRMS) Doppler Radar+gage product.  The MRMS is a national
!  level product, on an hourly interval, and supplements mainly
!  the which base forcing precipitation (e.g., NLDAS) is being used.
! 
!  The implementation in LIS has the derived data type {\tt mrms\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncol]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrow]
!    Number of rows (along the north south dimension) for the input data
!  \item[mrmsdir]
!    Directory containing the input data
!  \item[mrmstime]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[mi]
!    Number of points in the input grid
!  \item[rlat1]
!    Array containing the latitudes of the input grid for each corresponding
!    grid point in LIS (to be used for bilinear interpolation)
!  \item[rlon1]
!    Array containing the longitudes of the input grid for each corresponding
!    grid point in LIS (to be used for bilinear interpolation)
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES: 
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_MRMS_grib      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: mrms_grib_struc

!EOP

  type, public ::  mrms_grib_type_dec
     real               :: ts
     integer            :: ncol                 ! Number of cols
     integer            :: nrow                 ! Number of rows
     character*80       :: mrms_grib_dir        ! MRMS Directory
     real*8             :: mrms_grib_time       ! Nearest hourly instance of incoming file
     integer            :: mrms_mask_opt        ! Flag for whether or not to use mask 1=Yes
     real*8             :: mrms_mask_thresh     ! Threshold for masking MRMS data
     character*150      :: mrms_mask_dir        ! Directory of MRMS masks
     integer            :: mi                   ! Number of points in the input grid

! == Arrays for Bilinear Interpolation option (=1)
     real, allocatable      :: rlat1(:)
     real, allocatable      :: rlon1(:)
     integer, allocatable   :: n111(:), n121(:)
     integer, allocatable   :: n211(:), n221(:)
     real, allocatable      :: w111(:), w121(:)
     real, allocatable      :: w211(:), w221(:)
! == Arrays for Budget Bilinear Interpolation option (=2)
     real,    allocatable   :: rlat2(:)             ! Array containing lats of the input grid
     real,    allocatable   :: rlon2(:)             ! Array containing lons of the input grid
     integer, allocatable   :: n112(:,:), n122(:,:) ! Spatial Interpolation weights
     integer, allocatable   :: n212(:,:), n222(:,:) ! " "
     real,    allocatable   :: w112(:,:), w122(:,:) ! " "
     real,    allocatable   :: w212(:,:), w222(:,:) ! " "
! == Arrays for Nearest Neighbor Interpolation Options !JE 
     integer, allocatable   :: n113(:)
! J. Erlingis 

     real, allocatable :: metdata1(:,:)
     real, allocatable :: metdata2(:,:)

! J.Case (ENSCO, Inc. 2/13/2015)
! == Array for upscale by averaging.
     integer, allocatable       :: n11(:)   ! Array that maps the location of each 
                                        ! input grid point in the output grid.

  end type mrms_grib_type_dec

  type(mrms_grib_type_dec), allocatable :: mrms_grib_struc(:) ! MRMS Main Pointer array (number of LIS_domains)

contains
  
!BOP
!
! !ROUTINE: init_MRMS_grib
! \label{init_MRMS_grib}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 25May2006: Kristi Arsenault; Data and code implementation
! 13 Feb 2015: Jonathan Case; Modified for MRMS QPE
! 07 Sep 2017: Jessica Erlingis; Modified for operational MRMS QPE
! 
! !INTERFACE:
  subroutine init_MRMS_grib(findex)

! !USES: 
   use LIS_coreMod, only: LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
   use LIS_logMod, only: LIS_logunit, LIS_endrun
   use LIS_FORC_AttributesMod

   implicit none
   integer,   intent(in) :: findex

! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for MRMS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_mrms_grib](\ref{readcrd_mrms_grib}) \\
!     reads the runtime options specified for operational MRMS data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \\
!     computes the neighbor, weights for bilinear interpolation
!   \item[LIS\_date2time](\ref{LIS_date2time}) \\
!     converts date to the real time format - time of grid change 
!  \end{description}
!
!EOP
    real :: gridDesci(LIS_rc%nnest, 50)
    integer :: n

    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt

! J.Case -- local variables
    character*50 :: dom
    real         :: res

!=============================================================================

! - MRMS Data structure -- Allocate for different LIS_domains
    allocate ( mrms_grib_struc(LIS_rc%nnest) )

! - Retrieve MRMS Forcing Dataset Directory Name Location from lis.config
    call readcrd_mrms_grib()

    do n=1, LIS_rc%nnest
       mrms_grib_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, mrms_grib_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 2

    gridDesci = 0

    do n=1,LIS_rc%nnest

    allocate(mrms_grib_struc(n)%metdata1(LIS_rc%met_nf(findex), &
        LIS_rc%ngrid(n)))
    allocate(mrms_grib_struc(n)%metdata2(LIS_rc%met_nf(findex), &
        LIS_rc%ngrid(n)))

    mrms_grib_struc(n)%metdata1 = 0
    mrms_grib_struc(n)%metdata2 = 0

   !-- Original MRMS Grid Domain and starting time
       mrms_grib_struc(n)%ncol = 7000
       mrms_grib_struc(n)%nrow = 3500
       yr1 = 2002
       mo1 = 01
       da1 = 01
       hr1 = 00; mn1 = 0; ss1 = 0

       gridDesci(n,:) = 0.0
       gridDesci(n,1) = 0                       ! Projection type (Lat/Lon)
       gridDesci(n,2) = mrms_grib_struc(n)%ncol ! X-dir amount of points
       gridDesci(n,3) = mrms_grib_struc(n)%nrow ! y-dir amount of points
       gridDesci(n,4) =  54.995                 ! Starting latitude point
       gridDesci(n,5) = -129.995                ! Starting longitude point
       gridDesci(n,6) = 128                     ! (not used)
       gridDesci(n,7) = 20.005                  ! Ending latitude point
       gridDesci(n,8) =  -60.005                ! Ending longitude point 
       gridDesci(n,9) =  0.01                   ! spatial resolution in W-E dirn (deg)
       gridDesci(n,10) = 0.01                   ! spatial resolution in S-N dirn (deg)
       gridDesci(n,20) = 64                     ! N-S ordering (number divisible by 32; same as in NLDAS2)

       mrms_grib_struc(n)%mi = mrms_grib_struc(n)%ncol * mrms_grib_struc(n)%nrow

       dom=LIS_rc%lis_map_proj
       res=LIS_rc%gridDesc(n,9)

! J.Case (2/13/2015) -- Use interpolation only if LIS resolution is < MRMS resolution.
!                       Otherwise, use upscaling by default.
       if ( (((dom.eq."latlon").or.(dom.eq."gaussian")) .and. (res.eq.0.01)) .or. & !JE Add equal case
            (((dom.eq."mercator").or.(dom.eq."lambert").or.                       &
              (dom.eq."polar").or.(dom.eq."UTM")) .and. (res.eq.1.0)) ) then

          allocate(mrms_grib_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n, gridDesci(n,:), &
               mrms_grib_struc(n)%n113)
          write(LIS_logunit,*) 'MRMS and LIS grid spacing are equal. '// &
               'Using nearest neighbor interpolation.'

       elseif ( (((dom.eq."latlon").or.(dom.eq."gaussian")) .and. (res.lt.0.01)) .or. & !Je change to le (?)
            (((dom.eq."mercator").or.(dom.eq."lambert").or.                       &
              (dom.eq."polar").or.(dom.eq."UTM")) .and. (res.lt.1.0)) ) then
! === BILINEAR INTERPOLATION ==== 
          if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
             allocate(mrms_grib_struc(n)%rlat1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%rlon1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!JE              call bilinear_interp_input( gridDesci(n,:), LIS_rc%gridDesc(n,:), &
!JE                  LIS_rc%lnc(n)*LIS_rc%lnr(n), mrms_grib_struc(n)%rlat1, &
!JE                  mrms_grib_struc(n)%rlon1, &
!JE                  mrms_grib_struc(n)%n111, mrms_grib_struc(n)%n121, &
!JE                  mrms_grib_struc(n)%n211, mrms_grib_struc(n)%n221, &
!JE                  mrms_grib_struc(n)%w111, mrms_grib_struc(n)%w121, &
!JE                  mrms_grib_struc(n)%w211, mrms_grib_struc(n)%w221 )
          call bilinear_interp_input(n, gridDesci(n,:),&
               mrms_grib_struc(n)%n111,mrms_grib_struc(n)%n121,&
               mrms_grib_struc(n)%n211,mrms_grib_struc(n)%n221,&
               mrms_grib_struc(n)%w111,mrms_grib_struc(n)%w121,&
               mrms_grib_struc(n)%w211,mrms_grib_struc(n)%w221)

          elseif ( trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear" ) then 
! === BUDGET BILINEAR INTERPOLATION ==== 
             allocate(mrms_grib_struc(n)%rlat2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%rlon2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(mrms_grib_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(mrms_grib_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(mrms_grib_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(mrms_grib_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(mrms_grib_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(mrms_grib_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(mrms_grib_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(mrms_grib_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
!JE             call conserv_interp_input( gridDesci(n,:), LIS_rc%gridDesc(n,:), &
!JE                  LIS_rc%lnc(n)*LIS_rc%lnr(n), mrms_grib_struc(n)%rlat2, &
!JE                  mrms_grib_struc(n)%rlon2, &
!JE                  mrms_grib_struc(n)%n112, mrms_grib_struc(n)%n122, &
!JE                  mrms_grib_struc(n)%n212, mrms_grib_struc(n)%n222, &
!JE                  mrms_grib_struc(n)%w112, mrms_grib_struc(n)%w122, &
!JE                  mrms_grib_struc(n)%w212, mrms_grib_struc(n)%w222 )

          call conserv_interp_input(n, gridDesci(n,:),&
               mrms_grib_struc(n)%n112,mrms_grib_struc(n)%n122,&
               mrms_grib_struc(n)%n212,mrms_grib_struc(n)%n222,&
               mrms_grib_struc(n)%w112,mrms_grib_struc(n)%w122,&
               mrms_grib_struc(n)%w212,mrms_grib_struc(n)%w222)

          end if   ! End interp option statement

       else !! use upscaling if LIS resolution exceeds data resolution
! === UPSCALING ==== 
          allocate(mrms_grib_struc(n)%n11(mrms_grib_struc(n)%mi))
          call upscaleByAveraging_input( gridDesci(n,:), LIS_rc%gridDesc(n,:), &
               mrms_grib_struc(n)%ncol*mrms_grib_struc(n)%nrow, LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               mrms_grib_struc(n)%n11)

       endif !! dom/res check
! J.Case (end interpolation mods)

    enddo

  end subroutine init_MRMS_grib

end module mrms_grib_forcingMod
