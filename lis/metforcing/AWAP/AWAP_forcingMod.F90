!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module AWAP_forcingMod

!BOP
! !MODULE: AWAP_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the AWAP precipitation data. 
! 
!  The implementation in LIS has the derived data type {\tt AWAP\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncol]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrow]
!    Number of rows (along the north south dimension) for the input data
!  \item[AWAPdir]
!    Directory containing the input data
!  \item[AWAPtime]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[mi]
!    Number of points in the input grid
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
! !REVISION HISTORY: 
! 30 Jan 2017: Sujay Kumar, Initial version
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_AWAP      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: AWAP_struc

!EOP

  type, public ::  AWAP_type_dec
     real               :: ts
     integer            :: ncol                 ! Number of cols
     integer            :: nrow                 ! Number of rows
     character(len=LIS_CONST_PATH_LEN) :: AWAPdir ! STAGE IV Directory
     real*8             :: AWAPtime             ! Nearest hourly instance of incoming file
     integer            :: mi                   ! Number of points in the input grid
     logical            :: interp_flag
! == Arrays for Bilinear Interpolation option (=1)
     integer, allocatable   :: n111(:), n121(:)
     integer, allocatable   :: n211(:), n221(:)
     real, allocatable      :: w111(:), w121(:)
     real, allocatable      :: w211(:), w221(:)
! == Arrays for Budget Bilinear Interpolation option (=2)
     integer, allocatable   :: n112(:,:), n122(:,:) ! Spatial Interpolation weights
     integer, allocatable   :: n212(:,:), n222(:,:) ! " "
     real,    allocatable   :: w112(:,:), w122(:,:) ! " "
     real,    allocatable   :: w212(:,:), w222(:,:) ! " "

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type AWAP_type_dec

  type(AWAP_type_dec), allocatable :: AWAP_struc(:)       ! STAGE IV Main Pointer array (number of LIS_domains)

contains
  
!BOP
!
! !ROUTINE: init_AWAP
! \label{init_AWAP}
!
! 
! !INTERFACE:
  subroutine init_AWAP(findex)

! !USES: 
   use LIS_coreMod
   use LIS_timeMgrMod
   use LIS_logMod, only : LIS_logunit, LIS_endrun

   implicit none
   integer,   intent(in) :: findex

! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for STAGE4
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_AWAP](\ref{readcrd_AWAP}) \newline
!     reads the runtime options specified for STAGE4 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!     computes the neighbor, weights for bilinear interpolation
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format - time of grid change 
!  \end{description}
!
!EOP
    real :: gridDesci(LIS_rc%nnest, 50)
    integer :: n

    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt


!=============================================================================

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the AWAP forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate ( AWAP_struc(LIS_rc%nnest) )

    call readcrd_AWAP()

    do n=1, LIS_rc%nnest
       AWAP_struc(n)%ts = 86400
       call LIS_update_timestep(LIS_rc, n, AWAP_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 1

    gridDesci = 0

    do n=1,LIS_rc%nnest

       allocate(AWAP_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(AWAP_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       AWAP_struc(n)%ncol = 886
       AWAP_struc(n)%nrow = 691

       gridDesci(n,1) = 0.0                 ! Projection type (UPS)
       gridDesci(n,2) = AWAP_struc(n)%ncol  ! X-dir amount of points
       gridDesci(n,3) = AWAP_struc(n)%nrow  ! y-dir amount of points
       gridDesci(n,4) = -44.50
       gridDesci(n,5) = 112.00
       gridDesci(n,6) = 128
       gridDesci(n,7) = -10.00
       gridDesci(n,8) = 156.25
       gridDesci(n,9) = 0.05
       gridDesci(n,10) = 0.05
       gridDesci(n,20) = 64.0

       AWAP_struc(n)%mi = AWAP_struc(n)%ncol * AWAP_struc(n)%nrow

       if ( LIS_isatAfinerResolution(n,gridDesci(n,9)) ) then
          AWAP_struc(n)%interp_flag = .true. 
! === BILINEAR INTERPOLATION ==== 
          if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
             allocate(AWAP_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWAP_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWAP_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWAP_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWAP_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWAP_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWAP_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWAP_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             
             call bilinear_interp_input(n, gridDesci(n,:),&
                  AWAP_struc(n)%n111, AWAP_struc(n)%n121, &
                  AWAP_struc(n)%n211, AWAP_struc(n)%n221, &
                  AWAP_struc(n)%w111, AWAP_struc(n)%w121, &
                  AWAP_struc(n)%w211, AWAP_struc(n)%w221 )
             
          elseif ( trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear" ) then 
             ! === BUDGET BILINEAR INTERPOLATION ==== 
             allocate(AWAP_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWAP_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWAP_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWAP_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWAP_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWAP_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWAP_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWAP_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             
             call conserv_interp_input(n, gridDesci(n,:), &
                  AWAP_struc(n)%n112, AWAP_struc(n)%n122, &
                  AWAP_struc(n)%n212, AWAP_struc(n)%n222, &
                  AWAP_struc(n)%w112, AWAP_struc(n)%w122, &
                  AWAP_struc(n)%w212, AWAP_struc(n)%w222 )
             
          end if   ! End interp option statement
       else

          AWAP_struc(n)%interp_flag = .false.
          allocate(AWAP_struc(n)%n111(AWAP_struc(n)%mi))

          call upscaleByAveraging_input(&
               gridDesci(n,:),              &
               LIS_rc%gridDesc(n,:),        &
               AWAP_struc(n)%mi,            &
               LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               AWAP_struc(n)%n111)

       endif
    enddo

  end subroutine init_AWAP

end module AWAP_forcingMod
