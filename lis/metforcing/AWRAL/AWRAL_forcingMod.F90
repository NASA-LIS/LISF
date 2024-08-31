!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module AWRAL_forcingMod

!BOP
! !MODULE: AWRAL_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the AWRAL precipitation data. 
! 
!  The implementation in LIS has the derived data type {\tt AWRAL\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncol]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrow]
!    Number of rows (along the north south dimension) for the input data
!  \item[AWRALdir]
!    Directory containing the input data
!  \item[AWRALtime]
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
  public :: init_AWRAL      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: AWRAL_struc

!EOP

  type, public ::  AWRAL_type_dec
     real               :: ts
     integer            :: ncol                 ! Number of cols
     integer            :: nrow                 ! Number of rows
     real 		:: gridDesci(50)
     character(len=LIS_CONST_PATH_LEN) :: AWRALdir ! STAGE IV Directory
     real*8             :: AWRALtime             ! Nearest daily instance of incoming file
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

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type AWRAL_type_dec

  type(AWRAL_type_dec), allocatable :: AWRAL_struc(:)       ! STAGE IV Main Pointer array (number of LIS_domains)

contains
  
!BOP
!
! !ROUTINE: init_AWRAL
! \label{init_AWRAL}
!
! 
! !INTERFACE:
  subroutine init_AWRAL(findex)

! !USES: 
   use LIS_coreMod
   use LIS_timeMgrMod
   use LIS_logMod, only : LIS_logunit, LIS_endrun
   use LIS_constantsMod, only : LIS_CONST_CDAY

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
!   \item[readcrd\_AWRAL](\ref{readcrd_AWRAL}) \newline
!     reads the runtime options specified for STAGE4 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!     computes the neighbor, weights for bilinear interpolation
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format - time of grid change 
!  \end{description}
!
!EOP
    integer :: n

    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt


!=============================================================================

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the AWRAL forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate ( AWRAL_struc(LIS_rc%nnest) )

    call readcrd_AWRAL()
	
    do n=1, LIS_rc%nnest
       AWRAL_struc(n)%ts = LIS_CONST_CDAY
       call LIS_update_timestep(LIS_rc, n, AWRAL_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 6

    do n=1,LIS_rc%nnest

       allocate(AWRAL_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(AWRAL_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
   

       AWRAL_struc(n)%ncol = 841
       AWRAL_struc(n)%nrow = 681
       AWRAL_struc(n)%gridDesci = 0

       ! Define parameters for lat/lon projection
       AWRAL_struc(n)%gridDesci(1) = 0      ! indicates lat/lon projection
       AWRAL_struc(n)%gridDesci(2) = AWRAL_struc(n)%ncol  ! number of columns in the domain
       AWRAL_struc(n)%gridDesci(3) = AWRAL_struc(n)%nrow  ! number of rows in the domain
       AWRAL_struc(n)%gridDesci(4) = -44.00 ! latitude of the lower left corner CELL CENTER of the domain
       AWRAL_struc(n)%gridDesci(5) = 112.00 ! longitude of the lower left corner CELL CENTER of the domain
       AWRAL_struc(n)%gridDesci(7) = -10.00 ! latitude of the upper right corner CELL CENTER of the domain
       AWRAL_struc(n)%gridDesci(8) = 154.00 ! longitude of the upper right corner CELL CENTER of the domain
       AWRAL_struc(n)%gridDesci(9) = 0.05   ! spatial resolution (in degrees) along the E-W dimension
       AWRAL_struc(n)%gridDesci(10) = 0.05  ! spatial resolution (in degrees) along the N-S dimension
       AWRAL_struc(n)%gridDesci(20) = 64.0  ! used to specify the ordering of data
       ! (Non-divisible by 32 indicates E-W ordering, else N-S ordering)
       !  Set now to 255 per interp/get_fieldpos explanation:
       !   E-W ordering indicates elements located in row-order array, i.e., one row, then next and so on

       AWRAL_struc(n)%mi = AWRAL_struc(n)%ncol * AWRAL_struc(n)%nrow


       if ( LIS_isatAfinerResolution(n, AWRAL_struc(n)%gridDesci(9)) ) then
          AWRAL_struc(n)%interp_flag = .true. 
! === BILINEAR INTERPOLATION ==== 
          if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
             allocate(AWRAL_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWRAL_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWRAL_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWRAL_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWRAL_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWRAL_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWRAL_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(AWRAL_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             
             call bilinear_interp_input(n, AWRAL_struc(n)%gridDesci(:),&
                  AWRAL_struc(n)%n111, AWRAL_struc(n)%n121, &
                  AWRAL_struc(n)%n211, AWRAL_struc(n)%n221, &
                  AWRAL_struc(n)%w111, AWRAL_struc(n)%w121, &
                  AWRAL_struc(n)%w211, AWRAL_struc(n)%w221 )
             
          elseif ( trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear" ) then 
             ! === BUDGET BILINEAR INTERPOLATION ==== 
             allocate(AWRAL_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWRAL_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWRAL_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWRAL_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWRAL_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWRAL_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWRAL_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             allocate(AWRAL_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
             
             call conserv_interp_input(n, AWRAL_struc(n)%gridDesci(:), &
                  AWRAL_struc(n)%n112, AWRAL_struc(n)%n122, &
                  AWRAL_struc(n)%n212, AWRAL_struc(n)%n222, &
                  AWRAL_struc(n)%w112, AWRAL_struc(n)%w122, &
                  AWRAL_struc(n)%w212, AWRAL_struc(n)%w222 )
             
          end if   ! End interp option statement
       else

          AWRAL_struc(n)%interp_flag = .false.
          allocate(AWRAL_struc(n)%n111(AWRAL_struc(n)%mi))

          call upscaleByAveraging_input(&
               AWRAL_struc(n)%gridDesci(:),              &
               LIS_rc%gridDesc(n,:),        &
               AWRAL_struc(n)%mi,            &
               LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               AWRAL_struc(n)%n111)

       endif
    enddo

  end subroutine init_AWRAL

end module AWRAL_forcingMod
