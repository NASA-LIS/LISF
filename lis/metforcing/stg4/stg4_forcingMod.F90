!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module stg4_forcingMod

!BOP
! !MODULE: stg4_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  National Center for Environmental Prediction (NCEP) Stage IV  
!  (STAGE4) Doppler Radar+gage product.  The Stage IV is a national
!  level product, on an hourly interval, and supplements mainly
!  the which base forcing precipitation (e.g., NLDAS) is being used.
! 
!  The implementation in LIS has the derived data type {\tt stg4\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncol]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrow]
!    Number of rows (along the north south dimension) for the input data
!  \item[stg4dir]
!    Directory containing the input data
!  \item[stg4time]
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
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_STG4      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: stg4_struc

!EOP

  type, public ::  stg4_type_dec
     real               :: ts
     integer            :: ncol                 ! Number of cols
     integer            :: nrow                 ! Number of rows
     character(len=LIS_CONST_PATH_LEN) :: stg4dir ! STAGE IV Directory
     real*8             :: stg4time             ! Nearest hourly instance of incoming file
     real*8             :: griduptime1          ! Designated time of STAGEIV grid change
     logical            :: gridchange1          ! Flag for when grid change occurs
     integer            :: mi                   ! Number of points in the input grid

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

  end type stg4_type_dec

  type(stg4_type_dec), allocatable :: stg4_struc(:)       ! STAGE IV Main Pointer array (number of LIS_domains)

contains
  
!BOP
!
! !ROUTINE: init_STG4
! \label{init_STG4}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 25May2006: Kristi Arsenault; Data and code implementation
! 
! !INTERFACE:
  subroutine init_STG4(findex)

! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
   use LIS_logMod,     only : LIS_logunit, LIS_endrun
   use LIS_FORC_AttributesMod

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
!   \item[readcrd\_stg4](\ref{readcrd_stg4}) \newline
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
       write(LIS_logunit,*) '[ERR] Currently the STAGE-IV precip forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    ! Stage IV Data structure -- Allocate for different LIS_domains
    allocate ( stg4_struc(LIS_rc%nnest) )

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

    ! Retrieve Stage IV Forcing Dataset Directory Name Location from lis.config
    call readcrd_stg4()

    do n=1, LIS_rc%nnest
       stg4_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, stg4_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 2

    gridDesci = 0

    do n=1,LIS_rc%nnest

       allocate(stg4_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(stg4_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       stg4_struc(n)%metdata1 = 0
       stg4_struc(n)%metdata2 = 0

   !-- Original STAGE IV starting time
       yr1 = 2002      
       mo1 = 01
       da1 = 01
       hr1 = 00; mn1 = 0; ss1 = 0

       stg4_struc(n)%ncol = 1121
       stg4_struc(n)%nrow = 881

       gridDesci(n,1) = 5.0                 ! Projection type (UPS)
       gridDesci(n,2) = stg4_struc(n)%ncol  ! X-dir amount of points
       gridDesci(n,3) = stg4_struc(n)%nrow  ! y-dir amount of points
       gridDesci(n,4) = 23.117000           ! Starting latitude point
       gridDesci(n,5) = -119.017000         ! Starting longitude point
       gridDesci(n,6) = 8.0         
       gridDesci(n,7) = 0                   ! Orientation
       gridDesci(n,8) = 4.7625              ! X-spacing length (kms) 
       gridDesci(n,9) = 4.7625              ! Y-spacing length (kms)
       gridDesci(n,10) = 60.0               ! True lat  (???)
       gridDesci(n,11) = -105.0             ! Standard longitude (???)
       gridDesci(n,20) = 64.0

       stg4_struc(n)%mi = stg4_struc(n)%ncol * stg4_struc(n)%nrow


! === BILINEAR INTERPOLATION ==== 
       if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
          allocate(stg4_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(stg4_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(stg4_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(stg4_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(stg4_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(stg4_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(stg4_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(stg4_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n, gridDesci(n,:),&
               stg4_struc(n)%n111, stg4_struc(n)%n121, &
               stg4_struc(n)%n211, stg4_struc(n)%n221, &
               stg4_struc(n)%w111, stg4_struc(n)%w121, &
               stg4_struc(n)%w211, stg4_struc(n)%w221 )

       elseif ( trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear" ) then 
! === BUDGET BILINEAR INTERPOLATION ==== 
          allocate(stg4_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(stg4_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(stg4_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(stg4_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(stg4_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(stg4_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(stg4_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(stg4_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n, gridDesci(n,:), &
                 stg4_struc(n)%n112, stg4_struc(n)%n122, &
                 stg4_struc(n)%n212, stg4_struc(n)%n222, &
                 stg4_struc(n)%w112, stg4_struc(n)%w122, &
                 stg4_struc(n)%w212, stg4_struc(n)%w222 )

       end if   ! End interp option statement

   !-- 1st CHANGE IN STAGE IV Grid Domain
       yr1 = 2004      !grid update time
       mo1 = 5 
       da1 = 10 
       hr1 = 19; mn1 = 0; ss1 = 0

       call LIS_date2time( stg4_struc(n)%griduptime1, updoy, upgmt,&
            yr1, mo1, da1, hr1, mn1, ss1 )
       stg4_struc(n)%gridchange1 = .true.

    enddo

  end subroutine init_STG4

end module stg4_forcingMod
