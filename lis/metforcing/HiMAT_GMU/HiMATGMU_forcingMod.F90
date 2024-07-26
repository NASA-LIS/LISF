!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HiMATGMU_forcingMod

!BOP
! !MODULE: HiMATGMU_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the downscaled MERRA2 precipitation data 
!  over the High Mountain Asia domain from GMU. 
! 
!  The implementation in LIS has the derived data type {\tt HiMATGMU\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncol]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrow]
!    Number of rows (along the north south dimension) for the input data
!  \item[HiMATGMUdir]
!    Directory containing the input data
!  \item[HiMATGMUtime]
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
  public :: init_HIMATGMU      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: HiMATGMU_struc

!EOP

  type, public ::  HiMATGMU_type_dec
     real               :: ts
     integer            :: ncol                 ! Number of cols
     integer            :: nrow                 ! Number of rows
     character(len=LIS_CONST_PATH_LEN) :: HiMATGMUdir ! STAGE IV Directory
     real*8             :: HiMATGMUtime             ! Nearest hourly instance of incoming file
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

  end type HiMATGMU_type_dec

  type(HiMATGMU_type_dec), allocatable :: HiMATGMU_struc(:)       ! STAGE IV Main Pointer array (number of LIS_domains)

contains
  
!BOP
!
! !ROUTINE: init_HIMATGMU
! \label{init_HIMATGMU}
!
! !REVISION HISTORY: 
! 28 July 2017: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_HIMATGMU(findex)

! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
   use LIS_logMod,     only : LIS_logunit, LIS_endrun
   use LIS_FORC_AttributesMod

   implicit none
   integer,   intent(in) :: findex

! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for the HiMAT 
!  GMU data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_HiMATGMU](\ref{readcrd_HiMATGMU}) \newline
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
       write(LIS_logunit,*) '[ERR] Currently the HiMAT GMU precip forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    ! Stage IV Data structure -- Allocate for different LIS_domains
    allocate ( HiMATGMU_struc(LIS_rc%nnest) )

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

    ! Retrieve Stage IV Forcing Dataset Directory Name Location from lis.config
    call readcrd_HiMATGMU()

    do n=1, LIS_rc%nnest
       HiMATGMU_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, HiMATGMU_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 2

    gridDesci = 0

    do n=1,LIS_rc%nnest

       allocate(HiMATGMU_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(HiMATGMU_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       HiMATGMU_struc(n)%metdata1 = 0
       HiMATGMU_struc(n)%metdata2 = 0

   !-- Original STAGE IV starting time
       yr1 = 2002      
       mo1 = 01
       da1 = 01
       hr1 = 00; mn1 = 0; ss1 = 0

       HiMATGMU_struc(n)%ncol = 3000
       HiMATGMU_struc(n)%nrow = 1500

       gridDesci(n,1) = 0.0                 ! Projection type (UPS)
       gridDesci(n,2) = HiMATGMU_struc(n)%ncol  ! X-dir amount of points
       gridDesci(n,3) = HiMATGMU_struc(n)%nrow  ! y-dir amount of points
       gridDesci(n,4) = 25.0
       gridDesci(n,5) = 75.0
       gridDesci(n,6) = 128.0         
       gridDesci(n,7) = 40.0 
       gridDesci(n,8) = 105.0
       gridDesci(n,9) =  0.01
       gridDesci(n,10) = 0.01
       gridDesci(n,20) = 64.0

       HiMATGMU_struc(n)%mi = HiMATGMU_struc(n)%ncol * HiMATGMU_struc(n)%nrow


! === BILINEAR INTERPOLATION ==== 
       if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
          allocate(HiMATGMU_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(HiMATGMU_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(HiMATGMU_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(HiMATGMU_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(HiMATGMU_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(HiMATGMU_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(HiMATGMU_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(HiMATGMU_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n, gridDesci(n,:),&
               HiMATGMU_struc(n)%n111, HiMATGMU_struc(n)%n121, &
               HiMATGMU_struc(n)%n211, HiMATGMU_struc(n)%n221, &
               HiMATGMU_struc(n)%w111, HiMATGMU_struc(n)%w121, &
               HiMATGMU_struc(n)%w211, HiMATGMU_struc(n)%w221 )

       elseif ( trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear" ) then 
! === BUDGET BILINEAR INTERPOLATION ==== 
          allocate(HiMATGMU_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(HiMATGMU_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(HiMATGMU_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(HiMATGMU_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(HiMATGMU_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(HiMATGMU_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(HiMATGMU_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(HiMATGMU_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n, gridDesci(n,:), &
                 HiMATGMU_struc(n)%n112, HiMATGMU_struc(n)%n122, &
                 HiMATGMU_struc(n)%n212, HiMATGMU_struc(n)%n222, &
                 HiMATGMU_struc(n)%w112, HiMATGMU_struc(n)%w122, &
                 HiMATGMU_struc(n)%w212, HiMATGMU_struc(n)%w222 )

       end if   ! End interp option statement

    enddo

  end subroutine init_HIMATGMU

end module HiMATGMU_forcingMod
