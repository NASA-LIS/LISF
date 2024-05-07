!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module stg2_forcingMod

!BOP
! !MODULE: stg2_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  National Center for Environmental Prediction (NCEP) Stage II  
!  (STAGE2) Doppler Radar+gage product.  The Stage II is a national
!  level product, on an hourly interval, and supplements mainly
!  the which base forcing precipitation (e.g., NLDAS) is being used.
! 
!  The implementation in LDT has the derived data type {\tt stg2\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncol]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrow]
!    Number of rows (along the north south dimension) for the input data
!  \item[stg2dir]
!    Directory containing the input data
!  \item[stg2time]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LDT, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LDT, for conservative interpolation.
!  \end{description}
!
! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_STG2      !defines the native resolution of 
                           !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: stg2_struc

!EOP

  type, public :: stg2_type_dec
     real               :: ts
     integer            :: ncol                 ! Number of cols
     integer            :: nrow                 ! Number of rows
     character(len=LDT_CONST_PATH_LEN)       :: stg2dir              ! STAGE II Directory
     real*8             :: stg2time             ! Nearest hourly instance of incoming file
     real*8             :: griduptime1          ! Designated time of STAGEII grid change
     logical            :: gridchange1          ! Flag for when grid change occurs
     integer            :: mi                   ! Number of points in the input grid

! == Arrays for Bilinear Interpolation option (=1)
     integer, allocatable  :: n111(:), n121(:)
     integer, allocatable  :: n211(:), n221(:)
     real, allocatable     :: w111(:), w121(:)
     real, allocatable     :: w211(:), w221(:)
! == Arrays for Budget Bilinear Interpolation option (=2)
     integer, allocatable  :: n112(:,:), n122(:,:) ! Spatial Interpolation weights
     integer, allocatable  :: n212(:,:), n222(:,:) ! " "
     real,    allocatable  :: w112(:,:), w122(:,:) ! " "
     real,    allocatable  :: w212(:,:), w222(:,:) ! " "

  end type stg2_type_dec

  type(stg2_type_dec), allocatable :: stg2_struc(:)       
! STAGE II Main Pointer array (number of LDT_domains)

contains
  
!BOP
!
! !ROUTINE: init_STG2
! \label{init_STG2}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 25May2006: Kristi Arsenault; Data and code implementation
! 
! !INTERFACE:
  subroutine init_STG2(findex)

! !USES: 
   use LDT_coreMod,    only : LDT_rc, LDT_domain
   use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep
   use LDT_logMod,     only : LDT_logunit, LDT_endrun

   implicit none
   integer,   intent(in) :: findex

! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for STAGE2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_stg2](\ref{readcrd_stg2}) \newline
!     reads the runtime options specified for STAGE2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!     computes the neighbor, weights for bilinear interpolation
!   \item[LDT\_date2time](\ref{LDT_date2time}) \newline
!     converts date to the real time format - time of grid change 
!  \end{description}
!
!EOP
    integer :: n
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    real    :: gridDesci(20)
! _____________________________________________

! - Stage II Data structure -- Allocate for different LDT_domains
    allocate ( stg2_struc(LDT_rc%nnest) )

    write(LDT_logunit,fmt=*)"MSG: Initializing Stage-II forcing grid ... "

! - Retrieve Stage II Forcing Dataset Directory Name Location from ldt.config
    call readcrd_stg2()

    LDT_rc%met_nf(findex) = 2
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .false.
    LDT_rc%met_proj(findex)  = "polar"

    stg2_struc%ncol = 1160.
    stg2_struc%nrow = 880.
    LDT_rc%met_nc(findex) = stg2_struc(1)%ncol
    LDT_rc%met_nr(findex) = stg2_struc(1)%nrow

 !- STAGE II Grid description:
    LDT_rc%met_gridDesc(findex,1)  = 5       ! Projection type (UPS)
    LDT_rc%met_gridDesc(findex,2)  = stg2_struc(1)%ncol  ! X-dir amount of points
    LDT_rc%met_gridDesc(findex,3)  = stg2_struc(1)%nrow  ! Y-dir amount of points
    LDT_rc%met_gridDesc(findex,4)  = 22.7736000          ! Starting latitude point
    LDT_rc%met_gridDesc(findex,5)  = -120.376            ! Starting longitude point
    LDT_rc%met_gridDesc(findex,6)  = 8.0
    LDT_rc%met_gridDesc(findex,7)  = 0                   ! Orientation
    LDT_rc%met_gridDesc(findex,8)  = 4.7625              ! X-spacing length (kms)
    LDT_rc%met_gridDesc(findex,9)  = 4.7625              ! Y-spacing length (kms)
    LDT_rc%met_gridDesc(findex,10) = 60.0                ! True lat  (???)
    LDT_rc%met_gridDesc(findex,11) = -105.0              ! Standard longitude (???)
    LDT_rc%met_gridDesc(findex,20) = 64

    gridDesci(:) = LDT_rc%met_gridDesc(findex,:)  

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

! - Initial Dataset Date/Time:
    yr1 = 1996
    mo1 = 05
    da1 = 01
    hr1 = 00; mn1 = 0; ss1 = 0

! - Set time step:
    do n=1, LDT_rc%nnest
       stg2_struc(n)%ts = 3600
       call LDT_update_timestep(LDT_rc, n, stg2_struc(n)%ts)
    enddo

    do n=1,LDT_rc%nnest

       stg2_struc(n)%mi = stg2_struc(n)%ncol * stg2_struc(n)%nrow

       select case( LDT_rc%met_gridtransform(findex) )

         ! = BILINEAR INTERPOLATION =
         case( "bilinear" )
          allocate(stg2_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(stg2_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(stg2_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(stg2_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(stg2_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(stg2_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(stg2_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(stg2_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input(n, gridDesci(:), &
               stg2_struc(n)%n111, stg2_struc(n)%n121, &
               stg2_struc(n)%n211, stg2_struc(n)%n221, &
               stg2_struc(n)%w111, stg2_struc(n)%w121, &
               stg2_struc(n)%w211, stg2_struc(n)%w221 )

         ! = Budget-Bilinear Interpolation =
         case( "budget-bilinear" )
          allocate(stg2_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(stg2_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(stg2_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(stg2_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(stg2_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(stg2_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(stg2_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(stg2_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

          call conserv_interp_input( n, gridDesci(:), &
               stg2_struc(n)%n112, stg2_struc(n)%n122, &
               stg2_struc(n)%n212, stg2_struc(n)%n222, &
               stg2_struc(n)%w112, stg2_struc(n)%w122, &
               stg2_struc(n)%w212, stg2_struc(n)%w222 )
          
       end select   ! End interp option statement

   !-- 1st CHANGE IN STAGE II Grid Domain
       yr1 = 2001      !grid update time
       mo1 = 12 
       da1 = 31 
       hr1 = 23; mn1 = 0; ss1 = 0
       call LDT_date2time( stg2_struc(n)%griduptime1, updoy, upgmt, &
                yr1, mo1, da1, hr1, mn1, ss1 )

       stg2_struc(n)%gridchange1 = .true.

    enddo

  end subroutine init_STG2

end module stg2_forcingMod
