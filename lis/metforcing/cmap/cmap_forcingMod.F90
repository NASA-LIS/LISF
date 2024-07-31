!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module cmap_forcingMod
!BOP
! !MODULE: cmap_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  Climate Prediction Center (CPC)'s Merged Analysis of Precipitation
!  (CMAP). CMAP merges gauge measurements and satellite estimates
!  including GPI, OPI, SSM/I to produce a global 2.5 degree pentad
!  precipitation analysis. The 6-hourly product obtained by disaggregating 
!  CMAP using the GDAS forcing is used in this implementation. 
! 
!   upto 2000/1/24          :   T126 (384x190)  grid
!   2001/1/24  - 2002/10/29 :   T170 (512x256)  grid
!   2002/10/29 - 2005/5/31  :   T254 (768x384)  grid
!   2005/5/31  - 2012/9/30  :   T382 (1152x576) grid ~~ CMAP only
!   2012/10/01 onwards      :   T574 (1760x880) grid ~~ CMAP only
!  Original GDAS grid change:
!   2010/7/28  onwards      :   T574 (1760x880) grid ~~ Original GDAS
!
!  The implementation in LIS has the derived data type {\tt cmap\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation.
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[cmapdir]
!    Directory containing the input data
!  \item[cmaptime]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution to T170
!  \item[mi]
!    Number of points in the input grid
!  \item[n11,n121,n211,n221]
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
  public :: init_CMAP      !defines the native resolution of 
                           !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: cmap_struc

!EOP
 
  type, public ::  cmap_type_dec
     real                   :: ts
     integer                :: ncold
     integer                :: nrold  
     character(len=LIS_CONST_PATH_LEN) :: cmapdir  
     character*50           :: met_interp
     real*8                 :: cmaptime
     real*8                 :: griduptime1
     real*8                 :: griduptime2
     real*8                 :: griduptime3
     real*8                 :: griduptime4
     real*8                 :: griduptime5
     logical                :: gridchange1
     logical                :: gridchange2
     logical                :: gridchange3
     logical                :: gridchange4
     logical                :: gridchange5
     integer                :: mi

   ! Bilinear interp weights and indices
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)

   ! Budget-Bilinear (conserve) interp weights and indices
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real,    allocatable   :: w112(:,:),w122(:,:)
     real,    allocatable   :: w212(:,:),w222(:,:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type cmap_type_dec

  type(cmap_type_dec), allocatable :: cmap_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_CMAP
! \label{init_CMAP}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 02Dec2014: KR Arsenault: Added new grid change update (~2012)
! 
! !INTERFACE:
  subroutine init_CMAP(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod

    implicit none
    integer, intent(in) :: findex

! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for CMAP
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_cmap](\ref{readcrd_cmap}) \newline
!     reads the runtime options specified for CMAP data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!  \end{description}
!
!EOP

    real    :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    integer :: n 

    allocate(cmap_struc(LIS_rc%nnest))

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the CMAP forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

!- Read in LIS Config file entries:
    call readcrd_cmap()

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

    do n=1, LIS_rc%nnest
       cmap_struc(n)%ts = 21600
       call LIS_update_timestep(LIS_rc, n, cmap_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 2

 !- Set interp arrays for reinterpolation later:
    do n=1,LIS_rc%nnest

       allocate(cmap_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(cmap_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       cmap_struc(n)%metdata1 = 0
       cmap_struc(n)%metdata2 = 0
!------------------------------------------------------------------------ 
! CMAP product start time to use JonG's data. 1979010100 for 00Z-06Z.
!   This grid is good for some time in the 1980's.
!   Look up the exact dates.
!------------------------------------------------------------------------
!- Initialize first available grid domain/res for CMAP/GDAS:
!   2001/1/24 - 2002/10/29 : T170 (512x256)  grid
!- CMAP/GDAS Grid description:
       cmap_struc(n)%ncold = 512
       cmap_struc(n)%nrold = 256

       gridDesci = 0
       gridDesci(1) = 4
       gridDesci(2) = 512
       gridDesci(3) = 256
       gridDesci(4) = 89.463
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) = -89.463
       gridDesci(8) = -0.703
       gridDesci(9) = 0.703
       gridDesci(10) = 128
       gridDesci(20) = 64
       cmap_struc(n)%mi = gridDesci(2)*gridDesci(3)

     ! This grid is good for some time in the 1990's.
     ! Look up the exact dates.
       yr1 = 1991     !grid update time
       mo1 = 01
       da1 = 01
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(cmap_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! == 2001/1/24 - 2002/10/29 : T170 (512x256)  grid
       yr1 = 2000     !grid update time
       mo1 = 01
       da1 = 24
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(cmap_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! ==  2002/10/29 - 2005/5/31 : T254 (768x384)  grid
       yr1 = 2002     !grid update time
       mo1 = 10
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(cmap_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )
       
! == 2005/5/31 - 2012/9/30  : T382 (1152x576) grid ~~ CMAP only
       yr1 = 2005     !grid update time
       mo1 = 05
       da1 = 31
       hr1 = 0
       mn1 = 0; ss1 = 0
       call LIS_date2time(cmap_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! == 2012/10/01 onwards  : T574 (1760x880) grid ~~ CMAP only
       yr1 = 2012     !grid update time
       mo1 = 09
       da1 = 30
       hr1 = 18      ! -- ORIGINAL HOUR SET IN CODE
!       hr1 = 0      ! When actual grid change occurs
       mn1 = 0; ss1 = 0
       call LIS_date2time(cmap_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

       cmap_struc(n)%gridchange1 = .true.
       cmap_struc(n)%gridchange2 = .true.
       cmap_struc(n)%gridchange3 = .true.
       cmap_struc(n)%gridchange4 = .true.
       cmap_struc(n)%gridchange5 = .true.

       ! Setting up weights for Interpolation
       call cmap_reset_interp_input(n, findex, gridDesci)
    enddo
  end subroutine init_CMAP

end module cmap_forcingMod
