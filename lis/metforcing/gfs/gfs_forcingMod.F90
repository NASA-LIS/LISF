!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gfs_forcingMod
!BOP
! !MODULE: gfs_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Forecast
!  System (GFS) developed at the Environmental Modeling
!  Center (EMC) of NOAA/NCEP. GFS forcing variables are produced
!  on a quadratic gaussian grid. LIS uses the 00, 03, and 06 forecast
!  files. The forecasts are produced at 6 hr intervals. 
!  The resolution of GFS forcing varies as follows:
!   
!   upto 2000/1/24          :   T126 (384x190)  grid \newline
!   2001/1/24  - 2002/10/29 :   T170 (512x256)  grid \newline
!   2002/10/29 - 2005/5/31  :   T254 (768x384)  grid \newline
!   2005/5/31  - 2010/7/28  :   T382 (1152x576) grid \newline
!   2010/7/28  - 2015/1/14  :   T574 (1760x880) grid \newline
!   2015/1/14  - current    :  T1534 (3072x1536) grid \newline
!
!  The implementation in LIS has the derived data type {\tt gfs\_struc} that
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the GFS data
!  \item[gfstime1]
!    The nearest, previous 3 hour instance of the incoming 
!    data (as a real time). 
!  \item[gfstime2]
!    The nearest, next 3 hour instance of the incoming 
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch GFS resolution to T126
!  \item[griduptime2]
!    The time to switch GFS resolution to T170
!  \item[griduptime3]
!    The time to switch GFS resolution to T254
!  \item[griduptime4]
!    The time to switch GFS resolution to T382
!  \item[griduptime5]
!    The time to switch GFS resolution to T574
!  \item[griduptime6]
!    The time to switch GFS resolution to T1534
!  \item[gfsdir]
!    Directory containing the input data
!  \item[elevfile]
!    File with the elevation definition for the input data. 
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
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \end{description}
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GFS      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gfs_struc
!EOP

  type, public :: gfs_type_dec
     real         :: ts
     integer      :: ncold, nrold   !AWIPS 212 dimensions
     integer      :: nmif
     character(len=LIS_CONST_PATH_LEN) :: gfsdir  !GFS Forcing Directory
     character(len=LIS_CONST_PATH_LEN) :: elevfile
     character(len=LIS_CONST_PATH_LEN) :: t126elevfile
     character(len=LIS_CONST_PATH_LEN) :: t170elevfile
     character(len=LIS_CONST_PATH_LEN) :: t254elevfile
     character(len=LIS_CONST_PATH_LEN) :: t382elevfile
     character(len=LIS_CONST_PATH_LEN) :: t574elevfile
     character(len=LIS_CONST_PATH_LEN) :: t1534elevfile
     real*8       :: gfstime1,gfstime2
     real*8       :: griduptime1,griduptime2,griduptime3,griduptime4,griduptime5,griduptime6
     logical      :: gridchange1,gridchange2,gridchange3,gridchange4,gridchange5,gridchange6
     integer      :: mi

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
     integer                :: findtime1, findtime2 

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type gfs_type_dec
  
  type(gfs_type_dec), allocatable :: gfs_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_GFS
!  \label{init_GFS}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_GFS(findex)
! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
   use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
! !AGRUMENTS: 
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GFS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}). Based on the GFS data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights. The dates of the GFS resolution switches are also 
!  defined in this routine. 
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gfs](\ref{readcrd_gfs}) \newline
!     reads the runtime options specified for GFS data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP

    integer :: n 
    real    :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1,tdoy
    real    :: upgmt, tgmt

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GFS forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(gfs_struc(LIS_rc%nnest))

    call readcrd_gfs()

    do n=1, LIS_rc%nnest
       gfs_struc(n)%ts = 6*3600 
       call LIS_update_timestep(LIS_rc, n, gfs_struc(n)%ts)
    enddo
    
    LIS_rc%met_nf(findex) = 9

    do n=1,LIS_rc%nnest

       allocate(gfs_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gfs_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gfs_struc(n)%metdata1 = 0
       gfs_struc(n)%metdata2 = 0

       gfs_struc(n)%findtime1 = 0 
       gfs_struc(n)%findtime2 = 0 
       ! This grid is good for some time in the 1980's.
       ! Look up the exact dates.
       gridDesci = 0 
       gridDesci(1) = 4
       gridDesci(2) = 192
       gridDesci(3) = 94
       gridDesci(4) = 88.542
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) =  -88.542
       gridDesci(8) = -1.875
       gridDesci(9) = 1.875
       gridDesci(10) = 47
       gridDesci(20) = 0
       gfs_struc(n)%mi = gfs_struc(n)%ncold*gfs_struc(n)%nrold
       
       ! This grid is good for some time in the 1990's.
       ! Look up the exact dates.
       yr1 = 1991
       mo1 = 01
       da1 = 01
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time( gfs_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2000
       mo1 = 01
       da1 = 24
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time( gfs_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
       
       yr1 = 2002     !grid update time ~ 0.469
       mo1 = 10
       da1 = 29
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gfs_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2005     !grid update time ~ 0.313
       mo1 = 05
       da1 = 31
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gfs_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2010     !grid update time ~ 0.205
       mo1 = 07
       da1 = 28
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gfs_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       yr1 = 2015     !grid update time ~ 0.117
       mo1 = 01
       da1 = 14
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gfs_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       gfs_struc(n)%gridchange1 = .true.
       gfs_struc(n)%gridchange2 = .true.
       gfs_struc(n)%gridchange3 = .true.
       gfs_struc(n)%gridchange4 = .true.
       gfs_struc(n)%gridchange5 = .true.
       gfs_struc(n)%gridchange6 = .true.

!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
          allocate(gfs_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               gfs_struc(n)%n111,gfs_struc(n)%n121,&
               gfs_struc(n)%n211,gfs_struc(n)%n221,&
               gfs_struc(n)%w111,gfs_struc(n)%w121,&
               gfs_struc(n)%w211,gfs_struc(n)%w221)
       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(gfs_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gfs_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               gfs_struc(n)%n111,gfs_struc(n)%n121,&
               gfs_struc(n)%n211,gfs_struc(n)%n221,&
               gfs_struc(n)%w111,gfs_struc(n)%w121,&
               gfs_struc(n)%w211,gfs_struc(n)%w221)

          allocate(gfs_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gfs_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gfs_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gfs_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gfs_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gfs_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gfs_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gfs_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               gfs_struc(n)%n112,gfs_struc(n)%n122,&
               gfs_struc(n)%n212,gfs_struc(n)%n222,&
               gfs_struc(n)%w112,gfs_struc(n)%w122,&
               gfs_struc(n)%w212,gfs_struc(n)%w222)
       endif
    enddo
  end subroutine init_GFS
end module gfs_forcingMod

