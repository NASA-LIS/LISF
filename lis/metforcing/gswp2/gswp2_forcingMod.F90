!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gswp2_forcingMod
!BOP
! !MODULE: gswp2_forcingMod
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Soil
!  Wetness Project (GSWP2). GSWP2 forcing variables are produced
!  on a latlon 1degree grid at 3 hour intervals. 
!
!  The implemenatation in LIS has the derived data type {\tt gswp2\_struc} that
!  includes the variables that specify the runtime options
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[gswp2time1]
!    The nearest, previous 3 hour instance of the incoming 
!    data (as a real time). 
!  \item[gswp2time2]
!    The nearest, next 3 hour instance of the incoming 
!    data (as a real time).
!  \item[tair]
!    Directory containing the 2m air temperature data
!  \item[qair]
!    Directory containing the 2m specific humidity data
!  \item[psurf]
!    Directory containing the surface pressure data
!  \item[wind]
!    Directory containing the wind data
!  \item[rainf]
!    Directory containing the total precipitation data
!  \item[snowf]
!    Directory containing the total snowfall data
!  \item[swdown]
!    Directory containing the downward shortwave radiation data
!  \item[swdown]
!    Directory containing the downward longwave radiation data
!  \item[mi]
!    Number of points in the input grid
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
  public :: init_GSWP2      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gswp2_struc
!EOP

  type, public ::  gswp2_type_dec
     real     :: ts
     integer  :: ncold, nrold, vector_len   !AWIPS 212 dimensions
     real*8   :: gswp2time1,gswp2time2
     character(len=LIS_CONST_PATH_LEN) :: mfile
     character(len=LIS_CONST_PATH_LEN) :: tair
     character(len=LIS_CONST_PATH_LEN) :: qair
     character(len=LIS_CONST_PATH_LEN) :: psurf
     character(len=LIS_CONST_PATH_LEN) :: wind
     character(len=LIS_CONST_PATH_LEN) :: rainf
     character(len=LIS_CONST_PATH_LEN) :: snowf
     character(len=LIS_CONST_PATH_LEN) :: swdown
     character(len=LIS_CONST_PATH_LEN) :: lwdown
     character(len=LIS_CONST_PATH_LEN) :: rainf_c

     integer, allocatable   :: gindex(:,:)

     integer :: mi
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
     real, allocatable      ::  w112(:,:),w122(:,:)
     real, allocatable      ::  w212(:,:),w222(:,:)

     integer, allocatable   :: smask1(:,:)
     logical            :: fillflag1

     integer            :: findtime1, findtime2
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
     real, allocatable :: metdata3(:,:) 
  end type gswp2_type_dec
  
  type(gswp2_type_dec), allocatable :: gswp2_struc(:)
!EOP
contains
  
!BOP
!
! !ROUTINE: init_GSWP2
! \label{init_GSWP2}
! 
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_GSWP2(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !USES: 
    integer, intent(in)  :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GSWP2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gswp2](\ref{readcrd_gswp2}) \newline
!     reads the runtime options specified for GSWP2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[gswp2\_mask](\ref{gswp2_mask}) \newline
!    reads the GSWP2 mask
!  \end{description}
!EOP


    integer :: n 
    real  :: gridDesci(LIS_rc%nnest, 50)

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GSWP-2 forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif
    
    allocate(gswp2_struc(LIS_rc%nnest))
    call readcrd_gswp2()

    do n=1, LIS_rc%nnest
       gswp2_struc(n)%ts = 3*3600 
       call LIS_update_timestep(LIS_rc, n, gswp2_struc(n)%ts)
    enddo

    gridDesci = 0 

    LIS_rc%met_nf(findex) = 10 

    gswp2_struc(:)%ncold = 360
    gswp2_struc(:)%nrold = 150

    do n=1,LIS_rc%nnest

       allocate(gswp2_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gswp2_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gswp2_struc(n)%metdata3(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gswp2_struc(n)%metdata1 = 0
       gswp2_struc(n)%metdata2 = 0
       gswp2_struc(n)%metdata3 = 0

       gridDesci(n,1) = 0
       gridDesci(n,2) = gswp2_struc(n)%ncold
       gridDesci(n,3) = gswp2_struc(n)%nrold
       gridDesci(n,4) = -59.500
       gridDesci(n,5) = -179.500
       gridDesci(n,7) = 89.500
       gridDesci(n,8) = 179.500
       gridDesci(n,6) = 128
       gridDesci(n,9) = 1.000
       gridDesci(n,10) = 1.000
       gridDesci(n,20) = 0.0
       gswp2_struc(n)%mi = gswp2_struc(n)%ncold*gswp2_struc(n)%nrold
       gswp2_struc(n)%gswp2time1 = 3000.0
       gswp2_struc(n)%gswp2time2 = 0.0

!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(gswp2_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               gswp2_struc(n)%n111,gswp2_struc(n)%n121,&
               gswp2_struc(n)%n211,gswp2_struc(n)%n221,&
               gswp2_struc(n)%w111,gswp2_struc(n)%w121,&
               gswp2_struc(n)%w211,gswp2_struc(n)%w221)

          allocate(gswp2_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          gswp2_struc(n)%fillflag1 = .true. 

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(gswp2_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp2_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               gswp2_struc(n)%n111,gswp2_struc(n)%n121,&
               gswp2_struc(n)%n211,gswp2_struc(n)%n221,&
               gswp2_struc(n)%w111,gswp2_struc(n)%w121,&
               gswp2_struc(n)%w211,gswp2_struc(n)%w221)

          allocate(gswp2_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp2_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp2_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp2_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp2_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp2_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp2_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp2_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               gswp2_struc(n)%n112,gswp2_struc(n)%n122,&
               gswp2_struc(n)%n212,gswp2_struc(n)%n222,&
               gswp2_struc(n)%w112,gswp2_struc(n)%w122,&
               gswp2_struc(n)%w212,gswp2_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
          write(LIS_logunit,*) 'Neighbor interpolation is not supported'
          write(LIS_logunit,*) 'for GSWP2 forcing... Program stopping..'
          call LIS_endrun()
       endif
    enddo
    
    call gswp2_mask
    
  end subroutine init_GSWP2
end module gswp2_forcingMod
