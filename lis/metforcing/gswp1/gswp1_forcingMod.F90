!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !MODULE: gswp1_forcingMod
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Soil
!  Wetness Project (GSWP). GSWP forcing variables are produced
!  on a latlon 1degree grid at 3 hour intervals. 
!
!  The implemenatation in LIS has the derived data type {\tt gswp\_struc} that
!  includes the variables that specify the runtime options
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[gswptime1]
!    The nearest, previous 3 hour instance of the incoming 
!    data (as a real time). 
!  \item[gswptime2]
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
!EOP
module gswp1_forcingMod

  use LIS_constantsMod, only : LIS_CONST_PATH_LEN 

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_gswp1      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gswp1_struc

  type, public :: gswp1_type_dec
     real          :: ts
     integer       :: mi 
     integer       :: ncold
     integer       :: nrold
     character(len=LIS_CONST_PATH_LEN) :: gswp1dir
     real*8        :: gswp1time1
     real*8        :: gswp1time2

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
  end type gswp1_type_dec

  type(gswp1_type_dec), allocatable :: gswp1_struc(:)

contains

!BOP
!
! !ROUTINE: init_gswp1
! \label{init_gswp1}
!
! !REVISION HISTORY:
! 11 Dec 2003: Sujay Kumar, Initial Specification
!
! !INTERFACE:
  subroutine init_gswp1(findex)
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !AGRUMENTS: 
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GSWP1
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gswp1](\ref{readcrd_gswp1}) \newline
!     reads the runtime options specified for GSWP1 data
!  \end{description}
!EOP

    real     :: gridDesci(LIS_rc%nnest, 50)
    integer  :: n 

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GSWP-1 forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(gswp1_struc(LIS_rc%nnest))    
    call readcrd_gswp1()

    do n=1, LIS_rc%nnest
       gswp1_struc(n)%ts = 3600 
       call LIS_update_timestep(LIS_rc, n, gswp1_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 9

    gswp1_struc(:)%ncold = 360
    gswp1_struc(:)%nrold = 150

    gridDesci = 0 
    
    do n=1,LIS_rc%nnest

       allocate(gswp1_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gswp1_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gswp1_struc(n)%metdata1 = 0
       gswp1_struc(n)%metdata2 = 0

       gridDesci(n,1) = 0
       gridDesci(n,2) = gswp1_struc(n)%ncold
       gridDesci(n,3) = gswp1_struc(n)%nrold
       gridDesci(n,4) = -59.500
       gridDesci(n,5) = -179.500
       gridDesci(n,7) = 89.500
       gridDesci(n,8) = 179.500
       gridDesci(n,6) = 128
       gridDesci(n,9) = 1.000
       gridDesci(n,10) = 1.000
       gridDesci(n,20) = 0.0
       gswp1_struc(n)%mi = gswp1_struc(n)%ncold*gswp1_struc(n)%nrold
       gswp1_struc(n)%gswp1time1 = 1000.0
       gswp1_struc(n)%gswp1time2 = 0.0

!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(gswp1_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               gswp1_struc(n)%n111,gswp1_struc(n)%n121,&
               gswp1_struc(n)%n211,gswp1_struc(n)%n221,&
               gswp1_struc(n)%w111,gswp1_struc(n)%w121,&
               gswp1_struc(n)%w211,gswp1_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(gswp1_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gswp1_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n,gridDesci(n,:),&
               gswp1_struc(n)%n111,gswp1_struc(n)%n121,&
               gswp1_struc(n)%n211,gswp1_struc(n)%n221,&
               gswp1_struc(n)%w111,gswp1_struc(n)%w121,&
               gswp1_struc(n)%w211,gswp1_struc(n)%w221)

          allocate(gswp1_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp1_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp1_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp1_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp1_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp1_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp1_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gswp1_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               gswp1_struc(n)%n112,gswp1_struc(n)%n122,&
               gswp1_struc(n)%n212,gswp1_struc(n)%n222,&
               gswp1_struc(n)%w112,gswp1_struc(n)%w122,&
               gswp1_struc(n)%w212,gswp1_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
          write(LIS_logunit,*) 'Neighbor interpolation is not supported'
          write(LIS_logunit,*) 'for GSWP1 forcing... Program stopping..'
          call LIS_endrun()
       endif
    enddo
    
  end subroutine init_gswp1
end module gswp1_forcingMod

