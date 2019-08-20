!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module geos_forcingMod
!BOP
! !MODULE: geos_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the GEOS forcing data. 
!  The data is global 1 degree dataset in latlon
!  projection, and at 3 hourly intervals. The derived
!  data type {\tt geos\_struc}
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the ECMWF data
!  \item[geostime1]
!    The nearest, previous 3 hour instance of the incoming 
!    data (as a real time). 
!  \item[geostime2]
!    The nearest, next 3 hour instance of the incoming 
!    data (as a real time).
!  \item[griduptime]
!    The time to switch GEOS resolution 
!  \item[geosdir]
!    Directory containing the input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
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
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_geos      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: geos_struc

!EOP
  type, public ::  geos_type_dec
     real    :: ts
     integer :: ncold, nrold   !AWIPS 212 dimensions
     integer :: nmif
     character*40 :: geosdir   !GEOS Forcing Directory
     real*8  :: geostime1,geostime2
     real*8  :: griduptime1, griduptime2
     logical :: gridchange1, gridchange2

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
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)
     integer            :: findtime1, findtime2
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type geos_type_dec
  
  type(geos_type_dec), allocatable :: geos_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_geos
! \label{init_geos}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_geos(findex)
! !USES: 
   use LIS_coreMod,    only:  LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
   use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS: 
    integer,  intent(in)    :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GEOS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_geos](\ref{readcrd_geos}) \newline
!     reads the runtime options specified for GEOS data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real :: gridDesci(LIS_rc%nnest,50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GEOS (older version) forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif


    allocate(geos_struc(LIS_rc%nnest))
    call readcrd_geos()

    do n=1, LIS_rc%nnest
       geos_struc(n)%ts = 3*3600 
       call LIS_update_timestep(LIS_rc, n, geos_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 10

    geos_struc(:)%ncold = 360
    geos_struc(:)%nrold = 181
    geos_struc(:)%nmif  = 13

    gridDesci = 0 

    do n=1,LIS_rc%nnest

       allocate(geos_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(geos_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       geos_struc(n)%metdata1 = 0
       geos_struc(n)%metdata2 = 0

       gridDesci(n,1) = 0
       gridDesci(n,2) = geos_struc(n)%ncold
       gridDesci(n,3) = geos_struc(n)%nrold
       gridDesci(n,4) = -90.000
       gridDesci(n,5) = -180.000
       gridDesci(n,7) = 90.000
       gridDesci(n,8) = 179.000
       gridDesci(n,6) = 128
       gridDesci(n,9) = 1.000
       gridDesci(n,10) = 1.000
       gridDesci(n,20) = 0
       geos_struc(n)%mi = geos_struc(n)%ncold*geos_struc(n)%nrold

       yr1 = 2002  !grid update time
       mo1 = 10
       da1 = 01
       hr1 = 0; mn1 = 0; ss1 = 0
       call LIS_date2time(geos_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)

       yr1 = 2007  !grid update time
       mo1 = 3
       da1 = 2
       hr1 = 0; mn1 = 0; ss1 = 0
       call LIS_date2time(geos_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)

       geos_struc(n)%gridchange1 = .true.
       geos_struc(n)%gridchange2 = .true. 
!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(geos_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))


          call bilinear_interp_input(n,gridDesci(n,:),&
               geos_struc(n)%n111,geos_struc(n)%n121,&
               geos_struc(n)%n211,geos_struc(n)%n221,&
               geos_struc(n)%w111,geos_struc(n)%w121,&
               geos_struc(n)%w211,geos_struc(n)%w221)
       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(geos_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(geos_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))


          call bilinear_interp_input(n,gridDesci(n,:),&
               geos_struc(n)%n111,geos_struc(n)%n121,&
               geos_struc(n)%n211,geos_struc(n)%n221,&
               geos_struc(n)%w111,geos_struc(n)%w121,&
               geos_struc(n)%w211,geos_struc(n)%w221)

          allocate(geos_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(geos_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(geos_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(geos_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(geos_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(geos_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(geos_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(geos_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               geos_struc(n)%n112,geos_struc(n)%n122,&
               geos_struc(n)%n212,geos_struc(n)%n222,&
               geos_struc(n)%w112,geos_struc(n)%w122,&
               geos_struc(n)%w212,geos_struc(n)%w222)
       endif
    enddo
  end subroutine init_geos
end module geos_forcingMod
