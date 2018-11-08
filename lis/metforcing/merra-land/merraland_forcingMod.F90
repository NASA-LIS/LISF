!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module merraland_forcingMod
!BOP
! !MODULE: merraland_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the MERRA-Land forcing data.
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt merraland\_struc}
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
!  \item[merralandtime1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[merralandtime2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[merralanddir]
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
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for n. neighbor interpolation.
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
  public :: init_MERRALAND      !defines the native resolution of
                                !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: merraland_struc

!EOP
  type, public ::  merraland_type_dec
     real    :: ts
     integer :: ncold, nrold
     character*40 :: merralanddir   !MERRA-Land Forcing Directory
     real*8  :: merralandtime1,merralandtime2

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

     integer, allocatable   :: n113(:)
     integer            :: findtime1, findtime2
     logical            :: startFlag, dayFlag
     real, allocatable      :: merraforc1(:,:,:), merraforc2(:,:,:)

     integer            :: nvars
     integer            :: uselml

     real*8            :: ringtime
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:)      

  end type merraland_type_dec

  type(merraland_type_dec), allocatable :: merraland_struc(:)

contains

!BOP
!
! !ROUTINE: init_MERRALAND
! \label{init_MERRALAND}
!
! !REVISION HISTORY:
! 12 Oct 2009: Eric Kemp
! 22 Jul 2010: David Mocko, changed to hourly forcing
!
! !INTERFACE:
  subroutine init_MERRALAND(findex)

! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for MERRA-Land
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_merraland](\ref{readcrd_merraland}) \newline
!     reads the runtime options specified for MERRA-Land data
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
       write(LIS_logunit,*) '[ERR] Currently the MERRA-Land forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(merraland_struc(LIS_rc%nnest))
    call readcrd_merraland()

    do n=1, LIS_rc%nnest
       merraland_struc(n)%ts = 3600  !check
       call LIS_update_timestep(LIS_rc, n, merraland_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 14

    gridDesci = 0

    do n=1,LIS_rc%nnest

       allocate(merraland_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(merraland_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       merraland_struc(n)%metdata1 = 0
       merraland_struc(n)%metdata2 = 0

       gridDesci(n,1) = 0
       gridDesci(n,2) = merraland_struc(n)%ncold
       gridDesci(n,3) = merraland_struc(n)%nrold
       gridDesci(n,4) = -90.000
       gridDesci(n,5) = -180.000
       gridDesci(n,6) = 128
       gridDesci(n,7) = 90.000
       gridDesci(n,8) = 179.333333
       gridDesci(n,9) = 0.66666666667
       gridDesci(n,10) = 0.5
       gridDesci(n,20) = 0

       merraland_struc(n)%mi = merraland_struc(n)%ncold*merraland_struc(n)%nrold
! Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(merraland_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               merraland_struc(n)%n111,merraland_struc(n)%n121,&
               merraland_struc(n)%n211,merraland_struc(n)%n221,&
               merraland_struc(n)%w111,merraland_struc(n)%w121,&
               merraland_struc(n)%w211,merraland_struc(n)%w221)
       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then

          allocate(merraland_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merraland_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               merraland_struc(n)%n111,merraland_struc(n)%n121,&
               merraland_struc(n)%n211,merraland_struc(n)%n221,&
               merraland_struc(n)%w111,merraland_struc(n)%w121,&
               merraland_struc(n)%w211,merraland_struc(n)%w221)

          allocate(merraland_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merraland_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merraland_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merraland_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merraland_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merraland_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merraland_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merraland_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               merraland_struc(n)%n112,merraland_struc(n)%n122,&
               merraland_struc(n)%n212,merraland_struc(n)%n222,&
               merraland_struc(n)%w112,merraland_struc(n)%w122,&
               merraland_struc(n)%w212,merraland_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then

          allocate(merraland_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci(n,:),&
               merraland_struc(n)%n113)
       else
          write(LIS_logunit,*) 'Interpolation option '// &
               trim(LIS_rc%met_interp(findex))//&
               ' for MERRA-Land forcing is not supported'
          call LIS_endrun()

       endif

       call LIS_registerAlarm("MERRA-Land forcing alarm",&
            86400.0,86400.0)
       merraland_struc(n)%startFlag = .true.
       merraland_struc(n)%dayFlag = .true.

       merraland_struc(n)%nvars = 14

       allocate(merraland_struc(n)%merraforc1(&
            merraland_struc(n)%nvars, 24, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(merraland_struc(n)%merraforc2(&
            merraland_struc(n)%nvars, 24, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       merraland_struc(n)%merraforc1 = LIS_rc%udef
       merraland_struc(n)%merraforc2 = LIS_rc%udef

    enddo
  end subroutine init_MERRALAND
end module merraland_forcingMod

