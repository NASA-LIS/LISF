!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module pildas_forcingMod
!BOP
! !MODULE: pildas_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the PILDAS forcing data. 
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt pildas\_struc}
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
!  \item[pildastime1]
!    The nearest, previous 1 hour instance of the incoming 
!    data (as a real time). 
!  \item[pildastime2]
!    The nearest, next 1 hour instance of the incoming 
!    data (as a real time).
!  \item[pildasdir]
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
  public :: init_pildas     ! defines the native resolution of 
                            ! the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: pildas_struc

!EOP
  type, public ::  pildas_type_dec
     real    :: ts
     integer :: nc, nr
     integer :: version
     integer :: uselml
     character*40 :: pildasdir
     real*8 :: pildastime1,pildastime2
     real   :: gmt1, gmt2

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
     logical            :: startFlag

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type pildas_type_dec
  
  type(pildas_type_dec), allocatable :: pildas_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_pildas
! \label{init_pildas}
!
! !REVISION HISTORY: 
! 25 Apr 2013: Sujay Kumar, initial specification
! 14 Jul 2016: Mahdi Navari, Modified for PILDAS
! 
! !INTERFACE:
  subroutine init_pildas(findex)

! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod

    implicit none
! !AGRUMENTS: 
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for PILDAS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_pildas](\ref{readcrd_pildas}) \newline
!     reads the runtime options specified for PILDAS data
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
       write(LIS_logunit,*) '[ERR] Currently the PILDAS forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(pildas_struc(LIS_rc%nnest))
    call readcrd_pildas()
    
    do n=1, LIS_rc%nnest
       pildas_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, pildas_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 10

    gridDesci = 0 

    do n=1,LIS_rc%nnest

       allocate(pildas_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(pildas_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       pildas_struc(n)%metdata1 = 0
       pildas_struc(n)%metdata2 = 0

       gridDesci(n,1) = 0
       gridDesci(n,2) = pildas_struc(n)%nc
       gridDesci(n,3) = pildas_struc(n)%nr
       gridDesci(n,4) = 33.1875
       gridDesci(n,5) = -103.3125
       gridDesci(n,6) = 128
       gridDesci(n,7) = 39.3125
       gridDesci(n,8) = -94.6875
       gridDesci(n,9) = 0.125
       gridDesci(n,10) = 0.125
       gridDesci(n,20) = 0

       pildas_struc(n)%mi = pildas_struc(n)%nc*pildas_struc(n)%nr
!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
          allocate(pildas_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               pildas_struc(n)%n111,pildas_struc(n)%n121,&
               pildas_struc(n)%n211,pildas_struc(n)%n221,&
               pildas_struc(n)%w111,pildas_struc(n)%w121,&
               pildas_struc(n)%w211,pildas_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
          allocate(pildas_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(pildas_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               pildas_struc(n)%n111,pildas_struc(n)%n121,&
               pildas_struc(n)%n211,pildas_struc(n)%n221,&
               pildas_struc(n)%w111,pildas_struc(n)%w121,&
               pildas_struc(n)%w211,pildas_struc(n)%w221)

          allocate(pildas_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(pildas_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(pildas_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(pildas_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(pildas_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(pildas_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(pildas_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(pildas_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               pildas_struc(n)%n112,pildas_struc(n)%n122,&
               pildas_struc(n)%n212,pildas_struc(n)%n222,&
               pildas_struc(n)%w112,pildas_struc(n)%w122,&
               pildas_struc(n)%w212,pildas_struc(n)%w222)

       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 

          allocate(pildas_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci(n,:),&
               pildas_struc(n)%n113)
       else
          write(LIS_logunit,*) 'Interpolation option '// &
               trim(LIS_rc%met_interp(findex))//&
               ' for PILDAS forcing is not supported'
          call LIS_endrun()

       endif

       call LIS_registerAlarm("PILDAS forcing alarm",&
            86400.0,86400.0)
       pildas_struc(n)%startFlag = .true. 

    enddo
  end subroutine init_pildas
end module pildas_forcingMod
