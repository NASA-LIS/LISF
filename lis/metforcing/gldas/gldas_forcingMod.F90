!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gldas_forcingMod
!BOP
! !MODULE: gldas_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the GLDAS forcing data produced by merging
!  GEOS, AGRMET radiation and CMAP precipitation
!  
!  The forcing data consists of the following meteorological variables. 
!  \begin{description}
!   \item[2m air temperature]
!   \item[2m specific humidity]
!   \item[incident downward shortwave radiation]
!   \item[incident downward longwave radiation]
!   \item[eastward wind]
!   \item[northward wind]
!   \item[surface pressure]
!   \item[rainfall rate]
!   \item[snowfall rate]
!   \item[convective rainfall rate]
!  \end{description}
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_gldas      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gldas_struc
!EOP

  type, public :: gldas_type_dec
     real         :: ts
     integer      :: ncold, nrold   !AWIPS 212 dimensions
     character(len=LIS_CONST_PATH_LEN) :: gldasdir       !GLDAS Forcing Directory
     real*8       :: gldastime1,gldastime2
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

     integer, allocatable   :: n113(:)

     integer, allocatable   :: smask1(:,:)
     integer, allocatable   :: colfill(:,:)
     integer, allocatable   :: rowfill(:,:)
     logical            :: fillflag1
     integer            :: findtime1, findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type gldas_type_dec
  
  type(gldas_type_dec), allocatable :: gldas_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_gldas
!  \label{init_gldas}
!
! !REVISION HISTORY: 
!  19 Sept 2008: Sujay Kumar: Initial Implementation
! 
! !INTERFACE:
  subroutine init_gldas(findex)
! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_update_timestep
   use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GLDAS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gldas](\ref{readcrd_gldas}) \newline
!     reads the runtime options specified for GLDAS data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[neighbor\_interp\_input](\ref{neighbor_interp_input}) \newline
!    computes the computational weights for neighbor interpolation
!  \end{description}
!EOP

    integer :: n 
    real :: gridDesci(50)

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GLDAS (older version) forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(gldas_struc(LIS_rc%nnest))

    call readcrd_gldas()

    do n=1, LIS_rc%nnest
       gldas_struc(n)%ts = 3*3600 
       call LIS_update_timestep(LIS_rc, n, gldas_struc(n)%ts)
    enddo
    
    LIS_rc%met_nf(findex) = 10 

    gldas_struc(:)%ncold = 144
    gldas_struc(:)%nrold = 76
    
    gridDesci = 0

    do n=1,LIS_rc%nnest

       allocate(gldas_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gldas_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gldas_struc(n)%metdata1 = 0
       gldas_struc(n)%metdata2 = 0

       gridDesci(1)  = 0
       gridDesci(2)  = gldas_struc(n)%ncold
       gridDesci(3)  = gldas_struc(n)%nrold
       gridDesci(4)  = -60.0
       gridDesci(5)  = -180.0
       gridDesci(6)  = 128
       gridDesci(7)  = 90.0
       gridDesci(8)  = 177.5 ! 180.00
       gridDesci(9)  = 2.50
       gridDesci(10) = 2.00
       gridDesci(20) = 0

       gldas_struc(n)%mi = gldas_struc(n)%ncold*gldas_struc(n)%nrold
       
!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(gldas_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call bilinear_interp_input(n,gridDesci,&
               gldas_struc(n)%n111,gldas_struc(n)%n121,&
               gldas_struc(n)%n211,gldas_struc(n)%n221,&
               gldas_struc(n)%w111,gldas_struc(n)%w121,&
               gldas_struc(n)%w211,gldas_struc(n)%w221)

          allocate(gldas_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%colfill(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%rowfill(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          gldas_struc(n)%fillflag1 = .true. 

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(gldas_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci,&
               gldas_struc(n)%n111,gldas_struc(n)%n121,&
               gldas_struc(n)%n211,gldas_struc(n)%n221,&
               gldas_struc(n)%w111,gldas_struc(n)%w121,&
               gldas_struc(n)%w211,gldas_struc(n)%w221)

          allocate(gldas_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gldas_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gldas_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gldas_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gldas_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gldas_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gldas_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gldas_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci,&
               gldas_struc(n)%n112,gldas_struc(n)%n122,&
               gldas_struc(n)%n212,gldas_struc(n)%n222,&
               gldas_struc(n)%w112,gldas_struc(n)%w122,&
               gldas_struc(n)%w212,gldas_struc(n)%w222)

          allocate(gldas_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%colfill(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%rowfill(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          gldas_struc(n)%fillflag1 = .true. 

       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
          allocate(gldas_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci,&
               gldas_struc(n)%n113)

          allocate(gldas_struc(n)%smask1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%colfill(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(gldas_struc(n)%rowfill(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          gldas_struc(n)%fillflag1 = .true. 
       endif
    enddo
  end subroutine init_gldas
end module gldas_forcingMod

