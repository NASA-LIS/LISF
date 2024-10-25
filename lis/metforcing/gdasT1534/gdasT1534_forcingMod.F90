!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdasT1534_forcingMod
!BOP
! !MODULE: gdasT1534_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the Global Data
!  Assimilation System (GDAST1534) developed at the Environmental Modeling
!  Center (EMC) of NOAA/NCEP. GDAST1534 forcing variables are produced
!  on a quadratic gaussian grid. 
!  The implementation in LIS has the derived data type {\tt gdasT1534\_struc} that
!  includes the variables that specify the runtime options, and the 
!  weights and neighbor information to be used for spatial interpolation. 
!  They are described below: 
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the GDAST1534 data
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for 
!   temporal interpolation.
!  \item[gdasT1534dir]
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
!  \end{description}
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDAST1534      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdasT1534_struc
!EOP

  type, public :: gdasT1534_type_dec
     real          :: ts
     integer       :: ncold, nrold   !AWIPS 212 dimensions
     integer       :: nmif
     character(len=LIS_CONST_PATH_LEN) :: gdasT1534dir !GDAST1534 Forcing Directory
     real*8        :: gdasT1534time1,gdasT1534time2
     integer       :: mi
     integer       :: findtime1, findtime2

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
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type gdasT1534_type_dec
  
  type(gdasT1534_type_dec), allocatable :: gdasT1534_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_GDAST1534
!  \label{init_GDAST1534}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_GDAST1534(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_logMod,     only : LIS_logunit, LIS_endrun
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    
    implicit none
! !ARGUMENTS: 
    integer,  intent(in)    :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GDAST1534
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes \ref{interp}. Based on the GDAST1534 data map projection
!  and resolution, this routine sets up the spatial interpolation
!  weights. The dates of the GDAST1534 resolution switches are also 
!  defined in this routine. 
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gdasT1534](\ref{readcrd_gdasT1534}) \newline
!     reads the runtime options specified for GDAST1534 data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n 

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GDAS-T1534 forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(gdasT1534_struc(LIS_rc%nnest))
    
    do n=1, LIS_rc%nnest
       gdasT1534_struc(n)%ts = 21600
       call LIS_update_timestep(LIS_rc, n, gdasT1534_struc(n)%ts)
    enddo

    call readcrd_gdasT1534()
    
    gdasT1534_struc(:)%nmif    = 16 
    LIS_rc%met_nf(findex) = 16 !number of met variables in GDAST1534 forcing
    
    do n=1,LIS_rc%nnest

       allocate(gdasT1534_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gdasT1534_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gdasT1534_struc(n)%metdata1 = 0
       gdasT1534_struc(n)%metdata2 = 0

       gdasT1534_struc(n)%findtime1 = 0 
       gdasT1534_struc(n)%findtime2 = 0 

       ! This grid is good for some time in the 1980's.
       ! Look up the exact dates.
       gridDesci = 0 
       gridDesci(1) = 4
       gridDesci(2) = 3072
       gridDesci(3) = 1536
       gridDesci(4) = 89.91153
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) =  -89.91153
       gridDesci(8) = -0.117187
       gridDesci(9) = 0.117187
       gridDesci(10) = 768.0
       gridDesci(20) = 0
       gdasT1534_struc(n)%ncold = 3072
       gdasT1534_struc(n)%nrold = 1536
       gdasT1534_struc(n)%mi = gdasT1534_struc(n)%ncold*gdasT1534_struc(n)%nrold
       
!Setting up weights for Interpolation
       if(LIS_rc%met_interp(findex).eq."bilinear") then 

          allocate(gdasT1534_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n,gridDesci,&
               gdasT1534_struc(n)%n111,gdasT1534_struc(n)%n121,&
               gdasT1534_struc(n)%n211,gdasT1534_struc(n)%n221,&
               gdasT1534_struc(n)%w111,gdasT1534_struc(n)%w121,&
               gdasT1534_struc(n)%w211,gdasT1534_struc(n)%w221)
       elseif(LIS_rc%met_interp(findex).eq."budget-bilinear") then 
          allocate(gdasT1534_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasT1534_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n,gridDesci,&
               gdasT1534_struc(n)%n111,gdasT1534_struc(n)%n121,&
               gdasT1534_struc(n)%n211,gdasT1534_struc(n)%n221,&
               gdasT1534_struc(n)%w111,gdasT1534_struc(n)%w121,&
               gdasT1534_struc(n)%w211,gdasT1534_struc(n)%w221)
          allocate(gdasT1534_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasT1534_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasT1534_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasT1534_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasT1534_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasT1534_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasT1534_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gdasT1534_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          call conserv_interp_input(n,gridDesci,&
               gdasT1534_struc(n)%n112,gdasT1534_struc(n)%n122,&
               gdasT1534_struc(n)%n212,gdasT1534_struc(n)%n222,&
               gdasT1534_struc(n)%w112,gdasT1534_struc(n)%w122,&
               gdasT1534_struc(n)%w212,gdasT1534_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
          allocate(gdasT1534_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n,gridDesci,&
               gdasT1534_struc(n)%n113)

       else
          write(LIS_logunit,*) 'The specified spatial interpolation option not'
          write(LIS_logunit,*) 'supported for GDAST1534..'
          write(LIS_logunit,*) 'program stopping..'
          call LIS_endrun()
       endif
    enddo
  end subroutine init_GDAST1534
end module gdasT1534_forcingMod

