!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module hmaens_forcingMod
!BOP
! !MODULE: hmaens_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the HMAENS forcing data.
!  The data is over HMA 0.05 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt hmaens\_struc}
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
!  \item[hmaenstime1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[hmaenstime2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[hmaensdir]
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
  public :: init_hmaens      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: hmaens_struc

!EOP
  type, public ::  hmaens_type_dec

     integer      :: npts
     real         :: ts
     integer      :: ncold, nrold
     character*40 :: hmaensdir   !HMAENS Forcing Directory
     character*40 :: mapfile
     real*8       :: hmaenstime1,hmaenstime2
     logical      :: reset_flag
     integer      :: mo1,mo2

     integer, allocatable :: G2P(:,:)
     integer                :: mi
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
     integer                :: findtime1, findtime2
     logical                :: startFlag, dayFlag

     real, allocatable      :: tair1(:,:)
     real, allocatable      :: qair1(:,:)
     real, allocatable      :: wind1(:,:)
     real, allocatable      :: ps1(:,:)
     real, allocatable      :: prec1(:,:)
     real, allocatable      :: dirswd1(:,:)
     real, allocatable      :: difswd1(:,:)
     real, allocatable      :: swd1(:,:)
     real, allocatable      :: lwd1(:,:)

     real, allocatable      :: tair2(:,:)
     real, allocatable      :: qair2(:,:)
     real, allocatable      :: wind2(:,:)
     real, allocatable      :: ps2(:,:)
     real, allocatable      :: prec2(:,:)
     real, allocatable      :: dirswd2(:,:)
     real, allocatable      :: difswd2(:,:)
     real, allocatable      :: swd2(:,:)
     real, allocatable      :: lwd2(:,:)

     integer            :: nvars
     integer            :: uselml

     real*8             :: ringtime
     
     integer            :: nIter, st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type hmaens_type_dec

  type(hmaens_type_dec), allocatable :: hmaens_struc(:)

contains

!BOP
!
! !ROUTINE: init_hmaens
! \label{init_hmaens}
!
! !REVISION HISTORY:
! 23 Dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:
  subroutine init_hmaens(findex)

! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use LIS_forecastMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for HMAENS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_hmaens](\ref{readcrd_hmaens}) \newline
!     reads the runtime options specified for HMAENS data
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
    integer :: ftn
!    integer :: G2Pid


    allocate(hmaens_struc(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       
       hmaens_struc(n)%ncold = 1020
       hmaens_struc(n)%nrold = 520
!       hmaens_struc(n)%npts = 340819
!       hmaens_struc(n)%mo1 = -1
!       hmaens_struc(n)%mo2 = -1

         
       allocate(hmaens_struc(n)%tair1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%qair1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%wind1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%ps1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%prec1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%swd1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%lwd1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))

       allocate(hmaens_struc(n)%tair2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%qair2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%wind2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%ps2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%prec2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%swd2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(hmaens_struc(n)%lwd2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))

    enddo

    call readcrd_hmaens()
    LIS_rc%met_nf(findex) = 7

    hmaens_struc%reset_flag = .false.

    do n=1, LIS_rc%nnest
       hmaens_struc(n)%ts = 3600  !check
       call LIS_update_timestep(LIS_rc, n, hmaens_struc(n)%ts)
    enddo

    gridDesci = 0

    do n=1,LIS_rc%nnest
       gridDesci(n,1) = 0
       gridDesci(n,2) = hmaens_struc(n)%ncold
       gridDesci(n,3) = hmaens_struc(n)%nrold
       gridDesci(n,4) =20.025 
       gridDesci(n,5) = 69.025
       gridDesci(n,6) = 128
       gridDesci(n,7) = 45.975
       gridDesci(n,8) = 110.975
       gridDesci(n,9) = 0.05
       gridDesci(n,10) = 0.05
       gridDesci(n,20) = 0

       hmaens_struc(n)%mi = hmaens_struc(n)%ncold*hmaens_struc(n)%nrold

       ! Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(hmaens_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(n,:),&
               hmaens_struc(n)%n111,hmaens_struc(n)%n121,&
               hmaens_struc(n)%n211,hmaens_struc(n)%n221,&
               hmaens_struc(n)%w111,hmaens_struc(n)%w121,&
               hmaens_struc(n)%w211,hmaens_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(hmaens_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmaens_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(n,:),&
               hmaens_struc(n)%n111,hmaens_struc(n)%n121,&
               hmaens_struc(n)%n211,hmaens_struc(n)%n221,&
               hmaens_struc(n)%w111,hmaens_struc(n)%w121,&
               hmaens_struc(n)%w211,hmaens_struc(n)%w221)

          allocate(hmaens_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmaens_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmaens_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmaens_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmaens_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmaens_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmaens_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmaens_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          call conserv_interp_input(n, gridDesci(n,:),&
               hmaens_struc(n)%n112,hmaens_struc(n)%n122,&
               hmaens_struc(n)%n212,hmaens_struc(n)%n222,&
               hmaens_struc(n)%w112,hmaens_struc(n)%w122,&
               hmaens_struc(n)%w212,hmaens_struc(n)%w222)

       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(hmaens_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n, gridDesci(n,:),&
               hmaens_struc(n)%n113)

       else
          write(LIS_logunit,*) '[ERR] Interpolation option '// &
               trim(LIS_rc%met_interp(findex))//&
               ' for HMAENS forcing is not supported'
          call LIS_endrun()
       endif

       call LIS_registerAlarm("HMAENS forcing alarm",&
            86400.0,86400.0)
       hmaens_struc(n)%startFlag = .true.
       hmaens_struc(n)%dayFlag = .true.

       hmaens_struc(n)%nvars = 8

       ! Forecast mode:
       if(LIS_rc%forecastMode.eq.1) then 
          
          if(mod(LIS_rc%nensem(n),&
               LIS_forecast_struc(1)%niterations).ne.0) then 
             write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
             write(LIS_logunit,*) '[ERR] of the number of iterations '
             write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
             write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
             call LIS_endrun()
          endif

          hmaens_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          hmaens_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          hmaens_struc(n)%nIter = LIS_forecast_struc(1)%niterations
          
          allocate(hmaens_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(hmaens_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       ! Regular retrospective or non-forecast mode:
       else

          hmaens_struc(n)%st_iterid = 1
          hmaens_struc(n)%en_iterId = 1
          hmaens_struc(n)%nIter = 1
          
          allocate(hmaens_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(hmaens_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       endif

       hmaens_struc(n)%metdata1 = 0
       hmaens_struc(n)%metdata2 = 0


    enddo   ! End nest loop
    
    
  end subroutine init_hmaens
end module hmaens_forcingMod
