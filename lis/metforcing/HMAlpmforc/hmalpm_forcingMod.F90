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
module hmalpm_forcingMod
!BOP
! !MODULE: hmalpm_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the HMALPM forcing data.
!  The data is over HMA 0.05 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt hmalpm\_struc}
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
!  \item[hmalpmtime1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[hmalpmtime2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[hmalpmdir]
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
  public :: init_hmalpm      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: hmalpm_struc

!EOP
  type, public ::  hmalpm_type_dec

     integer      :: npts
     real         :: ts
     integer      :: ncold, nrold
     character*100 :: hmalpmdir   !HMALPM Forcing Directory
     character*100 :: mapfile
     real*8       :: hmalpmtime1,hmalpmtime2
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

     real, allocatable      :: tair1(:)
     real, allocatable      :: qair1(:)
     real, allocatable      :: wind1(:)
     real, allocatable      :: ps1(:)
     real, allocatable      :: prec1(:)
     real, allocatable      :: dirswd1(:)
     real, allocatable      :: difswd1(:)
     real, allocatable      :: swd1(:)
     real, allocatable      :: lwd1(:)

     real, allocatable      :: tair2(:)
     real, allocatable      :: qair2(:)
     real, allocatable      :: wind2(:)
     real, allocatable      :: ps2(:)
     real, allocatable      :: prec2(:)
     real, allocatable      :: dirswd2(:)
     real, allocatable      :: difswd2(:)
     real, allocatable      :: swd2(:)
     real, allocatable      :: lwd2(:)

     integer            :: nvars
     integer            :: uselml

     real*8             :: ringtime
     
     integer            :: nIter, st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type hmalpm_type_dec

  type(hmalpm_type_dec), allocatable :: hmalpm_struc(:)

contains

!BOP
!
! !ROUTINE: init_hmalpm
! \label{init_hmalpm}
!
! !REVISION HISTORY:
! 23 Dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:
  subroutine init_hmalpm(findex)

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
!  Defines the native resolution of the input forcing for HMALPM
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_hmalpm](\ref{readcrd_hmalpm}) \newline
!     reads the runtime options specified for HMALPM data
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


    allocate(hmalpm_struc(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       
       hmalpm_struc(n)%ncold = 1020
       hmalpm_struc(n)%nrold = 520
!       hmalpm_struc(n)%npts = 340819
!       hmalpm_struc(n)%mo1 = -1
!       hmalpm_struc(n)%mo2 = -1

         
       allocate(hmalpm_struc(n)%tair1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%qair1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%wind1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%ps1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%prec1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%swd1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%lwd1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

       allocate(hmalpm_struc(n)%tair2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%qair2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%wind2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%ps2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%prec2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%swd2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(hmalpm_struc(n)%lwd2(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

    enddo

    call readcrd_hmalpm()
    LIS_rc%met_nf(findex) = 7

    hmalpm_struc%reset_flag = .false.

    do n=1, LIS_rc%nnest
       hmalpm_struc(n)%ts = 3600  !check
       call LIS_update_timestep(LIS_rc, n, hmalpm_struc(n)%ts)
    enddo

    gridDesci = 0

    do n=1,LIS_rc%nnest
       gridDesci(n,1) = 0
       gridDesci(n,2) = hmalpm_struc(n)%ncold
       gridDesci(n,3) = hmalpm_struc(n)%nrold
       gridDesci(n,4) =20.025 
       gridDesci(n,5) = 60.025
       gridDesci(n,6) = 128
       gridDesci(n,7) = 45.975
       gridDesci(n,8) = 110.975
       gridDesci(n,9) = 0.05
       gridDesci(n,10) = 0.05
       gridDesci(n,20) = 0

       hmalpm_struc(n)%mi = hmalpm_struc(n)%ncold*hmalpm_struc(n)%nrold

       ! Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(hmalpm_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(n,:),&
               hmalpm_struc(n)%n111,hmalpm_struc(n)%n121,&
               hmalpm_struc(n)%n211,hmalpm_struc(n)%n221,&
               hmalpm_struc(n)%w111,hmalpm_struc(n)%w121,&
               hmalpm_struc(n)%w211,hmalpm_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(hmalpm_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(hmalpm_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(n,:),&
               hmalpm_struc(n)%n111,hmalpm_struc(n)%n121,&
               hmalpm_struc(n)%n211,hmalpm_struc(n)%n221,&
               hmalpm_struc(n)%w111,hmalpm_struc(n)%w121,&
               hmalpm_struc(n)%w211,hmalpm_struc(n)%w221)

          allocate(hmalpm_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmalpm_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmalpm_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmalpm_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmalpm_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmalpm_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmalpm_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(hmalpm_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          call conserv_interp_input(n, gridDesci(n,:),&
               hmalpm_struc(n)%n112,hmalpm_struc(n)%n122,&
               hmalpm_struc(n)%n212,hmalpm_struc(n)%n222,&
               hmalpm_struc(n)%w112,hmalpm_struc(n)%w122,&
               hmalpm_struc(n)%w212,hmalpm_struc(n)%w222)

       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(hmalpm_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n, gridDesci(n,:),&
               hmalpm_struc(n)%n113)

       else
          write(LIS_logunit,*) '[ERR] Interpolation option '// &
               trim(LIS_rc%met_interp(findex))//&
               ' for HMALPM forcing is not supported'
          call LIS_endrun()
       endif

       call LIS_registerAlarm("HMALPM forcing alarm",&
            86400.0,86400.0)
       hmalpm_struc(n)%startFlag = .true.
       hmalpm_struc(n)%dayFlag = .true.

       hmalpm_struc(n)%nvars = 8

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

          hmalpm_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          hmalpm_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          hmalpm_struc(n)%nIter = LIS_forecast_struc(1)%niterations
          
          allocate(hmalpm_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(hmalpm_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       ! Regular retrospective or non-forecast mode:
       else

          hmalpm_struc(n)%st_iterid = 1
          hmalpm_struc(n)%en_iterId = 1
          hmalpm_struc(n)%nIter = 1
          
          allocate(hmalpm_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(hmalpm_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       endif

       hmalpm_struc(n)%metdata1 = 0
       hmalpm_struc(n)%metdata2 = 0


    enddo   ! End nest loop
    
    
  end subroutine init_hmalpm
end module hmalpm_forcingMod
