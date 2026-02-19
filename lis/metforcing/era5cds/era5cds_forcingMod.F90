!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module era5cds_forcingMod
!BOP
! !MODULE: era5cds_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the ERA5 forcing data from Climate Data Store.
!  The data is global 0.25 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt era5cds\_struc}
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
!  \item[era5cdstime1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[era5cdstime2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[era5cdsdir]
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

  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_era5cds      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: era5cds_struc

!EOP
  type, public ::  era5cds_type_dec

     real         :: ts
     integer      :: ncold, nrold
     character(len=LIS_CONST_PATH_LEN) :: era5cdsdir   !ERA5CDS Forcing Directory
     character(len=LIS_CONST_PATH_LEN) :: era5cdsalt_file ! ERA5CDS altitude file
     real*8       :: era5cdstime1,era5cdstime2
     character*50 :: met_interp
     logical      :: reset_flag
     integer      :: mon

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
   ! Subset parameters for "none" reprojection type:
     integer, allocatable  :: lat_line(:,:)
     integer, allocatable  :: lon_line(:,:)
     real                  :: subset_gridDesc(20)
     integer               :: subset_nc, subset_nr

     integer                :: findtime1, findtime2
     logical                :: startFlag, dayFlag
     real*8                 :: validstart

     real, allocatable      :: tair(:,:)
     real, allocatable      :: qair(:,:)
     real, allocatable      :: uwind(:,:)
     real, allocatable      :: vwind(:,:)
     real, allocatable      :: ps(:,:)
     real, allocatable      :: rainf(:,:)
     real, allocatable      :: crainf(:,:)
     real, allocatable      :: swd(:,:)
     real, allocatable      :: lwd(:,:)

     real, allocatable      :: prev_rainf(:,:)
     real, allocatable      :: prev_crainf(:,:)
     real, allocatable      :: prev_swd(:,:)
     real, allocatable      :: prev_lwd(:,:)

     integer            :: nvars
     integer            :: uselml

     real*8             :: ringtime
     
     integer            :: nIter, st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type era5cds_type_dec

  type(era5cds_type_dec), allocatable :: era5cds_struc(:)

contains

!BOP
!
! !ROUTINE: init_era5cds
! \label{init_era5cds}
!
! !REVISION HISTORY:
! 23 Dec 2019: Sujay Kumar, initial code 
! 04 Mar 2025: Hiroko Beudoing, adopted ERA5 routines for the public CDS
!                               data format
!
! !INTERFACE:
  subroutine init_era5cds(findex)

! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use LIS_forecastMod
    use LIS_gridmappingMod, only : LIS_RunDomainPts
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for ERA5CDS
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_era5cds](\ref{readcrd_era5cds}) \newline
!     reads the runtime options specified for ERA5CDS data
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
    integer      :: glpnc, glpnr


    allocate(era5cds_struc(LIS_rc%nnest))

    do n=1, LIS_rc%nnest

       era5cds_struc(n)%ncold = 1440
       era5cds_struc(n)%nrold = 720
       era5cds_struc(n)%mon = -1

       allocate(era5cds_struc(n)%tair(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(era5cds_struc(n)%qair(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(era5cds_struc(n)%uwind(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(era5cds_struc(n)%vwind(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(era5cds_struc(n)%ps(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(era5cds_struc(n)%rainf(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(era5cds_struc(n)%crainf(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(era5cds_struc(n)%swd(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(era5cds_struc(n)%lwd(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))

       ! save accumulation fields
       ! previous month's 17-23Z on last day = 0-6Z on 1st next month
       allocate(era5cds_struc(n)%prev_rainf(LIS_rc%lnc(n)*LIS_rc%lnr(n),7))
       allocate(era5cds_struc(n)%prev_crainf(LIS_rc%lnc(n)*LIS_rc%lnr(n),7))
       allocate(era5cds_struc(n)%prev_swd(LIS_rc%lnc(n)*LIS_rc%lnr(n),7))
       allocate(era5cds_struc(n)%prev_lwd(LIS_rc%lnc(n)*LIS_rc%lnr(n),7))

    enddo

    call readcrd_era5cds()
    LIS_rc%met_nf(findex) = 9

    era5cds_struc%reset_flag = .false.

    do n=1, LIS_rc%nnest
       era5cds_struc(n)%ts = 3600  !hourly
       call LIS_update_timestep(LIS_rc, n, era5cds_struc(n)%ts)
    enddo

    gridDesci = 0
    LIS_rc%met_proj(findex) = "latlon"

    do n=1,LIS_rc%nnest
       gridDesci(n,1) = 0
       gridDesci(n,2) = era5cds_struc(n)%ncold
       gridDesci(n,3) = era5cds_struc(n)%nrold
       gridDesci(n,4) = -89.875
       gridDesci(n,5) = -179.875
       gridDesci(n,6) = 128
       gridDesci(n,7) = 89.875
       gridDesci(n,8) = 179.875
       gridDesci(n,9) = 0.25
       gridDesci(n,10) = 0.25
       gridDesci(n,20) = 0

       ! ERA5 accumulation data starts at 7z on 1940-01-01
       yr1 = 1940
       mo1 = 01
       da1 = 01
       hr1 = 07
       mn1 = 0; ss1 = 0
       call LIS_date2time( era5cds_struc(n)%validstart,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       era5cds_struc(n)%mi = era5cds_struc(n)%ncold*era5cds_struc(n)%nrold

       ! Check resolution and set up weights for Interpolation
       if ( LIS_isatAfinerResolution(n,gridDesci(n,9)) )  then

         era5cds_struc(n)%met_interp = LIS_rc%met_interp(findex)

         write(LIS_logunit,*) 'MSG: The ERA5CDS forcing resolution is ' // &
                              ' coaser than the running domain.'
         write(LIS_logunit,*) '     Interpolating with the ' // &
                               trim(era5cds_struc(n)%met_interp) // ' method.'
       ! Setting up weights for Interpolation
         if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(era5cds_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(n,:),&
               era5cds_struc(n)%n111,era5cds_struc(n)%n121,&
               era5cds_struc(n)%n211,era5cds_struc(n)%n221,&
               era5cds_struc(n)%w111,era5cds_struc(n)%w121,&
               era5cds_struc(n)%w211,era5cds_struc(n)%w221)

         elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(era5cds_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(era5cds_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(n,:),&
               era5cds_struc(n)%n111,era5cds_struc(n)%n121,&
               era5cds_struc(n)%n211,era5cds_struc(n)%n221,&
               era5cds_struc(n)%w111,era5cds_struc(n)%w121,&
               era5cds_struc(n)%w211,era5cds_struc(n)%w221)

          allocate(era5cds_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(era5cds_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(era5cds_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(era5cds_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(era5cds_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(era5cds_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(era5cds_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(era5cds_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          call conserv_interp_input(n, gridDesci(n,:),&
               era5cds_struc(n)%n112,era5cds_struc(n)%n122,&
               era5cds_struc(n)%n212,era5cds_struc(n)%n222,&
               era5cds_struc(n)%w112,era5cds_struc(n)%w122,&
               era5cds_struc(n)%w212,era5cds_struc(n)%w222)

         elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(era5cds_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n, gridDesci(n,:),&
               era5cds_struc(n)%n113)

         endif
       ! Running domain at 0.25 degree, no need for interpolation
       else if (gridDesci(n,9)  == LIS_rc%gridDesc(n,9) .and. &
                gridDesci(n,10) == LIS_rc%gridDesc(n,10).and. &
                LIS_rc%gridDesc(n,1) == 0 ) then
          write(LIS_logunit,*) '[INFO] ERA5CDS and LIS resolutions match ' // &
                               'no spatial trainsorm of input forcing (only subset may occur). '
          era5cds_struc(n)%met_interp = "none"
          call LIS_RunDomainPts( n, LIS_rc%met_proj(findex), gridDesci(n,:), &
               glpnc, glpnr, era5cds_struc(n)%subset_nc,                     &
               era5cds_struc(n)%subset_nr, era5cds_struc(n)%subset_gridDesc, &
               era5cds_struc(n)%lat_line, era5cds_struc(n)%lon_line )

       else
          era5cds_struc(n)%met_interp = LIS_rc%met_upscale(findex)
          write(LIS_logunit,*) 'MSG: The ERA5CDS forcing resolution is finer ' // &
                               'than the running domain.'
          write(LIS_logunit,*) '     Upscaling with the ' // &
                               trim(era5cds_struc(n)%met_interp) // ' method.'

          select case( era5cds_struc(n)%met_interp )
           case( "average" )
            allocate(era5cds_struc(n)%n111(era5cds_struc(n)%mi))

            call upscaleByAveraging_input(gridDesci,                   &
                                          LIS_rc%gridDesc(n,:),        &
                                          era5cds_struc(n)%mi,           &
                                          LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                                          era5cds_struc(n)%n111)
           case default
            write(LIS_logunit,*) '[ERR] Interpolation option '// &
                 trim(LIS_rc%met_interp(findex))//&
                 ' for ERA5CDS forcing is not supported'
            call LIS_endrun()
          end select
       endif

       call LIS_registerAlarm("ERA5CDS forcing alarm",&
            86400.0,86400.0)
       era5cds_struc(n)%startFlag = .true.
       era5cds_struc(n)%dayFlag = .true.

       era5cds_struc(n)%nvars = 9

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

          era5cds_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          era5cds_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          era5cds_struc(n)%nIter = LIS_forecast_struc(1)%niterations
          
          allocate(era5cds_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(era5cds_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       ! Regular retrospective or non-forecast mode:
       else

          era5cds_struc(n)%st_iterid = 1
          era5cds_struc(n)%en_iterId = 1
          era5cds_struc(n)%nIter = 1
          
          allocate(era5cds_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(era5cds_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       endif

       era5cds_struc(n)%metdata1 = 0
       era5cds_struc(n)%metdata2 = 0

       era5cds_struc(n)%findtime1 = 0
       era5cds_struc(n)%findtime2 = 0

       ! Set up precipitation climate downscaling:
       if(LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,&
               era5cds_struc(n)%ncold,&
               era5cds_struc(n)%nrold)
       endif
       
       if ( LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
            LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" ) then

          call read_era5cds_elev(n,findex)
       endif

    enddo   ! End nest loop
    
    
  end subroutine init_era5cds
end module era5cds_forcingMod

