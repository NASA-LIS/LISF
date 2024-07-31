!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module WRF_AKdom_forcingMod
!BOP
! !MODULE: WRF_AKdom_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data extracted from WRF output
!  files. 
!
!  The implementation in LIS has the derived data type {\tt WRFAK\_struc} that
!  includes the variables that specify the runtime options.
!  They are described below: 
!  \begin{description}
!  \item[nest\_id]
!    Value of the WRF nest
!  \item[WRFAKdir]
!    Directory containing the input data
!  \item[ts]
!    Frequency in seconds of the forcing data
!  \item[WRFouttime1]
!    The nearest, previous 1 hour instance of the incoming 
!    data (as a real time). 
!  \item[WRFouttime2]
!    The nearest, next 1 hour instance of the incoming 
!    data (as a real time).
!  \item[findtime1, findtime2]
!    boolean flags to indicate which time is to be read for 
!    temporal interpolation.
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
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

! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_WRF_AKdom      !defines the native resolution of 
                                !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: WRFAK_struc
!EOP

  type, public    :: WRFAK_type_dec

     integer      :: nest_id
     real         :: ts
     integer      :: nc, nr
     integer      :: mi
     character(len=LIS_CONST_PATH_LEN) :: WRFAKdir
     real*8       :: WRFouttime1,WRFouttime2
     integer      :: findtime1,findtime2
     integer      :: nIter, st_iterid, en_iterid

     real         :: gridDesci(50)

     integer, allocatable  :: n111(:)
     integer, allocatable  :: n121(:)
     integer, allocatable  :: n211(:)
     integer, allocatable  :: n221(:)
     real, allocatable     :: w111(:),w121(:)
     real, allocatable     :: w211(:),w221(:)
     integer, allocatable  :: n112(:,:)
     integer, allocatable  :: n122(:,:)
     integer, allocatable  :: n212(:,:)
     integer, allocatable  :: n222(:,:)
     real, allocatable     :: w112(:,:),w122(:,:)
     real, allocatable     :: w212(:,:),w222(:,:)

     ! Neighbor
     integer, allocatable  :: n113(:)

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type WRFAK_type_dec

  type(WRFAK_type_dec), allocatable :: WRFAK_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_WRF_AKdom
! \label{init_WRF_AKdom}
!
! !REVISION HISTORY: 
! 21 Jun 2021; K.R. Arsenault, Updated for different WRF output files
! 
! !INTERFACE:
  subroutine init_WRF_AKdom(findex)
! !USES: 
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain, LIS_localPet,&
        LIS_ews_ind, LIS_ewe_ind,&
        LIS_nss_ind, LIS_nse_ind,&
        LIS_ews_halo_ind,LIS_ewe_halo_ind,&
        LIS_nss_halo_ind, LIS_nse_halo_ind
  use LIS_timeMgrMod, only : LIS_update_timestep, LIS_calendar
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, LIS_verify
  use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
  use LIS_forecastMod
  use map_utils  

   implicit none
! !ARGUMENTS:  
   integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for LIS output
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  Note that no interpolation is performed for this forcing. 
!  The data is expected to be in the same map projection
!  and resolution as that of the current LIS run.
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readconfig\_WRF_AKdom](\ref{readconfig_WRF_AKdom}) \newline
!    reads the runtime options specified for WRF output data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    
   integer :: n
   integer :: status
   integer :: doy
   real    :: gmt
   type(ESMF_Time)  :: DatastartTime
   type(ESMF_Time)  :: LISstartTime

   type(proj_info)   :: proj_temp        
   real              :: lat_str,lon_str  
   integer           :: c,r              
   real, allocatable :: xlat(:,:),xlon(:,:) 
    
    allocate(WRFAK_struc(LIS_rc%nnest))

    LIS_rc%met_nf(findex) = 8

    ! Read inputs from LIS config file:
    call readconfig_WRF_AKdom()

    ! Start date of WRF data availability: Sept 1, 2002
    call ESMF_TimeSet(DatastartTime, yy=2002, &
           mm = 9, dd = 1, h=0, m = 0, s=0, calendar=LIS_calendar, &
           rc=status)
    call LIS_verify(status, &
            'WRF_AKdom_forcingMod: ESMF_TimeSet DatastartTime')
    call ESMF_TimeSet(LISstartTime, yy=LIS_rc%syr, &
             mm = LIS_rc%smo, dd = LIS_rc%sda, h=LIS_rc%shr, &
             m = LIS_rc%smn, s=LIS_rc%sss, calendar=LIS_calendar, &
             rc=status)
    call LIS_verify(status,'WRF_AKdom_forcingMod: ESMF_TimeSet LISstartTime')

    ! Check if LIS start date/time begins before Data start date/time :
    if( LISstartTime < DatastartTime ) then
      write(LIS_logunit,*) &
            '[ERR] LIS start date+time is before WRF-AK data start date+time ...'
      write(LIS_logunit,*) ' WRF-AK data is available from Oct. 1, 2000, and onwards.'
      write(LIS_logunit,*) ' Please update your lis.config file ...'
      call LIS_endrun()
    endif

    ! Number of x- and y-dir points:
    WRFAK_struc(:)%nc = 659
    WRFAK_struc(:)%nr = 609

    do n=1,LIS_rc%nnest

       ! WRF data timestep: 1 hour  == 3600 sec (60 min * 60 sec)
       WRFAK_struc(n)%ts = 60*60
       call LIS_update_timestep(LIS_rc, n, WRFAK_struc(n)%ts)

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

          allocate(WRFAK_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(WRFAK_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

          WRFAK_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          WRFAK_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations

       ! Regular retrospective or non-forecast mode:
       else
          allocate(WRFAK_struc(n)%metdata1(1,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(WRFAK_struc(n)%metdata2(1,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

          WRFAK_struc(n)%st_iterid = 1
          WRFAK_struc(n)%en_iterId = 1
          WRFAK_struc(n)%nIter = 1

       endif

       WRFAK_struc(n)%metdata1 = 0
       WRFAK_struc(n)%metdata2 = 0

      ! See below link for more information on polar stereographic grids:
      !  https://apps.ecmwf.int/codes/grib/format/grib2/templates/3/20

       WRFAK_struc(n)%gridDesci = 0
       WRFAK_struc(n)%gridDesci(1) = 5          ! Polar-stereographic
       WRFAK_struc(n)%gridDesci(2) = WRFAK_struc(n)%nc
       WRFAK_struc(n)%gridDesci(3) = WRFAK_struc(n)%nr
       WRFAK_struc(n)%gridDesci(4) = 51.5418    ! latitude of origin -- LL Lat
       WRFAK_struc(n)%gridDesci(5) = -168.175   ! longitude of origin -- LL Lon
       WRFAK_struc(n)%gridDesci(6) = 0          ! ** Matches NAM242 AK grid entry
       WRFAK_struc(n)%gridDesci(7) = -150.0     ! Orientation (LoV)
       WRFAK_struc(n)%gridDesci(8) = 4.0        ! grid spacing in km
       WRFAK_struc(n)%gridDesci(9) = 4.0        ! grid spacing in km
       WRFAK_struc(n)%gridDesci(10) = 64.0      ! true lat1
       WRFAK_struc(n)%gridDesci(11) = -150.0    ! standard long
       WRFAK_struc(n)%gridDesci(13) = 0         ! regional domain
       WRFAK_struc(n)%gridDesci(20) = 0         ! ** Matches NAM242 AK grid entry

       ! CEN_LAT = 64.f 
       ! CEN_LON = -150.f 
       ! TRUELAT1 = 64.f 
       ! TRUELAT2 = 64.f 
       ! MOAD_CEN_LAT = 64.f 
       ! STAND_LON = -150.f 
       ! POLE_LAT = 90.f 
       ! POLE_LON = 0.f 

       WRFAK_struc(n)%mi = WRFAK_struc(n)%nc*WRFAK_struc(n)%nr

     ! Setting up weights for Interpolation
       select case( LIS_rc%met_interp(findex) )

         case( "bilinear" )
          allocate(WRFAK_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,WRFAK_struc(n)%gridDesci,&
               WRFAK_struc(n)%n111,WRFAK_struc(n)%n121,&
               WRFAK_struc(n)%n211,WRFAK_struc(n)%n221,&
               WRFAK_struc(n)%w111,WRFAK_struc(n)%w121,&
               WRFAK_struc(n)%w211,WRFAK_struc(n)%w221)

         case( "budget-bilinear" )
          allocate(WRFAK_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(WRFAK_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,WRFAK_struc(n)%gridDesci,&
               WRFAK_struc(n)%n111,WRFAK_struc(n)%n121,&
               WRFAK_struc(n)%n211,WRFAK_struc(n)%n221,&
               WRFAK_struc(n)%w111,WRFAK_struc(n)%w121,&
               WRFAK_struc(n)%w211,WRFAK_struc(n)%w221)

          allocate(WRFAK_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(WRFAK_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,WRFAK_struc(n)%gridDesci,&
               WRFAK_struc(n)%n112,WRFAK_struc(n)%n122,&
               WRFAK_struc(n)%n212,WRFAK_struc(n)%n222,&
               WRFAK_struc(n)%w112,WRFAK_struc(n)%w122,&
               WRFAK_struc(n)%w212,WRFAK_struc(n)%w222)

         case( "neighbor" )

          allocate(WRFAK_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,WRFAK_struc(n)%gridDesci,&
                        WRFAK_struc(n)%n113)

       case default
         write(LIS_logunit,*) "[ERR] User-input issue with WRF AK domain forcing ..."
         write(LIS_logunit,*) " -- Currently only supported interpolation options include:"
         write(LIS_logunit,*) "  - bilinear, budget-bilinear, or neighbor - "
         call LIS_endrun
       end select

       ! Read in WRF-forcing terrain height file
       if ( LIS_rc%met_ecor(findex) .ne. "none" ) then 
          call read_WRFAK_elev(n,findex)
       endif

       ! Set up precipitation climate downscaling:
       if(LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,&
               WRFAK_struc(n)%nc,&
               WRFAK_struc(n)%nr)
       endif

    enddo  ! End nest loop

  end subroutine init_WRF_AKdom

end module WRF_AKdom_forcingMod
