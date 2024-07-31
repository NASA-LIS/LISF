!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module RFE2Daily_forcingMod
!BOP
! !MODULE: RFE2Daily_forcingMod
!
! !DESCRIPTION:
!
!  This module contains variables, data structures and subroutines used for 
!  the implementation of the RFE2.0 rainfall forcing data from the CPC site 
!  ftp://ftp.cpc.ncep.noaa.gov/fews/newalgo\_est/ and used in USAID/FEWS-NET. 
!  This precip is a daily precipitation analysis at 0.1 deg lat x 0.1 deg lon
!  produced by merging GTS gauge observations and  3 kinds of satellite 
!  estimates (GPI,SSM/I and AMSU). Units are in millimeters (mm).
!  Data is in big endian binary with coverage -40.00S to 40.00N Northward
!  (801 grid points in south - north direction)  and 20.00W to 55.00E 
!  Eastward  (751 grid points in east - west direction)  
!
!  The implementation in LIS has the derived data type {\tt RFE2Daily\_struc}
!  that includes the variables to specify the runtime options, and
!  the weights and neighbor information for upscaling or spatial interpolation
!
!  They are desribed below:
!  \begin{description}
!  \item[RFE2DailyDir]
!    Directory containing the input data
!  \item[RFE2DailyEndTime]
!    The nearest, daily instance of the incoming
!    data (as a real time).
!  \item[startTime]
!    The earliest possible daily availability time of 
!    the incoming data (as an ESMF real time).
!  \item[st\_real]
!    The earliest possible daily availability time of 
!    the incoming data (as a real time).
!  \item[timeStep]
!    Input data time step (as an ESMF real time).
!  \item[gridDesci]
!    Input grid description parameters
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
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for n. neighbor interpolation.
!  \end{description}
!
! !USES:
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_RFE2Daily      !defines the native resolution of
                                       !the input data

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: RFE2Daily_struc
!EOP

  type, public :: RFE2Daily_type_dec
     real                    :: ts
     character(len=LIS_CONST_PATH_LEN) :: RFE2DailyDir
     real*8                  :: RFE2DailyEndTime
     type(ESMF_Time)         :: startTime
     real*8                  :: st_real
     type(ESMF_TimeInterval) :: timeStep
     real                    :: gridDesci(50)
     integer                 :: mi
     integer                 :: hour_offset

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

     integer           :: nIter, st_iterid,en_iterid  ! Forecast parameters

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type RFE2Daily_type_dec

  type(RFE2Daily_type_dec), allocatable :: RFE2Daily_struc(:)

contains

!BOP
!
! !ROUTINE: init_RFE2Daily
!  \label{init_RFE2Daily}
!
! !REVISION HISTORY:
!  26 MAY 2010: Soni Yatheendradas; Initial LIS version for FEWSNET
!  20 Mar 2013; KR Arsenault, Cleaned up code and documentation
!
! !INTERFACE:
  subroutine init_RFE2Daily(findex)
! !USES:
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_logMod,     only : LIS_logunit, LIS_endrun
   use LIS_timeMgrMod, only : LIS_update_timestep
   use LIS_FORC_AttributesMod
   use LIS_forecastMod

   implicit none
   integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for RFE2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  and upscaling schemes \ref{interp}
!
!  This routine performs required initializations, including 
!  allocating memory, initializing data structures and setting
!  up interpolation weights and upscaling indices.
!
!  The arguments are:
!  \begin{description}
!  \item[findex]
!    index of the supplemental forcing source
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_RFE2Daily](\ref{readcrd_RFE2Daily}) \newline
!     reads the runtime options specified for RFE2Daily data
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!     computes the neighbor information for upscaling data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!     computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!     computes the neighbor, weights for conservative interpolation
!   \item[neighbor\_interp\_input](\ref{neighbor_interp_input}) \newline
!     computes the neighbor for neighbor interpolation
!  \end{description}
!
!EOP
   
    integer              :: n
 
    allocate(RFE2Daily_struc(LIS_rc%nnest))

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

 !- Read in daily RFE2 lis.config inputs:
    call readcrd_RFE2Daily()
    LIS_rc%met_nf(findex) = 2  ! number of met variables in RFE2 forcing

    do n=1, LIS_rc%nnest

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

         RFE2Daily_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
         RFE2Daily_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
         RFE2Daily_struc(n)%nIter = LIS_forecast_struc(1)%niterations

         allocate(RFE2Daily_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(RFE2Daily_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))

       ! Regular retrospective or non-forecast mode:
       else
         RFE2Daily_struc(n)%st_iterid = 1
         RFE2Daily_struc(n)%en_iterId = 1
         RFE2Daily_struc(n)%nIter = 1

         allocate(RFE2Daily_struc(n)%metdata1(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
         allocate(RFE2Daily_struc(n)%metdata2(1,&
              LIS_rc%met_nf(findex),&
              LIS_rc%ngrid(n)))
       endif

       RFE2Daily_struc(n)%metdata1 = 0
       RFE2Daily_struc(n)%metdata2 = 0

       RFE2Daily_struc(n)%ts = 24*60*60 
       call LIS_update_timestep(LIS_rc, n, RFE2Daily_struc(n)%ts)
    enddo

    do n=1,LIS_rc%nnest

       RFE2Daily_struc(n)%gridDesci = 0

     ! Define parameters for geographic projection:
       RFE2Daily_struc(n)%gridDesci(1) = 0      ! lat/lon projection
       RFE2Daily_struc(n)%gridDesci(4) = -40.05 ! lat of lower left grid-cell
                                                ! (upper left corner)
       RFE2Daily_struc(n)%gridDesci(5) = -19.95 ! long of lower left grid-cell
                                                ! (upper left corner)
       RFE2Daily_struc(n)%gridDesci(7) = 39.95  ! lat of upper right grid-cell
                                                ! (upper left corner)
       RFE2Daily_struc(n)%gridDesci(8) = 55.05  ! long of upper right grid-cell
                                                ! (upper left corner)
       RFE2Daily_struc(n)%gridDesci(9) = 0.1    ! spatial res (in degrees)
                                                ! along E-W
       RFE2Daily_struc(n)%gridDesci(10)= 0.1    ! spatial res (in degrees)
                                                ! along N-S
     ! Number of columns in the domain: 751
       RFE2Daily_struc(n)%gridDesci(2) = NINT((RFE2Daily_struc(n)%gridDesci(8)-&
            RFE2Daily_struc(n)%gridDesci(5))/RFE2Daily_struc(n)%gridDesci(9))+1 
     ! Number of rows in the domain: 801
       RFE2Daily_struc(n)%gridDesci(3) = NINT((RFE2Daily_struc(n)%gridDesci(7)-&
            RFE2Daily_struc(n)%gridDesci(4))/RFE2Daily_struc(n)%gridDesci(10))+1 
       RFE2Daily_struc(n)%gridDesci(6) = 128
       RFE2Daily_struc(n)%gridDesci(20)= 64 
     ! Used to specify the ordering of data 
     !  (non-divisible by 32 indicates E-W ordering, else N-S ordering) ? 
     ! Set now to 64 and not 255 of JimG's initial interp/get_fieldpos explanation: 
     !  E-W ordering have elements are located in array row-wise, i.e., one row, then next and so on. 
 
       write(LIS_logunit,*) 'Number of RFE2Daily grid columns = ', &
                            NINT(RFE2Daily_struc(n)%gridDesci(2))
       write(LIS_logunit,*) 'Number of RFE2Daily grid rows = ', &
                            NINT(RFE2Daily_struc(n)%gridDesci(3))
   
       RFE2Daily_struc(n)%mi = NINT(RFE2Daily_struc(n)%gridDesci(2))* &
                               NINT(RFE2Daily_struc(n)%gridDesci(3))

       write(LIS_logunit,*) "RFE2Daily reprojection choice selected is ",&
                            trim(LIS_rc%met_interp(findex))," (interpolation)"

    !- Determine if upscaling or downscaling option:
       select case( LIS_rc%met_interp(findex) )
   
         case( "average" )   ! Upscaling 

           IF( NINT(LIS_rc%gridDesc(n,1)) .NE. 0) THEN ! SY
             write(LIS_logunit,*) "RFE2Daily precip implementation supports"
             write(LIS_logunit,*) "upscaling only to lat/lon run domain grids ... "
             write(LIS_logunit,*) "Program stopping ..."
             call LIS_endrun()
           ENDIF
           IF( RFE2Daily_struc(n)%gridDesci(10) .GT. LIS_rc%gridDesc(n,10)  ) THEN
         ! Confirm if LIS run domain resolution really not less than RFE2Daily forcing
         ! resolution for lat/lon projection 
             write(LIS_logunit,*) "LIS lat/lon run domain resolution less than "
             write(LIS_logunit,*) "RFE2Daily forcing resolution in upscaling option."
             write(LIS_logunit,*) "Program stopping ..."
             call LIS_endrun()
           ENDIF

           allocate(RFE2Daily_struc(n)%n111(RFE2Daily_struc(n)%mi))
           call upscaleByAveraging_input(RFE2Daily_struc(n)%gridDesci(:),&
                LIS_rc%gridDesc(n,:), RFE2Daily_struc(n)%mi,&
                LIS_rc%lnc(n)*LIS_rc%lnr(n), RFE2Daily_struc(n)%n111)

         case( "bilinear" )
            allocate(RFE2Daily_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

            call bilinear_interp_input(n,RFE2Daily_struc(n)%gridDesci(:), &
                 RFE2Daily_struc(n)%n111,RFE2Daily_struc(n)%n121, &
                 RFE2Daily_struc(n)%n211,&
                 RFE2Daily_struc(n)%n221,RFE2Daily_struc(n)%w111, &
                 RFE2Daily_struc(n)%w121,&
                 RFE2Daily_struc(n)%w211,RFE2Daily_struc(n)%w221)

         case( "budget-bilinear" )
            allocate(RFE2Daily_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

            call conserv_interp_input(n,RFE2Daily_struc(n)%gridDesci(:), &
                 RFE2Daily_struc(n)%n112,RFE2Daily_struc(n)%n122,&
                 RFE2Daily_struc(n)%n212,RFE2Daily_struc(n)%n222,&
                 RFE2Daily_struc(n)%w112,RFE2Daily_struc(n)%w122,&
                 RFE2Daily_struc(n)%w212,RFE2Daily_struc(n)%w222)

         case( "neighbor" )
            allocate(RFE2Daily_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

            call neighbor_interp_input(n,RFE2Daily_struc(n)%gridDesci(:), &
                 RFE2Daily_struc(n)%n113)

         case default
           write(LIS_logunit,*) "This interp. option not supported for RFE2Daily Precip."
           write(LIS_logunit,*) "Program stopping ... "
           call LIS_endrun()

       end select ! end interp option

     enddo !  end nest loop
  
  end subroutine init_RFE2Daily

end module RFE2Daily_forcingMod

