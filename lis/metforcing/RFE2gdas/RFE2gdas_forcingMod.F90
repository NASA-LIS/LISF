!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module RFE2gdas_forcingMod
!BOP
! !MODULE: RFE2gdas_forcingMod
!
! !DESCRIPTION:
!
!  This module contains variables, data structures and subroutines used for 
!  the implementation of the RFE2gdas.0 rainfall forcing data from the CPC site 
!  ftp://ftp.cpc.ncep.noaa.gov/fews/newalgo\_est/ and used in USAID/FEWS-NET. 
!  This precip is a daily precipitation analysis at 0.1 deg lat x 0.1 deg lon
!  produced by merging GTS gauge observations and  3 kinds of satellite 
!  estimates (GPI,SSM/I and AMSU). Units are in millimeters (mm).
!  Data is in big endian binary with coverage -40.00S to 40.00N Northward
!  (801 grid points in south - north direction)  and 20.00W to 55.00E 
!  Eastward  (751 grid points in east - west direction)  
!
!  The implementation in LIS has the derived data type {\tt RFE2gdas\_struc}
!  that includes the variables to specify the runtime options, and
!  the weights and neighbor information for upscaling or spatial interpolation
!
!  They are desribed below:
!  \begin{description}
!  \item[RFE2gdasDir]
!    Directory containing the input data
!  \item[RFE2gdasEndTime]
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
!  \item[stc]
!    Starting index (along the east-west direction)
!    of the input grid for each output grid point
!  \item[str]
!    Starting index (along the north-south direction)
!    of the input grid for each output grid point
!  \item[enc]
!    Ending index (along the east-west direction)
!    of the input grid for each output grid point
!  \item[enr]
!    Ending index (along the north-south direction)
!    of the input grid for each output grid point
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
  public :: init_RFE2gdas      ! Defines the native resolution of
                               ! the input data

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: RFE2gdas_struc
!EOP

  type, public :: RFE2gdas_type_dec

     character(len=LIS_CONST_PATH_LEN) :: RFE2gdasDir
     real*8                   :: RFE2gdasEndTime
     type(ESMF_Time)          :: startTime
     real*8                   :: st_real
     type(ESMF_TimeInterval)  :: timeStep
     real                     :: gridDesci(50)
     integer                  :: mi
     real                     :: ts
     ! Upscaling arrays
     integer, allocatable :: stc(:,:)
     integer, allocatable :: str(:,:)
     integer, allocatable :: enc(:,:)
     integer, allocatable :: enr(:,:)

     ! Interpolation arrays

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

     integer            :: nIter, st_iterid,en_iterid
     real, allocatable      :: metdata1(:,:,:) 
     real, allocatable      :: metdata2(:,:,:) 

  end type RFE2gdas_type_dec

  type(RFE2gdas_type_dec), allocatable :: RFE2gdas_struc(:)

contains

!BOP
!
! !ROUTINE: init_RFE2gdas
!  \label{init_RFE2gdas}
!
! !REVISION HISTORY:
! 26 MAY 2010: Soni Yatheendradas; Initial LIS version for FEWSNET
!
! !INTERFACE:
  subroutine init_RFE2gdas(findex)
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use LIS_logMod,     only : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod
    use LIS_forecastMod

    implicit none

    integer,  intent(in)     :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for CMAP
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  and upscaling schemes (see Section~\ref{interp}).
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
!   \item[readRFE2gdascrd](\ref{readRFE2gdascrd}) \newline
!     reads the runtime options specified for RFE2gdas data
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbor information for upscaling data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[neighbor\_interp\_input](\ref{neighbor_interp_input}) \newline
!    computes the neighbor for neighbor interpolation
!  \end{description}
!
!EOP
   
    integer    :: n
 
    allocate(RFE2gdas_struc(LIS_rc%nnest))

    write(LIS_logunit,fmt=*)"[INFO] Initializing RFE2-GDAS forcing grid ... "

    ! Temporary note to alert users of issue with convective precip ratios:
    if( LIS_FORC_CRainf%selectOpt == 1 ) then
      write(LIS_logunit,*)"[WARN] At this time, convective rainfall is NOT constrained"
      write(LIS_logunit,*)"[WARN]  to match this supplemental observed rainfall dataset."
      write(LIS_logunit,*)" -- This feature will be applied in future LIS releases -- "
    endif

    call readRFE2gdascrd()

    do n=1, LIS_rc%nnest
       RFE2gdas_struc(n)%ts = 21600
       call LIS_update_timestep(LIS_rc, n, RFE2gdas_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 2

    do n=1,LIS_rc%nnest

       if(LIS_rc%forecastMode.eq.1) then 
          
          if(mod(LIS_rc%nensem(n),&
               LIS_forecast_struc(1)%niterations).ne.0) then 
             write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
             write(LIS_logunit,*) '[ERR] of the number of iterations '
             write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
             write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
             call LIS_endrun()
          endif
          allocate(RFE2gdas_struc(n)%metdata1(&
               LIS_forecast_struc(1)%niterations,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(RFE2gdas_struc(n)%metdata2(&
               LIS_forecast_struc(1)%niterations,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

          RFE2gdas_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          RFE2gdas_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          RFE2gdas_struc(n)%nIter = LIS_forecast_struc(1)%niterations

       else ! Non-forecast mode

          allocate(RFE2gdas_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(RFE2gdas_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))

          RFE2gdas_struc(n)%st_iterid = 1
          RFE2gdas_struc(n)%en_iterId = 1
          RFE2gdas_struc(n)%nIter = 1

       endif

       RFE2gdas_struc(n)%metdata1 = 0
       RFE2gdas_struc(n)%metdata2 = 0

       RFE2gdas_struc(n)%gridDesci = 0

     ! Define parameters for lat/lon projection
       RFE2gdas_struc(n)%gridDesci(1) = 0      ! indicates lat/lon projection 
       RFE2gdas_struc(n)%gridDesci(8) = 55.0   ! longitude of the upper right corner CELL CENTER of the domain?
       RFE2gdas_struc(n)%gridDesci(5) = -20.0  ! longitude of the lower left corner CELL CENTER of the domain?
       RFE2gdas_struc(n)%gridDesci(7) = 40.0   ! latitude of the upper right corner CELL CENTER of the domain?
       RFE2gdas_struc(n)%gridDesci(4) = -40.0  ! latitude of the lower left corner CELL CENTER of the domain?
       RFE2gdas_struc(n)%gridDesci(10) = 0.1   ! spatial resolution (in degrees) along the E-W dimension 
       RFE2gdas_struc(n)%gridDesci(9) = 0.1    ! spatial resolution (in degrees) along the N-S dimension 
       RFE2gdas_struc(n)%gridDesci(2) = NINT((RFE2gdas_struc(n)%gridDesci(8)-&
            RFE2gdas_struc(n)%gridDesci(5))/RFE2gdas_struc(n)%gridDesci(10))+1 
     ! number of columns in the domain: 751  
       RFE2gdas_struc(n)%gridDesci(3) = NINT((RFE2gdas_struc(n)%gridDesci(7)-&
            RFE2gdas_struc(n)%gridDesci(4))/RFE2gdas_struc(n)%gridDesci(9))+1 
     ! number of rows in the domain: 801 
       RFE2gdas_struc(n)%gridDesci(6) = 128 ! not used
       RFE2gdas_struc(n)%gridDesci(20) = 64 ! used to specify the ordering of data 
     ! (Non-divisible by 32 indicates E-W ordering, else N-S ordering) 
     !  Set now to 255 per interp/get_fieldpos explanation:  
     !   E-W ordering indicates elements located in row-order array, i.e., one row, then next and so on. 
 
       write(LIS_logunit,*)'[INFO] Number of RFE2gdas grid columns = ', &
                            NINT(RFE2gdas_struc(n)%gridDesci(2))
       write(LIS_logunit,*)'[INFO] Number of RFE2gdas grid rows = ', &
                            NINT(RFE2gdas_struc(n)%gridDesci(3))
   
       RFE2gdas_struc(n)%mi = NINT(RFE2gdas_struc(n)%gridDesci(2))* &
                              NINT(RFE2gdas_struc(n)%gridDesci(3))

    !- Determine if upscaling or downscaling option:
       select case( LIS_rc%met_interp(findex) )

         case( "average" )   ! Upscaling 

           IF( NINT(LIS_rc%gridDesc(n,1)) .NE. 0 ) THEN ! SY
             write(LIS_logunit,*) "[ERR] RFE2gdas precip implementation supports "
             write(LIS_logunit,*) "  upscaling only to lat/lon run domain grid ... "
             write(LIS_logunit,*) " Program stopping ..."
             call LIS_endrun()
           ENDIF
           IF ( RFE2gdas_struc(n)%gridDesci(10) .GT. LIS_rc%gridDesc(n,10)  ) THEN 
          ! Confirm if LIS run domain resolution is not less than RFE2gdas forcing 
          !  resolution for lat/lon projection 
             write(LIS_logunit,*) "[ERR] LIS lat/lon run domain resolution is less than"
             write(LIS_logunit,*) "  RFE2gdas forcing resolution in upscaling."
             write(LIS_logunit,*) " Program stopping ..."
             call LIS_endrun()
           ENDIF

           write(LIS_logunit,*) &
             "[INFO] RFE2-GDAS reprojection choice selected is average (upscaling)"
#if 0
! - Soni's old way of averaging code:
           allocate(RFE2gdas_struc(n)%stc(LIS_rc%lnc(n), LIS_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%str(LIS_rc%lnc(n), LIS_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%enc(LIS_rc%lnc(n), LIS_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%enr(LIS_rc%lnc(n), LIS_rc%lnr(n)))
           RFE2gdas_struc(n)%stc = 0
           RFE2gdas_struc(n)%str = 0
           RFE2gdas_struc(n)%enc = 0
           RFE2gdas_struc(n)%enr = 0
           call upscaleByAveraging_input(RFE2gdas_struc(n)%gridDesci(:),&
              LIS_rc%gridDesc(n,:), RFE2gdas_struc(n)%mi,&
              LIS_rc%lnc(n), LIS_rc%lnr(n),&
              RFE2gdas_struc(n)%stc, RFE2gdas_struc(n)%str, &
              RFE2gdas_struc(n)%enc, RFE2gdas_struc(n)%enr)
#endif
! - Updated way:
           allocate( RFE2gdas_struc(n)%n111(RFE2gdas_struc(n)%mi) )
           call upscaleByAveraging_input( RFE2gdas_struc(n)%gridDesci(:), &
                   LIS_rc%gridDesc(n,:), RFE2gdas_struc(n)%mi, &
                   LIS_rc%lnc(n)*LIS_rc%lnr(n), RFE2gdas_struc(n)%n111 )


       ! Setting up weights for interpolation:
         case( "bilinear" )
           write(LIS_logunit,*) &
               "[INFO] RFE2-GDAS reprojection choice selected is ",&
                trim(LIS_rc%met_interp(findex))," (interpolation)"

            allocate(RFE2gdas_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2gdas_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2gdas_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2gdas_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2gdas_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2gdas_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2gdas_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(RFE2gdas_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            
            call bilinear_interp_input(n,RFE2gdas_struc(n)%gridDesci(:), &
                 RFE2gdas_struc(n)%n111,RFE2gdas_struc(n)%n121, &
                 RFE2gdas_struc(n)%n211,&
                 RFE2gdas_struc(n)%n221,RFE2gdas_struc(n)%w111, &
                 RFE2gdas_struc(n)%w121,&
                 RFE2gdas_struc(n)%w211,RFE2gdas_struc(n)%w221)

         case( "budget-bilinear" )
           write(LIS_logunit,*) &
               "[INFO] RFE2-GDAS reprojection choice selected is ",&
                trim(LIS_rc%met_interp(findex))," (interpolation)"

            allocate(RFE2gdas_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

            call conserv_interp_input(n,RFE2gdas_struc(n)%gridDesci(:), &
                 RFE2gdas_struc(n)%n112,RFE2gdas_struc(n)%n122,&
                 RFE2gdas_struc(n)%n212,RFE2gdas_struc(n)%n222,&
                 RFE2gdas_struc(n)%w112,RFE2gdas_struc(n)%w122,&
                 RFE2gdas_struc(n)%w212,RFE2gdas_struc(n)%w222)

         case( "neighbor" )
           write(LIS_logunit,*) &
               "[INFO] RFE2-GDAS reprojection choice selected is ",&
                trim(LIS_rc%met_interp(findex))," (interpolation)"

            allocate(RFE2gdas_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            call neighbor_interp_input(n,RFE2gdas_struc(n)%gridDesci(:), &
                 RFE2gdas_struc(n)%n113)

         case default
           write(LIS_logunit,*) "[ERR] This interpolation option not "
           write(LIS_logunit,*) "  supported for RFE2gdas met forcing."
           write(LIS_logunit,*) " Program stopping ... "
           call LIS_endrun()

       end select
       
       if(LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,&
               nint(RFE2gdas_struc(n)%gridDesci(2)),&
               nint(RFE2gdas_struc(n)%gridDesci(3)))
          
       endif

    enddo !  end nest loop
  
  end subroutine init_RFE2gdas

end module RFE2gdas_forcingMod

