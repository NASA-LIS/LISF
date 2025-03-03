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
!  The implementation in LDT has the derived data type {\tt RFE2Daily\_struc}
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
!  \item[gridDesc]
!    Input grid description parameters
!  \item[mi]  
!    Number of points in the input grid
!  \item[n111,n121,n211,n221] 
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for n. neighbor interpolation.
!  \end{description}
!
! !USES:
  use ESMF
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

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
     real                     :: ts
     integer                  :: nc, nr
     character(len=LDT_CONST_PATH_LEN)             :: RFE2DailyDir
     real*8                   :: RFE2DailyEndTime
     type(ESMF_Time)          :: startTime
     real*8                   :: st_real
     type(ESMF_TimeInterval)  :: timeStep
     integer                  :: mi
     integer                  :: hour_offset

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

  end type RFE2Daily_type_dec

  type(RFE2Daily_type_dec), allocatable :: RFE2Daily_struc(:)

contains

!BOP
!
! !ROUTINE: init_RFE2Daily
!  \label{init_RFE2Daily}
!
! !REVISION HISTORY:
!  26 MAY 2010: Soni Yatheendradas; Initial LDT version for FEWSNET
!  20 Mar 2013; KR Arsenault, Cleaned up code and documentation
!
! !INTERFACE:
  subroutine init_RFE2Daily(findex)

! !USES:
   use LDT_coreMod,    only : LDT_rc
   use LDT_logMod,     only : LDT_verify, LDT_logunit, LDT_endrun
   use LDT_timeMgrMod, only : LDT_update_timestep

   implicit none
   integer,  intent(in)  :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for RFE2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
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
   
   integer  :: n
   integer  :: status
   real     :: gridDesci(20)
 
   allocate(RFE2Daily_struc(LDT_rc%nnest))

   write(LDT_logunit,fmt=*) "MSG: Initializing RFE2 Daily forcing grid ... "

 !- Read in daily RFE2 ldt.config inputs:
    call readcrd_RFE2Daily()

    LDT_rc%met_nf(findex) = 2  ! Number of RFE2 met variables 
    LDT_rc%met_ts(findex) = 24*60*60
    LDT_rc%met_validhr(findex) = 6       ! Daily RFE2 valid at 6Z
    LDT_rc%met_proj(findex)       = "latlon"

    RFE2Daily_struc%nc = 751
    RFE2Daily_struc%nr = 801
    LDT_rc%met_nc(findex) = RFE2Daily_struc(1)%nc
    LDT_rc%met_nr(findex) = RFE2Daily_struc(1)%nr

  ! Define parameters for geographic projection:
    LDT_rc%met_gridDesc(findex,1) = 0      ! lat/lon projection
    LDT_rc%met_gridDesc(findex,2) = RFE2Daily_struc(1)%nc
    LDT_rc%met_gridDesc(findex,3) = RFE2Daily_struc(1)%nr
    LDT_rc%met_gridDesc(findex,4) = -40.05 ! lat of lower left grid-cell
                                           ! (upper left corner)
    LDT_rc%met_gridDesc(findex,5) = -19.95 ! long of lower left grid-cell
                                             ! (upper left corner)
    LDT_rc%met_gridDesc(findex,6) = 128
    LDT_rc%met_gridDesc(findex,7) = 39.95  ! lat of upper right grid-cell
                                             ! (upper left corner)
    LDT_rc%met_gridDesc(findex,8) = 55.05  ! long of upper right grid-cell
                                             ! (upper left corner)
    LDT_rc%met_gridDesc(findex,9) = 0.1    ! spatial res (in degrees)
                                             ! along E-W
    LDT_rc%met_gridDesc(findex,10)= 0.1    ! spatial res (in degrees)
                                             ! along N-S
    LDT_rc%met_gridDesc(findex,20)= 64

    gridDesci(:) = LDT_rc%met_gridDesc(findex,:)

  ! Used to specify the ordering of data 
  !  (non-divisible by 32 indicates E-W ordering, else N-S ordering) 
  ! Set now to 64 and not 255 of JimG's initial interp/get_fieldpos explanation: 
  !  E-W ordering have elements are located in array row-wise, 
  !  i.e., one row, then next and so on. 

    write(LDT_logunit,*) 'Number of RFE2Daily grid columns = ', &
                         RFE2Daily_struc(1)%nc
    write(LDT_logunit,*) 'Number of RFE2Daily grid rows = ', &
                         RFE2Daily_struc(1)%nr

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

 !- Timestep set and update:
    do n=1, LDT_rc%nnest
       RFE2Daily_struc(n)%ts = 24*60*60
       call LDT_update_timestep(LDT_rc, n, RFE2Daily_struc(n)%ts)
    enddo

 !- Loop over nested domains:
    do n=1,LDT_rc%nnest

       RFE2Daily_struc(n)%mi = RFE2Daily_struc(1)%nc * &
                               RFE2Daily_struc(1)%nr 

       write(LDT_logunit,*) "RFE2Daily reprojection choice selected is ",&
                 trim(LDT_rc%met_gridtransform(findex))," (interpolation)"

    !- Determine if upscaling or downscaling option:
       select case( LDT_rc%met_gridtransform(findex) )

         case( "average" )   ! Upscaling 

           IF( NINT(LDT_rc%gridDesc(n,1)) .NE. 0) THEN ! SY
             write(LDT_logunit,*) "RFE2Daily precip implementation supports"
             write(LDT_logunit,*) "upscaling only to lat/lon run domain grids ... "
             write(LDT_logunit,*) "Program stopping ..."
             call LDT_endrun()
           ENDIF

         ! Confirm if LDT run domain resolution really not less than RFE2Daily forcing
         ! resolution for lat/lon projection 
           IF( gridDesci(10) .GT. LDT_rc%gridDesc(n,10)  ) THEN
             write(LDT_logunit,*) "LDT lat/lon run domain resolution less than "
             write(LDT_logunit,*) "RFE2Daily forcing resolution in upscaling option."
             write(LDT_logunit,*) "Program stopping ..."
             call LDT_endrun()
           ENDIF

           allocate(RFE2Daily_struc(n)%n111(RFE2Daily_struc(n)%mi))
           call upscaleByAveraging_input( gridDesci(:), &
                  LDT_rc%gridDesc(n,:), &
                  RFE2Daily_struc(n)%mi,&
                  LDT_rc%lnc(n)*LDT_rc%lnr(n), &
                  RFE2Daily_struc(n)%n111 )

         case( "bilinear" )

            allocate(RFE2Daily_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
            allocate(RFE2Daily_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

            call bilinear_interp_input(n, gridDesci(:), &
                 RFE2Daily_struc(n)%n111,RFE2Daily_struc(n)%n121, &
                 RFE2Daily_struc(n)%n211,&
                 RFE2Daily_struc(n)%n221,RFE2Daily_struc(n)%w111, &
                 RFE2Daily_struc(n)%w121,&
                 RFE2Daily_struc(n)%w211,RFE2Daily_struc(n)%w221)

         case( "budget-bilinear" )

            allocate(RFE2Daily_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2Daily_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))

            call conserv_interp_input(n, gridDesci(:), &
                 RFE2Daily_struc(n)%n112,RFE2Daily_struc(n)%n122,&
                 RFE2Daily_struc(n)%n212,RFE2Daily_struc(n)%n222,&
                 RFE2Daily_struc(n)%w112,RFE2Daily_struc(n)%w122,&
                 RFE2Daily_struc(n)%w212,RFE2Daily_struc(n)%w222)

         case( "neighbor" )

            allocate(RFE2Daily_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

            call neighbor_interp_input(n, gridDesci(:), &
                 RFE2Daily_struc(n)%n113)

         case default
           write(LDT_logunit,*) "This interp. option not supported for RFE2Daily Precip."
           write(LDT_logunit,*) "Program stopping ... "
           call LDT_endrun()

       end select ! end upscale/downscale option

     enddo !  end nest loop
  
  end subroutine init_RFE2Daily

end module RFE2Daily_forcingMod

