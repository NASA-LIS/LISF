!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! See RELEASE_NOTES.txt for more information.
!
! The LDT source code and documentation are not in the public domain
! and may not be freely distributed.  Only qualified entities may receive 
! the source code and documentation. 
!
! Qualified entities must be covered by a Software Usage Agreement. 
! The Software Usage Agreement contains all the terms and conditions
! regarding the release of the LDT software.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See the Software Usage Agreement for the full disclaimer of warranty.
!
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
!  The implementation in LDT has the derived data type {\tt RFE2gdas\_struc}
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
  public :: init_RFE2gdas      !defines the native resolution of
                               !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: RFE2gdas_struc

!EOP

  type, public :: RFE2gdas_type_dec
     real                     :: ts
     integer                  :: nc, nr
     character(len=LDT_CONST_PATH_LEN) :: RFE2gdasDir
     real*8                   :: RFE2gdasEndTime
     type(ESMF_Time)          :: startTime
     real*8                   :: st_real
     type(ESMF_TimeInterval)  :: timeStep
     integer                  :: mi
     real                     :: gridDesci(20)

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

  end type RFE2gdas_type_dec

  type(RFE2gdas_type_dec), allocatable :: RFE2gdas_struc(:)

contains

!BOP
!
! !ROUTINE: init_RFE2gdas
!  \label{init_RFE2gdas}
!
! !REVISION HISTORY:
! 26 MAY 2010: Soni Yatheendradas; Initial LDT version for FEWSNET
!
! !INTERFACE:
  subroutine init_RFE2gdas(findex)
! !USES:
    use LDT_coreMod,    only : LDT_rc
    use LDT_timeMgrMod, only : LDT_update_timestep
    use LDT_logMod,     only : LDT_logunit, LDT_endrun

    implicit none

    integer,  intent(in)     :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for CMAP
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
    integer :: n
 
    allocate(RFE2gdas_struc(LDT_rc%nnest))

    write(LDT_logunit,fmt=*)"[INFO] Initializing RFE2-GDAS forcing grid ... "

  ! Read entries from config file
    call readRFE2gdascrd()

    LDT_rc%met_nf(findex) = 2
    LDT_rc%met_ts(findex) = 21600

    RFE2gdas_struc%nc = 751
    RFE2gdas_struc%nr = 801
    LDT_rc%met_nc(findex) = RFE2gdas_struc(1)%nc
    LDT_rc%met_nr(findex) = RFE2gdas_struc(1)%nr

  ! Define parameters for geographic projection:
    LDT_rc%met_proj(findex)       = "latlon"
    LDT_rc%met_gridDesc(findex,1) = 0      ! lat/lon projection
    LDT_rc%met_gridDesc(findex,2) = RFE2gdas_struc(1)%nc
    LDT_rc%met_gridDesc(findex,3) = RFE2gdas_struc(1)%nr
    LDT_rc%met_gridDesc(findex,4) = -40.0  ! latitude of the lower left corner CELL CENTER of the domain?
    LDT_rc%met_gridDesc(findex,5) = -20.0  ! longitude of the lower left corner CELL CENTER of the domain?
    LDT_rc%met_gridDesc(findex,6) = 128    ! not used
    LDT_rc%met_gridDesc(findex,7) = 40.0   ! latitude of the upper right corner CELL CENTER of the domain?
    LDT_rc%met_gridDesc(findex,8) = 55.0   ! longitude of the upper right corner CELL CENTER of the domain?
    LDT_rc%met_gridDesc(findex,9) = 0.1    ! spatial resolution (in degrees) along the N-S dimension 
    LDT_rc%met_gridDesc(findex,10) = 0.1   ! spatial resolution (in degrees) along the E-W dimension 
    LDT_rc%met_gridDesc(findex,20) = 64   ! used to specify the ordering of data 

     ! (Non-divisible by 32 indicates E-W ordering, else N-S ordering) 
     !  Set now to 255 per interp/get_fieldpos explanation:  
     !   E-W ordering indicates elements located in row-order array, i.e., one row, then next and so on. 

    write(LDT_logunit,*)'[INFO] Number of RFE2gdas grid columns = ', &
                          RFE2gdas_struc(1)%nc
    write(LDT_logunit,*)'[INFO] Number of RFE2gdas grid rows = ', &
                          RFE2gdas_struc(1)%nr

    do n=1,LDT_rc%nnest

       RFE2gdas_struc(n)%gridDesci(:) = LDT_rc%met_gridDesc(findex,:)

       RFE2gdas_struc(n)%mi = NINT(LDT_rc%met_gridDesc(findex,2))* &
                              NINT(LDT_rc%met_gridDesc(findex,3))

       RFE2gdas_struc(n)%ts = 21600
       call LDT_update_timestep(LDT_rc, n, RFE2gdas_struc(n)%ts)


    !- Determine if upscaling or downscaling option:
       select case( LDT_rc%met_gridtransform(findex) )

         case( "average" )   ! Upscaling

         IF (NINT(LDT_rc%gridDesc(n,1)) .NE. 0) THEN ! SY
           write(LDT_logunit,*) "[ERR] RFE2gdas precip implementation supports "
           write(LDT_logunit,*) "  upscaling only to geographic run domain grid ... "
           write(LDT_logunit,*) " Program stopping ..."
           call LDT_endrun()
         ENDIF

         IF ( LDT_rc%met_gridDesc(findex,10) .GT. LDT_rc%gridDesc(n,10)  ) THEN 
        ! Confirm if LDT run domain resolution is not less than RFE2gdas forcing 
        !  resolution for geographic projection 
           write(LDT_logunit,*) "[ERR] LDT geographic run domain resolution is less than"
           write(LDT_logunit,*) "  RFE2gdas forcing resolution in upscaling."
           write(LDT_logunit,*) " Program stopping ..."
           call LDT_endrun()
         ENDIF

         write(LDT_logunit,*) &
             "[INFO] RFE2-GDAS reprojection choice selected is average (upscaling)"
#if 0
! - Soni's old way of averaging code:
         allocate(RFE2gdas_struc(n)%stc(LDT_rc%lnc(n), LDT_rc%lnr(n)))
         allocate(RFE2gdas_struc(n)%str(LDT_rc%lnc(n), LDT_rc%lnr(n)))
         allocate(RFE2gdas_struc(n)%enc(LDT_rc%lnc(n), LDT_rc%lnr(n)))
         allocate(RFE2gdas_struc(n)%enr(LDT_rc%lnc(n), LDT_rc%lnr(n)))
         RFE2gdas_struc(n)%stc = 0
         RFE2gdas_struc(n)%str = 0
         RFE2gdas_struc(n)%enc = 0
         RFE2gdas_struc(n)%enr = 0
         call upscaleByAveraging_input(LDT_rc%met_gridDesc(findex,:),&
              LDT_rc%gridDesc(n,:), RFE2gdas_struc(n)%mi,&
              LDT_rc%lnc(n), LDT_rc%lnr(n),&
              RFE2gdas_struc(n)%stc, RFE2gdas_struc(n)%str, &
              RFE2gdas_struc(n)%enc, RFE2gdas_struc(n)%enr)
#endif
! - Updated way:
           allocate( RFE2gdas_struc(n)%n111(RFE2gdas_struc(n)%mi) )
           call upscaleByAveraging_input( RFE2gdas_struc(n)%gridDesci(:), &
                   LDT_rc%gridDesc(n,:), RFE2gdas_struc(n)%mi, &
                   LDT_rc%lnc(n)*LDT_rc%lnr(n), RFE2gdas_struc(n)%n111 )

     ! Setting up weights for interpolation:
       case( "bilinear" )
           write(LDT_logunit,*) &
               "[INFO] RFE2-GDAS reprojection choice selected is ",&
                trim(LDT_rc%met_gridtransform(findex))," (interpolation)"

           allocate(RFE2gdas_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
           allocate(RFE2gdas_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

           call bilinear_interp_input(n, LDT_rc%met_gridDesc(findex,:), &
                 RFE2gdas_struc(n)%n111,RFE2gdas_struc(n)%n121, &
                 RFE2gdas_struc(n)%n211,&
                 RFE2gdas_struc(n)%n221,RFE2gdas_struc(n)%w111, &
                 RFE2gdas_struc(n)%w121,&
                 RFE2gdas_struc(n)%w211,RFE2gdas_struc(n)%w221)

         case( "budget-bilinear" )
           write(LDT_logunit,*) &
               "[INFO] RFE2-GDAS reprojection choice selected is ",&
                trim(LDT_rc%met_gridtransform(findex))," (interpolation)"

            allocate(RFE2gdas_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            allocate(RFE2gdas_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
            call conserv_interp_input(n, LDT_rc%met_gridDesc(findex,:), &
                 RFE2gdas_struc(n)%n112,RFE2gdas_struc(n)%n122,&
                 RFE2gdas_struc(n)%n212,RFE2gdas_struc(n)%n222,&
                 RFE2gdas_struc(n)%w112,RFE2gdas_struc(n)%w122,&
                 RFE2gdas_struc(n)%w212,RFE2gdas_struc(n)%w222)

         case( "neighbor" )
           write(LDT_logunit,*) &
               "[INFO] RFE2-GDAS reprojection choice selected is ",&
                trim(LDT_rc%met_gridtransform(findex))," (interpolation)"

            allocate(RFE2gdas_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
            call neighbor_interp_input(n, LDT_rc%met_gridDesc(findex,:), &
                 RFE2gdas_struc(n)%n113)

         case default
           write(LDT_logunit,*) "[ERR] This interpolation option not "
           write(LDT_logunit,*) " supported for RFE2gdas supp forcing."
           write(LDT_logunit,*) "Program stopping ... "
           call LDT_endrun()

        end select 

     enddo ! end nest loop
  
  end subroutine init_RFE2gdas

end module RFE2gdas_forcingMod
