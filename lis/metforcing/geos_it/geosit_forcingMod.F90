!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      module geosit_forcingMod
!BOP
! !MODULE: geosit_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the GEOS-IT forcing data.
!  The data is global 0.625-degree lon. by 0.5-degree lat, in latlon
!  projection, and at 1 hourly intervals. The derived data type
!  {\tt geosit\_struc}
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the data
!  \item[geosittime1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[geosittime2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[geositdir]
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
      public :: init_geosit    ! defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
      public :: geosit_struc

!EOP
      type, public :: geosit_type_dec
      real         :: ts
      integer      :: ncold, nrold
      character(len=LIS_CONST_PATH_LEN) :: geositdir   ! GEOS-IT Forcing Directory
      real*8       :: geosittime1,geosittime2
      logical      :: reset_flag

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
      integer                :: findtime1,findtime2
      logical                :: startFlag,dayFlag
      real, allocatable      :: geositforc1(:,:,:),geositforc2(:,:,:)

      integer            :: nvars
      integer            :: uselml
      real*8             :: ringtime
      integer            :: nIter,st_iterid,en_iterid

      real, allocatable :: metdata1(:,:,:)
      real, allocatable :: metdata2(:,:,:)

      integer                 :: use2mwind
      character(len=LIS_CONST_PATH_LEN) :: scaleffile
      integer, allocatable    :: rseed(:,:)
      end type geosit_type_dec

      type(geosit_type_dec), allocatable :: geosit_struc(:)

      contains

!BOP
!
! !ROUTINE: init_geosit
! \label{init_geosit}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
!
! !INTERFACE:
      subroutine init_geosit(findex)

! !USES:
      use LIS_coreMod
      use LIS_timeMgrMod
      use LIS_logMod
      use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native

      implicit none
! !AGRUMENTS:
      integer, intent(in) :: findex

! !DESCRIPTION:
!  Defines the native resolution of the input forcing for GEOS-IT
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_geosit](\ref{readcrd_geosit}) \newline
!     reads the runtime options specified for GEOS-IT data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
      real :: gridDesci(LIS_rc%nnest,50)
      integer :: updoy,yr1,mo1,da1,hr1,mn1,ss1
      real :: upgmt
      integer :: n

      allocate(geosit_struc(LIS_rc%nnest))

      do n = 1,LIS_rc%nnest
         geosit_struc(n)%ncold = 576
         geosit_struc(n)%nrold = 361
      enddo

      call readcrd_geosit()
      LIS_rc%met_nf(findex) = 14

      geosit_struc%reset_flag = .false.

      do n = 1, LIS_rc%nnest
         geosit_struc(n)%ts = 3600 !check
         call LIS_update_timestep(LIS_rc,n,geosit_struc(n)%ts)
      enddo

      gridDesci = 0

      do n = 1,LIS_rc%nnest
         gridDesci(n,1)  = 0
         gridDesci(n,2)  = geosit_struc(n)%ncold
         gridDesci(n,3)  = geosit_struc(n)%nrold
         gridDesci(n,4)  = -90.000
         gridDesci(n,5)  = -180.000
         gridDesci(n,6)  = 128
         gridDesci(n,7)  = 90.000
         gridDesci(n,8)  = 179.375
         gridDesci(n,9)  = 0.625
         gridDesci(n,10) = 0.5
         gridDesci(n,20) = 0

         geosit_struc(n)%mi = geosit_struc(n)%ncold*geosit_struc(n)%nrold

       ! Setting up weights for Interpolation
         if (trim(LIS_rc%met_interp(findex)).eq."bilinear") then
            allocate(geosit_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            call bilinear_interp_input(n,gridDesci(n,:),               &
                   geosit_struc(n)%n111,geosit_struc(n)%n121,          &
                   geosit_struc(n)%n211,geosit_struc(n)%n221,          &
                   geosit_struc(n)%w111,geosit_struc(n)%w121,          &
                   geosit_struc(n)%w211,geosit_struc(n)%w221)

         elseif (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
            allocate(geosit_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geosit_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            call bilinear_interp_input(n,gridDesci(n,:),               &
                   geosit_struc(n)%n111,geosit_struc(n)%n121,          &
                   geosit_struc(n)%n211,geosit_struc(n)%n221,          &
                   geosit_struc(n)%w111,geosit_struc(n)%w121,          &
                   geosit_struc(n)%w211,geosit_struc(n)%w221)

            allocate(geosit_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geosit_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geosit_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geosit_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geosit_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geosit_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geosit_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geosit_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            call conserv_interp_input(n,gridDesci(n,:),                &
                   geosit_struc(n)%n112,geosit_struc(n)%n122,          &
                   geosit_struc(n)%n212,geosit_struc(n)%n222,          &
                   geosit_struc(n)%w112,geosit_struc(n)%w122,          &
                   geosit_struc(n)%w212,geosit_struc(n)%w222)

         elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
            allocate(geosit_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            call neighbor_interp_input(n,gridDesci(n,:),               &
                                       geosit_struc(n)%n113)

         else
            write(LIS_logunit,*) '[ERR] Interpolation option '//       &
                                  trim(LIS_rc%met_interp(findex))//    &
                              ' for GEOS-IT forcing is not supported'
            call LIS_endrun()
         endif

         call LIS_registerAlarm("GEOS-IT forcing alarm",86400.0,86400.0)
         geosit_struc(n)%startFlag = .true.
         geosit_struc(n)%dayFlag = .true.

         geosit_struc(n)%nvars = 14

         allocate(geosit_struc(n)%geositforc1(1,                       &
                  geosit_struc(n)%nvars, LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(geosit_struc(n)%geositforc2(1,&
                  geosit_struc(n)%nvars, LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         geosit_struc(n)%st_iterid = 1
         geosit_struc(n)%en_iterId = 1
         geosit_struc(n)%nIter = 1

         allocate(geosit_struc(n)%metdata1(1,LIS_rc%met_nf(findex),    &
                  LIS_rc%ngrid(n)))
         allocate(geosit_struc(n)%metdata2(1,LIS_rc%met_nf(findex),    &
                  LIS_rc%ngrid(n)))

         geosit_struc(n)%metdata1 = 0
         geosit_struc(n)%metdata2 = 0

         geosit_struc(n)%geositforc1 = LIS_rc%udef
         geosit_struc(n)%geositforc2 = LIS_rc%udef

         if ((LIS_rc%met_ecor(findex).eq."lapse-rate").or.             &
             (LIS_rc%met_ecor(findex).eq."lapse-rate and slope-aspect")) then
            call read_geosit_elev(n,findex)
         endif

! Set up precipitation climate downscaling:
         if (LIS_rc%pcp_downscale(findex).ne.0) then
            call LIS_init_pcpclimo_native(n,findex,geosit_struc(n)%ncold,&
                                                   geosit_struc(n)%nrold)
         endif
      enddo                     ! End nest loop

      end subroutine init_geosit
      end module geosit_forcingMod

