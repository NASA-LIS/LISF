!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      module geositbias_forcingMod
!BOP
! !MODULE: geositbias_forcingMod
!
! !REVISION HISTORY:
! 02 Oct 2025: Fadji Maina, initial code (based on geos-it) 
! 07 Jan 2026: Kristen Whitney, initial code for using dynamic lapse rate
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the GEOS-ITbias forcing data.
!  The data is global 0.625-degree lon. by 0.5-degree lat, in latlon
!  projection, and at 1 hourly intervals. The derived data type
!  {\tt geositbias\_struc}
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
!  \item[geositbiastime1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[geositbiastime2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[geositbiasdir]
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
      public :: init_geositbias    ! defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
      public :: geositbias_struc

!EOP
      type, public :: geositbias_type_dec
      real         :: ts
      integer      :: ncold, nrold
      character(len=LIS_CONST_PATH_LEN) :: geositbiasdir   ! GEOS-ITbias Forcing Directory
      real*8       :: geositbiastime1,geositbiastime2
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
      real, allocatable      :: geositbiasforc1(:,:,:),geositbiasforc2(:,:,:)
      real, allocatable      :: lapserate1(:), lapserate2(:)
      integer            :: nvars
      integer            :: uselml
      real*8             :: ringtime
      integer            :: nIter,st_iterid,en_iterid

      real, allocatable :: metdata1(:,:,:)
      real, allocatable :: metdata2(:,:,:)

      integer                 :: use2mwind
      character(len=LIS_CONST_PATH_LEN) :: scaleffile
      integer, allocatable    :: rseed(:,:)
      integer                 :: usedynlapserate
      character(len=LIS_CONST_PATH_LEN) :: dynlapseratedir
      character(len=LIS_CONST_PATH_LEN) :: dynlapseratepfx
      character(len=LIS_CONST_PATH_LEN) :: dynlapseratesfx
      integer                 :: applydynlapseratecutoff
      real                    :: dynlapseratemincutoff
      real                    :: dynlapseratemaxcutoff
      end type geositbias_type_dec

      type(geositbias_type_dec), allocatable :: geositbias_struc(:)

      contains

!BOP
!
! !ROUTINE: init_geositbias
! \label{init_geositbias}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
!
! !INTERFACE:
      subroutine init_geositbias(findex)

! !USES:
      use LIS_coreMod
      use LIS_timeMgrMod
      use LIS_logMod
      use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native

      implicit none
! !AGRUMENTS:
      integer, intent(in) :: findex

! !DESCRIPTION:
!  Defines the native resolution of the input forcing for GEOS-ITbias
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_geositbias](\ref{readcrd_geositbias}) \newline
!     reads the runtime options specified for GEOS-ITbias data
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
      
      allocate(geositbias_struc(LIS_rc%nnest))

      do n = 1,LIS_rc%nnest
         geositbias_struc(n)%ncold = 187
         geositbias_struc(n)%nrold = 131
      enddo

      call readcrd_geositbias()
      LIS_rc%met_nf(findex) = 4

      geositbias_struc%reset_flag = .false.

      do n = 1, LIS_rc%nnest
         geositbias_struc(n)%ts = 3600 !check
         call LIS_update_timestep(LIS_rc,n,geositbias_struc(n)%ts)
      enddo
      
      gridDesci = 0
      do n = 1,LIS_rc%nnest
         gridDesci(n,1)  = 0
         gridDesci(n,2)  = geositbias_struc(n)%ncold
         gridDesci(n,3)  = geositbias_struc(n)%nrold
         gridDesci(n,4)  = 7.0
         gridDesci(n,5)  = -168.75
         gridDesci(n,6)  = 128
         gridDesci(n,7)  = 72.0
         gridDesci(n,8)  = -52.5
         gridDesci(n,9)  = 0.625
         gridDesci(n,10) = 0.5
         gridDesci(n,20) = 0

         geositbias_struc(n)%mi = geositbias_struc(n)%ncold*geositbias_struc(n)%nrold

       ! Setting up weights for Interpolation
         if (trim(LIS_rc%met_interp(findex)).eq."bilinear") then
            allocate(geositbias_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            call bilinear_interp_input(n,gridDesci(n,:),               &
                   geositbias_struc(n)%n111,geositbias_struc(n)%n121,          &
                   geositbias_struc(n)%n211,geositbias_struc(n)%n221,          &
                   geositbias_struc(n)%w111,geositbias_struc(n)%w121,          &
                   geositbias_struc(n)%w211,geositbias_struc(n)%w221)

         elseif (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
            allocate(geositbias_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            allocate(geositbias_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            call bilinear_interp_input(n,gridDesci(n,:),               &
                   geositbias_struc(n)%n111,geositbias_struc(n)%n121,          &
                   geositbias_struc(n)%n211,geositbias_struc(n)%n221,          &
                   geositbias_struc(n)%w111,geositbias_struc(n)%w121,          &
                   geositbias_struc(n)%w211,geositbias_struc(n)%w221)

            allocate(geositbias_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geositbias_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geositbias_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geositbias_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geositbias_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geositbias_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geositbias_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            allocate(geositbias_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
            call conserv_interp_input(n,gridDesci(n,:),                &
                   geositbias_struc(n)%n112,geositbias_struc(n)%n122,          &
                   geositbias_struc(n)%n212,geositbias_struc(n)%n222,          &
                   geositbias_struc(n)%w112,geositbias_struc(n)%w122,          &
                   geositbias_struc(n)%w212,geositbias_struc(n)%w222)

         elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
            allocate(geositbias_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
            call neighbor_interp_input(n,gridDesci(n,:),               &
                                       geositbias_struc(n)%n113)

         else
            write(LIS_logunit,*) '[ERR] Interpolation option '//       &
                                  trim(LIS_rc%met_interp(findex))//    &
                              ' for GEOS-ITbias forcing is not supported'
            call LIS_endrun()
         endif

         if (geositbias_struc(n)%usedynlapserate.eq.1) then
            allocate(geositbias_struc(n)%lapserate1(LIS_rc%ngrid(n)))
            allocate(geositbias_struc(n)%lapserate2(LIS_rc%ngrid(n)))
         endif

         call LIS_registerAlarm("GEOS-ITbias forcing alarm",86400.0,86400.0)
         geositbias_struc(n)%startFlag = .true.
         geositbias_struc(n)%dayFlag = .true.

         geositbias_struc(n)%nvars = 4

         allocate(geositbias_struc(n)%geositbiasforc1(1,                       &
                  geositbias_struc(n)%nvars, LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(geositbias_struc(n)%geositbiasforc2(1,&
                  geositbias_struc(n)%nvars, LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         geositbias_struc(n)%st_iterid = 1
         geositbias_struc(n)%en_iterId = 1
         geositbias_struc(n)%nIter = 1

         allocate(geositbias_struc(n)%metdata1(1,LIS_rc%met_nf(findex),    &
                  LIS_rc%ngrid(n)))
         allocate(geositbias_struc(n)%metdata2(1,LIS_rc%met_nf(findex),    &
                  LIS_rc%ngrid(n)))

         geositbias_struc(n)%metdata1 = 0
         geositbias_struc(n)%metdata2 = 0

         geositbias_struc(n)%geositbiasforc1 = LIS_rc%udef
         geositbias_struc(n)%geositbiasforc2 = LIS_rc%udef

         if (LIS_rc%met_ecor(findex).eq."lapse-rate" .or.             &
             LIS_rc%met_ecor(findex).eq."lapse-rate and slope-aspect" .or. &
             LIS_rc%met_ecor(findex) == "micromet" ) then
            call read_geositbias_elev(n,findex)
         endif

! Set up precipitation climate downscaling:
         if (LIS_rc%pcp_downscale(findex).ne.0) then
            call LIS_init_pcpclimo_native(n,findex,geositbias_struc(n)%ncold,&
                                                   geositbias_struc(n)%nrold)
         endif

      enddo                     ! End nest loop

      end subroutine init_geositbias
      end module geositbias_forcingMod

