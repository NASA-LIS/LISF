!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module merra2bias_forcingMod
!BOP
! !MODULE: merra2bias_forcingMod
!
! !REVISION HISTORY:
! 29 Jun 2026: Kristen Whitney, initial code (based on geos-itbias)
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the MERRA2bias forcing data.
!  MERRA2bias consists of bias-adjusted MERRA2 surface forcing
!  fields that are generated upstream of LIS. LIS does not perform
!  the MERRA2 to MERRA2bias adjustment; it reads the prepared
!  MERRA2bias files and applies the requested LIS spatial
!  interpolation and downscaling procedures.
!
!  The MERRA2bias files used by this reader contain four hourly
!  surface forcing fields: 2-m air temperature, 2-m specific humidity,
!  surface pressure, and downward longwave radiation. The fields are
!  bias-adjusted using climatological reference datasets generated
!  outside of LIS.
!
!  The data are on a 0.625-degree longitude by 0.5-degree latitude
!  grid, in latlon projection, and at 1 hourly intervals. The derived
!  data type {\tt merra2bias\_struc} includes the variables that
!  specify the runtime options, and the weights and neighbor
!  information to be used for spatial interpolation. They are
!  described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the data
!  \item[merra2biastime1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[merra2biastime2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[merra2biasdir]
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
  public :: init_merra2bias    ! defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: merra2bias_struc

!EOP
  type, public :: merra2bias_type_dec
     real         :: ts
     integer      :: ncold, nrold
     character(len=LIS_CONST_PATH_LEN) :: merra2biasdir   ! MERRA2bias Forcing Directory
     real*8       :: merra2biastime1,merra2biastime2
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
     real, allocatable      :: merra2biasforc1(:,:,:),merra2biasforc2(:,:,:)
     real, allocatable      :: lapserate1(:), lapserate2(:)
     integer            :: nvars
     integer            :: uselml
     real*8             :: ringtime
     integer            :: nIter,st_iterid,en_iterid

     real, allocatable :: metdata(:,:,:)

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
  end type merra2bias_type_dec

  type(merra2bias_type_dec), allocatable :: merra2bias_struc(:)

contains

!BOP
!
! !ROUTINE: init_merra2bias
! \label{init_merra2bias}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
! 17 Jun 2026: Kristen Whitney, clarified MERRA2bias dataset
!              description and simplified metadata storage.
!
! !INTERFACE:
  subroutine init_merra2bias(findex)

! !USES:
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native

    implicit none

! !AGRUMENTS:
    integer, intent(in) :: findex

! !DESCRIPTION:
!  Defines the native resolution of the input forcing for MERRA2bias
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_merra2bias](\ref{readcrd_merra2bias}) \newline
!     reads the runtime options specified for MERRA2bias data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real :: gridDesci(LIS_rc%nnest,50)
    integer :: n

    external :: readcrd_merra2bias
    external :: bilinear_interp_input
    external :: conserv_interp_input
    external :: neighbor_interp_input
    external :: read_merra2bias_elev

    allocate(merra2bias_struc(LIS_rc%nnest))

    do n = 1,LIS_rc%nnest
       merra2bias_struc(n)%ncold = 187
       merra2bias_struc(n)%nrold = 131
    enddo

    call readcrd_merra2bias()
    LIS_rc%met_nf(findex) = 4

    merra2bias_struc%reset_flag = .false.

    do n = 1, LIS_rc%nnest
       merra2bias_struc(n)%ts = 3600 !check
       call LIS_update_timestep(LIS_rc,n,merra2bias_struc(n)%ts)
    enddo

    gridDesci = 0
    do n = 1,LIS_rc%nnest
       gridDesci(n,1)  = 0
       gridDesci(n,2)  = merra2bias_struc(n)%ncold
       gridDesci(n,3)  = merra2bias_struc(n)%nrold
       gridDesci(n,4)  = 7.0
       gridDesci(n,5)  = -168.75
       gridDesci(n,6)  = 128
       gridDesci(n,7)  = 72.0
       gridDesci(n,8)  = -52.5
       gridDesci(n,9)  = 0.625
       gridDesci(n,10) = 0.5
       gridDesci(n,20) = 0

       merra2bias_struc(n)%mi = merra2bias_struc(n)%ncold * &
            merra2bias_struc(n)%nrold

       ! Setting up weights for Interpolation
       if (trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(merra2bias_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n,gridDesci(n,:),               &
               merra2bias_struc(n)%n111,merra2bias_struc(n)%n121,    &
               merra2bias_struc(n)%n211,merra2bias_struc(n)%n221,    &
               merra2bias_struc(n)%w111,merra2bias_struc(n)%w121,    &
               merra2bias_struc(n)%w211,merra2bias_struc(n)%w221)

       elseif (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(merra2bias_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(merra2bias_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n,gridDesci(n,:),               &
               merra2bias_struc(n)%n111,merra2bias_struc(n)%n121,    &
               merra2bias_struc(n)%n211,merra2bias_struc(n)%n221,    &
               merra2bias_struc(n)%w111,merra2bias_struc(n)%w121,    &
               merra2bias_struc(n)%w211,merra2bias_struc(n)%w221)

          allocate(merra2bias_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merra2bias_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merra2bias_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merra2bias_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merra2bias_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merra2bias_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merra2bias_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(merra2bias_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          call conserv_interp_input(n,gridDesci(n,:),                &
               merra2bias_struc(n)%n112,merra2bias_struc(n)%n122,    &
               merra2bias_struc(n)%n212,merra2bias_struc(n)%n222,    &
               merra2bias_struc(n)%w112,merra2bias_struc(n)%w122,    &
               merra2bias_struc(n)%w212,merra2bias_struc(n)%w222)

       elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(merra2bias_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n,gridDesci(n,:),               &
               merra2bias_struc(n)%n113)

       else
          write(LIS_logunit,*) '[ERR] Interpolation option '//       &
               trim(LIS_rc%met_interp(findex))//                     &
               ' for MERRA2bias forcing is not supported'
          call LIS_endrun()
       endif

       if (merra2bias_struc(n)%usedynlapserate.eq.1) then
          allocate(merra2bias_struc(n)%lapserate1(LIS_rc%ngrid(n)))
          allocate(merra2bias_struc(n)%lapserate2(LIS_rc%ngrid(n)))
       endif

       call LIS_registerAlarm("MERRA2bias forcing alarm",86400.0,86400.0)
       merra2bias_struc(n)%startFlag = .true.
       merra2bias_struc(n)%dayFlag = .true.

       merra2bias_struc(n)%nvars = 4

       allocate(merra2bias_struc(n)%merra2biasforc1(1,               &
            merra2bias_struc(n)%nvars, LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       allocate(merra2bias_struc(n)%merra2biasforc2(1,               &
            merra2bias_struc(n)%nvars, LIS_rc%lnc(n)*LIS_rc%lnr(n)))
       merra2bias_struc(n)%st_iterid = 1
       merra2bias_struc(n)%en_iterId = 1
       merra2bias_struc(n)%nIter = 1

       allocate(merra2bias_struc(n)%metdata(1,LIS_rc%met_nf(findex),    &
            LIS_rc%ngrid(n)))

       merra2bias_struc(n)%metdata = 0

       merra2bias_struc(n)%merra2biasforc1 = LIS_rc%udef
       merra2bias_struc(n)%merra2biasforc2 = LIS_rc%udef

       if (LIS_rc%met_ecor(findex).eq."lapse-rate" .or.             &
            LIS_rc%met_ecor(findex).eq."lapse-rate and slope-aspect" .or. &
            LIS_rc%met_ecor(findex) == "micromet" ) then
          call read_merra2bias_elev(n,findex)
       endif

! Set up precipitation climate downscaling:
       if (LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,merra2bias_struc(n)%ncold,&
               merra2bias_struc(n)%nrold)
       endif

    enddo                     ! End nest loop

  end subroutine init_merra2bias
end module merra2bias_forcingMod

