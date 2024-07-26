!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module merra2_forcingMod
!BOP
! !MODULE: merra2_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the MERRA2 forcing data.
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt merra2\_struc}
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[nc]
!    Number of columns (along the east west dimension) for the input data
!  \item[nr]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the ECMWF data
!  \item[merra2time1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[merra2time2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[merra2dir]
!    Directory containing the input data
!  \item[merra2hgt_file]
!    File with the terrain height definition for the input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for bilinear interpolation.
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for n. neighbor interpolation.
!  \item[findtime1, findtime2]
!    boolean flags to indicate which time is to be read for
!    temporal interpolation.
!  \end{description}
!
! !USES:
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_merra2      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: merra2_struc

!EOP
  type, public ::  merra2_type_dec
     real         :: ts
     integer      :: nc, nr
     character(len=LDT_CONST_PATH_LEN) :: merra2dir   !MERRA2 Forcing Directory
     real*8       :: merra2time1, merra2time2, ringtime
     logical      :: reset_flag

     integer                :: mi
     real                   :: gridDesc(20)
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
     real, allocatable      :: merraforc1(:,:,:), merraforc2(:,:,:)

     integer                :: uselml
     integer                :: usecorr

     character*140          :: merra2hgt_file

  end type merra2_type_dec

  type(merra2_type_dec), allocatable :: merra2_struc(:)

contains

!BOP
!
! !ROUTINE: init_merra2
! \label{init_merra2}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
!
! !INTERFACE:
  subroutine init_merra2(findex)

! !USES:
    use LDT_coreMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for MERRA2
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_merra2](\ref{readcrd_merra2}) \newline
!     reads the runtime options specified for MERRA2 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    integer :: n
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt

    allocate(merra2_struc(LDT_rc%nnest))
    write(LDT_logunit,fmt=*)"[INFO] Initializing MERRA-2 forcing grid ... "

 !- Read in config file entries:
    call readcrd_merra2()

    do n=1, LDT_rc%nnest
       merra2_struc(n)%ts = 3600  
       call LDT_update_timestep(LDT_rc, n, merra2_struc(n)%ts)
    enddo

    merra2_struc%reset_flag = .false.

    LDT_rc%met_nf(findex) = 14
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true.

    merra2_struc%nc = 576
    merra2_struc%nr = 361
    LDT_rc%met_nc(findex) = merra2_struc(1)%nc
    LDT_rc%met_nr(findex) = merra2_struc(1)%nr

    do n=1,LDT_rc%nnest
       merra2_struc(n)%nc = 576
       merra2_struc(n)%nr = 361

       merra2_struc(n)%gridDesc(:) = 0.

       merra2_struc(n)%gridDesc(1) = 0
       merra2_struc(n)%gridDesc(2) = merra2_struc(n)%nc
       merra2_struc(n)%gridDesc(3) = merra2_struc(n)%nr
       merra2_struc(n)%gridDesc(4) = -90.000
       merra2_struc(n)%gridDesc(5) = -180.000
       merra2_struc(n)%gridDesc(6) = 128
       merra2_struc(n)%gridDesc(7) = 90.000
       merra2_struc(n)%gridDesc(8) = 179.375
       merra2_struc(n)%gridDesc(9) = 0.625
       merra2_struc(n)%gridDesc(10) = 0.5
       merra2_struc(n)%gridDesc(20) = 0

       LDT_rc%met_gridDesc(findex,1:20) = merra2_struc(n)%gridDesc(1:20)

       merra2_struc(n)%mi = merra2_struc(n)%nc*merra2_struc(n)%nr

       ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )

        case( "bilinear" )
          allocate(merra2_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, merra2_struc(n)%gridDesc(:),&
               merra2_struc(n)%n111,merra2_struc(n)%n121,&
               merra2_struc(n)%n211,merra2_struc(n)%n221,&
               merra2_struc(n)%w111,merra2_struc(n)%w121,&
               merra2_struc(n)%w211,merra2_struc(n)%w221)

        case( "budget-bilinear" )
          allocate(merra2_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(merra2_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, merra2_struc(n)%gridDesc(:),&
               merra2_struc(n)%n111,merra2_struc(n)%n121,&
               merra2_struc(n)%n211,merra2_struc(n)%n221,&
               merra2_struc(n)%w111,merra2_struc(n)%w121,&
               merra2_struc(n)%w211,merra2_struc(n)%w221)

          allocate(merra2_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merra2_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merra2_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merra2_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merra2_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merra2_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merra2_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(merra2_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n, merra2_struc(n)%gridDesc(:),&
               merra2_struc(n)%n112,merra2_struc(n)%n122,&
               merra2_struc(n)%n212,merra2_struc(n)%n222,&
               merra2_struc(n)%w112,merra2_struc(n)%w122,&
               merra2_struc(n)%w212,merra2_struc(n)%w222)

        case( "neighbor" )
          allocate(merra2_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n, merra2_struc(n)%gridDesc(:),&
               merra2_struc(n)%n113)

        case( "none" )
          write(LDT_logunit,*) "[INFO] No interpolation applied for MERRA-2 ..."
          write(LDT_logunit,*) " -- Reading native grid to be written out to -- " 

        case default
          write(LDT_logunit,*) '[ERR] Interpolation option '// &
               trim(LDT_rc%met_gridtransform(findex))//&
               ' for MERRA2 forcing is not supported'
          call LDT_endrun()
      end select

      call LDT_registerAlarm("MERRA2 forcing alarm",&
           86400.0,86400.0)
      merra2_struc(n)%startFlag = .true.
      merra2_struc(n)%dayFlag = .true.

      allocate(merra2_struc(n)%merraforc1(&
           LDT_rc%met_nf(findex), 24, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      allocate(merra2_struc(n)%merraforc2(&
           LDT_rc%met_nf(findex), 24, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n)))

      merra2_struc(n)%merraforc1 = LDT_rc%udef
      merra2_struc(n)%merraforc2 = LDT_rc%udef
    enddo

    write(LDT_logunit,*)"[INFO] MERRA-2 time interp option :: ",&
       trim(LDT_rc%met_tinterp(findex))

  end subroutine init_merra2
end module merra2_forcingMod

