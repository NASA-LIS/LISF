!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
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
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \end{description}
!
! !USES:
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
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
     character(len=LDT_CONST_PATH_LEN) :: era5cdsdir   !ERA5 Forcing Directory
     character(len=LDT_CONST_PATH_LEN) :: era5cdshgt_file
     real*8       :: era5cdstime1,era5cdstime2
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

     integer            :: nmif
     integer            :: uselml
     integer            :: tdimsize

     real*8             :: ringtime
     
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
! 17 Apr 2025: Hiroko Beudoing, adopted ERA5 routines for the public CDS
!                               data format
!
! !INTERFACE:
  subroutine init_era5cds(findex)

! !USES:
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_logMod,     only : LDT_logunit, LDT_endrun
    use LDT_timeMgrMod
    use map_utils,      only : proj_latlon

#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for ERA5
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LDT interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_era5cds](\ref{readcrd_era5cds}) \newline
!     reads the runtime options specified for ERA5 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP

    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n
    integer :: ftn
    real :: gridDesc(20)

    allocate(era5cds_struc(LDT_rc%nnest))

    write(LDT_logunit,*) "Reading the ERA5CDS forcing data"

!-  Read in config ERA5CDS inputs:
    call readcrd_era5cds(findex)

    do n=1, LDT_rc%nnest
       era5cds_struc(n)%ts = 3600  !hour
       call LDT_update_timestep(LDT_rc, n, era5cds_struc(n)%ts)
    enddo

    era5cds_struc%reset_flag = .false.
    era5cds_struc(:)%nmif = 9

    ! Metforcing and parameter grid info:
    LDT_rc%met_proj(findex) = "latlon"

    era5cds_struc(:)%ncold = 1440
    era5cds_struc(:)%nrold = 720
    gridDesc = 0
    gridDesc(1) = 0
    gridDesc(2) = era5cds_struc(1)%ncold
    gridDesc(3) = era5cds_struc(1)%nrold
    gridDesc(4) = -89.875
    gridDesc(5) = -179.875
    gridDesc(6) = 128
    gridDesc(7) = 89.875
    gridDesc(8) = 179.875
    gridDesc(9) = 0.25
    gridDesc(10) = 0.25
    gridDesc(20) = 0

    LDT_rc%met_nc(findex) = era5cds_struc(1)%ncold
    LDT_rc%met_nr(findex) = era5cds_struc(1)%nrold
    LDT_rc%met_gridDesc(findex,1:20) = gridDesc(1:20)

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return

    LDT_rc%met_nf(findex) = 9
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true. 

    do n=1,LDT_rc%nnest

       era5cds_struc(n)%findtime1 = 0
       era5cds_struc(n)%findtime2 = 0

       era5cds_struc(n)%mi = era5cds_struc(n)%ncold*era5cds_struc(n)%nrold

       ! ERA5 accumulation data starts at 7z on 1940-01-01
       yr1 = 1940
       mo1 = 01
       da1 = 01
       hr1 = 07
       mn1 = 0; ss1 = 0
       call LDT_date2time( era5cds_struc(n)%validstart,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

       if( gridDesc(9)  == LDT_rc%gridDesc(n,9) .and. &
           gridDesc(10) == LDT_rc%gridDesc(n,10).and. &
           LDT_rc%gridDesc(n,1) == proj_latlon .and. &
           LDT_rc%met_gridtransform(findex) .ne. "neighbor" ) then
         write(LDT_logunit,*) "[ERR]  The ERA5CDS 0.25 deg grid was selected for the"
         write(LDT_logunit,*) "  LDT run domain; however, 'bilinear', 'budget-bilinear',"
         write(LDT_logunit,*) "  or some other unknown option was selected to spatially"
         write(LDT_logunit,*) "  downscale the grid, which will cause errors during runtime."
         write(LDT_logunit,*) "Program stopping ..."
         call LDT_endrun()
       endif

       ! Setting up weights for Interpolation
       select case( LDT_rc%met_gridtransform(findex) )
        case("bilinear") 
          allocate(era5cds_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesc(:),&
               era5cds_struc(n)%n111,era5cds_struc(n)%n121,&
               era5cds_struc(n)%n211,era5cds_struc(n)%n221,&
               era5cds_struc(n)%w111,era5cds_struc(n)%w121,&
               era5cds_struc(n)%w211,era5cds_struc(n)%w221)

        case("budget-bilinear") 
          allocate(era5cds_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5cds_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesc(:),&
               era5cds_struc(n)%n111,era5cds_struc(n)%n121,&
               era5cds_struc(n)%n211,era5cds_struc(n)%n221,&
               era5cds_struc(n)%w111,era5cds_struc(n)%w121,&
               era5cds_struc(n)%w211,era5cds_struc(n)%w221)

          allocate(era5cds_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5cds_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5cds_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5cds_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5cds_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5cds_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5cds_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5cds_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n, gridDesc(:),&
               era5cds_struc(n)%n112,era5cds_struc(n)%n122,&
               era5cds_struc(n)%n212,era5cds_struc(n)%n222,&
               era5cds_struc(n)%w112,era5cds_struc(n)%w122,&
               era5cds_struc(n)%w212,era5cds_struc(n)%w222)

        case("neighbor")
          allocate(era5cds_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n, gridDesc(:),&
               era5cds_struc(n)%n113)

        case default
          write(LDT_logunit,*) '[ERR] Interpolation option '// &
               trim(LDT_rc%met_gridtransform(findex))//&
               ' for ERA5CDS forcing is not supported'
          call LDT_endrun()
       end select

       call LDT_registerAlarm("ERA5CDS forcing alarm",&
            86400.0,86400.0)
       era5cds_struc(n)%startFlag = .true.
       era5cds_struc(n)%dayFlag = .true.

!       if( LDT_rc%met_ecor(findex).eq."lapse-rate")  then
!          call read_era5cdselev_ldtproc(n, findex, 0)
!       endif

       era5cds_struc(n)%mon = -1
       allocate(era5cds_struc(n)%tair(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5cds_struc(n)%qair(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5cds_struc(n)%uwind(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5cds_struc(n)%vwind(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5cds_struc(n)%ps(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5cds_struc(n)%rainf(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5cds_struc(n)%crainf(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5cds_struc(n)%swd(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5cds_struc(n)%lwd(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))

       allocate(era5cds_struc(n)%prev_rainf(LDT_rc%lnc(n)*LDT_rc%lnr(n),7))
       allocate(era5cds_struc(n)%prev_crainf(LDT_rc%lnc(n)*LDT_rc%lnr(n),7))
       allocate(era5cds_struc(n)%prev_swd(LDT_rc%lnc(n)*LDT_rc%lnr(n),7))
       allocate(era5cds_struc(n)%prev_lwd(LDT_rc%lnc(n)*LDT_rc%lnr(n),7))

    enddo   ! End nest loop
    
    
  end subroutine init_era5cds
end module era5cds_forcingMod

