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
module era5_forcingMod
!BOP
! !MODULE: era5_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the ERA5 forcing data.
!  The data is global 1 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt era5\_struc}
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
!  \item[era5time1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[era5time2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[era5dir]
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
  public :: init_era5      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: era5_struc

!EOP
  type, public ::  era5_type_dec

     integer      :: npts
     real         :: ts
     integer      :: ncold, nrold
     character(len=LDT_CONST_PATH_LEN) :: era5dir   !ERA5 Forcing Directory
     character(len=LDT_CONST_PATH_LEN) :: mapfile
     character(len=LDT_CONST_PATH_LEN) :: era5hgt_file
     real*8       :: era5time1,era5time2
     logical      :: reset_flag
     integer      :: mo1,mo2

     real :: gridDesc(50)
     integer, allocatable :: G2P(:,:)
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

     integer            :: nvars
     integer            :: uselml

     real*8             :: ringtime
     
     integer            :: nIter, st_iterid,en_iterid

     real, allocatable      :: tair1(:,:)
     real, allocatable      :: qair1(:,:)
     real, allocatable      :: wind1(:,:)
     real, allocatable      :: ps1(:,:)
     real, allocatable      :: rainf1(:,:)
     real, allocatable      :: snowf1(:,:)
     real, allocatable      :: dirswd1(:,:)
     real, allocatable      :: difswd1(:,:)
     real, allocatable      :: swd1(:,:)
     real, allocatable      :: lwd1(:,:)

     real, allocatable      :: tair2(:,:)
     real, allocatable      :: qair2(:,:)
     real, allocatable      :: wind2(:,:)
     real, allocatable      :: ps2(:,:)
     real, allocatable      :: rainf2(:,:)
     real, allocatable      :: snowf2(:,:)
     real, allocatable      :: dirswd2(:,:)
     real, allocatable      :: difswd2(:,:)
     real, allocatable      :: swd2(:,:)
     real, allocatable      :: lwd2(:,:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type era5_type_dec

  type(era5_type_dec), allocatable :: era5_struc(:)

contains

!BOP
!
! !ROUTINE: init_era5
! \label{init_era5}
!
! !REVISION HISTORY:
! 23 Dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:
  subroutine init_era5(findex)

! !USES:
    use LDT_coreMod
    use LDT_timeMgrMod
    use LDT_logMod

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
!   \item[readcrd\_era5](\ref{readcrd_era5}) \newline
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
    integer :: G2Pid


    allocate(era5_struc(LDT_rc%nnest))

    do n=1,LDT_rc%nnest
       era5_struc(n)%ncold = 1440
       era5_struc(n)%nrold = 720
       era5_struc(n)%npts = 340819
       era5_struc(n)%mo1 = -1
       era5_struc(n)%mo2 = -1

       LDT_rc%met_nc(findex) = era5_struc(n)%ncold
       LDT_rc%met_nr(findex) = era5_struc(n)%nrold
    
       allocate(era5_struc(n)%tair1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%qair1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%wind1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%ps1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%rainf1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%swd1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%lwd1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))

       allocate(era5_struc(n)%tair2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%qair2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%wind2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%ps2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%rainf2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%swd2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(era5_struc(n)%lwd2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))

    enddo

    call readcrd_era5()
    LDT_rc%met_nf(findex) = 8
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true. 


    era5_struc%reset_flag = .false.

    do n=1, LDT_rc%nnest
       era5_struc(n)%ts = 3600  !check
       call LDT_update_timestep(LDT_rc, n, era5_struc(n)%ts)
    enddo

    do n=1,LDT_rc%nnest
       era5_struc(n)%gridDesc = 0
       era5_struc(n)%gridDesc(1) = 0
       era5_struc(n)%gridDesc(2) = era5_struc(n)%ncold
       era5_struc(n)%gridDesc(3) = era5_struc(n)%nrold
       era5_struc(n)%gridDesc(4) = -89.875
       era5_struc(n)%gridDesc(5) = -179.875
       era5_struc(n)%gridDesc(6) = 128
       era5_struc(n)%gridDesc(7) = 89.875
       era5_struc(n)%gridDesc(8) = 179.875
       era5_struc(n)%gridDesc(9) = 0.25
       era5_struc(n)%gridDesc(10) = 0.25
       era5_struc(n)%gridDesc(20) = 0

       LDT_rc%met_gridDesc(findex,1:20) = era5_struc(n)%gridDesc(1:20)

       era5_struc(n)%mi = era5_struc(n)%ncold*era5_struc(n)%nrold

       ! Setting up weights for Interpolation
       if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then
          allocate(era5_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, era5_struc(n)%gridDesc(:),&
               era5_struc(n)%n111,era5_struc(n)%n121,&
               era5_struc(n)%n211,era5_struc(n)%n221,&
               era5_struc(n)%w111,era5_struc(n)%w121,&
               era5_struc(n)%w211,era5_struc(n)%w221)

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then
          allocate(era5_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(era5_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, era5_struc(n)%gridDesc(:),&
               era5_struc(n)%n111,era5_struc(n)%n121,&
               era5_struc(n)%n211,era5_struc(n)%n221,&
               era5_struc(n)%w111,era5_struc(n)%w121,&
               era5_struc(n)%w211,era5_struc(n)%w221)

          allocate(era5_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(era5_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n, era5_struc(n)%gridDesc(:),&
               era5_struc(n)%n112,era5_struc(n)%n122,&
               era5_struc(n)%n212,era5_struc(n)%n222,&
               era5_struc(n)%w112,era5_struc(n)%w122,&
               era5_struc(n)%w212,era5_struc(n)%w222)

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then
          allocate(era5_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n, era5_struc(n)%gridDesc(:),&
               era5_struc(n)%n113)

       else
          write(LDT_logunit,*) '[ERR] Interpolation option '// &
               trim(LDT_rc%met_gridtransform(findex))//&
               ' for ERA5 forcing is not supported'
          call LDT_endrun()
       endif

       call LDT_registerAlarm("ERA5 forcing alarm",&
            86400.0,86400.0)
       era5_struc(n)%startFlag = .true.
       era5_struc(n)%dayFlag = .true.

       era5_struc(n)%nvars = 8

       allocate(era5_struc(n)%metdata1(LDT_rc%met_nf(findex),&
            LDT_rc%ngrid(n)))
       allocate(era5_struc(n)%metdata2(LDT_rc%met_nf(findex),&
            LDT_rc%ngrid(n)))

       era5_struc(n)%metdata1 = 0
       era5_struc(n)%metdata2 = 0

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 

       allocate(era5_struc(n)%G2P(era5_struc(n)%ncold,&
               era5_struc(n)%nrold))
       
       call LDT_verify(nf90_open(path=trim(era5_struc(n)%mapfile), &
            mode=NF90_NOWRITE, &
          ncid=ftn), 'nf90_open failed for '//trim(era5_struc(n)%mapfile))

       call LDT_verify(nf90_inq_varid(ftn,'G2P',G2PId), &
            'nf90_inq_varid failed for G2P in read_era5')

       call LDT_verify(nf90_get_var(ftn,G2PId, era5_struc(n)%G2P),&
            'nf90_get_var failed for G2P in read_era5') 
       call LDT_verify(nf90_close(ftn))
#endif
    enddo   ! End nest loop
    
    
  end subroutine init_era5
end module era5_forcingMod

