!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module galwemge_forcingMod
!BOP
! !MODULE: galwemge_forcingMod

! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of GALWEM-GE 16-day forecast data used as forcing
!  within LIS.

! REVISION HISTORY:
! 09 May 2022; Yeosang Yoon; Initial Specification
! 05 Apr 2023; Yeosang Yoon; Update code to fit new format

! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_galwemge      !defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: galwemge_struc

  type, public ::  galwemge_type_dec
     real                              :: ts
     integer                           :: nc, nr, vector_len   
     real*8                            :: fcsttime1,fcsttime2
     character(len=LIS_CONST_PATH_LEN) :: odir     !GALWEM-GE forecast forcing Directory
     character*20                      :: runmode
     integer                           :: max_ens_members

     integer, allocatable   :: gindex(:,:)

     integer                :: mi
     !integer                :: day_check1
     !integer                :: day_check2
     
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
     integer                :: fcst_hour
     integer                :: init_yr, init_mo, init_da, init_hr
     real, allocatable      :: metdata1(:,:,:) 
     real, allocatable      :: metdata2(:,:,:)

     integer                :: nmodels   

  end type galwemge_type_dec

  type(galwemge_type_dec), allocatable :: galwemge_struc(:)
!EOP
contains

!BOP
!
! !ROUTINE: init_galwemge
! \label{init_galwemge}
! 
! !INTERFACE:
  subroutine init_galwemge(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !USES: 
    integer, intent(in)  :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GALWEM
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!EOP

    integer :: n
    real    :: gridDesci(LIS_rc%nnest,50)

    write(LIS_logunit,*) "[INFO] Initializing the GALWEM-GE forecast inputs "

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GALWEM-GE forecast forcing reader'
       write(LIS_logunit,*) '[ERR] is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR] LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(galwemge_struc(LIS_rc%nnest))
    call readcrd_galwemge()

    do n=1, LIS_rc%nnest
       galwemge_struc(n)%ts = 3600*3  !check
       call LIS_update_timestep(LIS_rc, n, galwemge_struc(n)%ts)
    enddo

    do n=1, LIS_rc%nnest
       galwemge_struc(:)%nc = 720  ! galwem-ge
       galwemge_struc(:)%nr = 361
    enddo

    ! 8 - key met field
    LIS_rc%met_nf(findex) = 8  

    do n=1,LIS_rc%nnest
     
       ! Check if starting hour of LIS run matches 00/12:
       if((LIS_rc%shr .ne.  0) .and. (LIS_rc%shr .ne. 12)) then
          write(LIS_logunit,*) "[ERR] GALWEM forecast type begins"
          write(LIS_logunit,*) "[ERR] at 00/12Z for a forecast window, so the "
          write(LIS_logunit,*) "[ERR] 'Starting hour:' should be set to 0/12 in"
          write(LIS_logunit,*) "[ERR]  your lis.config file.."
          call LIS_endrun()
       endif
      
       ! Allocate and initialize GALWEM-GE metforcing data structures:
       LIS_rc%met_nensem(findex) = galwemge_struc(n)%max_ens_members
 
       allocate(galwemge_struc(n)%metdata1(LIS_rc%met_nf(findex),&
                galwemge_struc(n)%max_ens_members,LIS_rc%ngrid(n)))
       allocate(galwemge_struc(n)%metdata2(LIS_rc%met_nf(findex),&
                galwemge_struc(n)%max_ens_members,LIS_rc%ngrid(n)))

       ! Initialize the forecast initial date-time and grib record:
       galwemge_struc(n)%init_yr = LIS_rc%syr
       galwemge_struc(n)%init_mo = LIS_rc%smo
       galwemge_struc(n)%init_da = LIS_rc%sda
       galwemge_struc(n)%init_hr = LIS_rc%shr

       galwemge_struc(n)%fcst_hour = 0
       galwemge_struc(n)%metdata1 = 0
       galwemge_struc(n)%metdata2 = 0
       gridDesci = 0
 
       gridDesci(n,1)  = 0
       gridDesci(n,2)  = galwemge_struc(n)%nc !gnc
       gridDesci(n,3)  = galwemge_struc(n)%nr !gnr
       gridDesci(n,4)  = -89.750    !lat(1,1)
       gridDesci(n,5)  = -179.750   !lon(1,1)
       gridDesci(n,6)  = 128
       gridDesci(n,7)  = 89.750     !lat(gnc,gnr)
       gridDesci(n,8)  = 179.750    !lon(gnc,gnr)
       gridDesci(n,9)  = 0.500      !dx
       gridDesci(n,10) = 0.500      !dy
       gridDesci(n,20) = 0           !for 0 to 360?

       galwemge_struc(n)%mi = galwemge_struc(n)%nc*galwemge_struc(n)%nr
       galwemge_struc(n)%fcsttime1 = 3000.0
       galwemge_struc(n)%fcsttime2 = 0.0
    enddo

    do n=1,LIS_rc%nnest
       !Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(galwemge_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               galwemge_struc(n)%n111,galwemge_struc(n)%n121,&
               galwemge_struc(n)%n211,galwemge_struc(n)%n221,&
               galwemge_struc(n)%w111,galwemge_struc(n)%w121,&
               galwemge_struc(n)%w211,galwemge_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(galwemge_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemge_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               galwemge_struc(n)%n111,galwemge_struc(n)%n121,&
               galwemge_struc(n)%n211,galwemge_struc(n)%n221,&
               galwemge_struc(n)%w111,galwemge_struc(n)%w121,&
               galwemge_struc(n)%w211,galwemge_struc(n)%w221)

          allocate(galwemge_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemge_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemge_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemge_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemge_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemge_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemge_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemge_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               galwemge_struc(n)%n112,galwemge_struc(n)%n122,&
               galwemge_struc(n)%n212,galwemge_struc(n)%n222,&
               galwemge_struc(n)%w112,galwemge_struc(n)%w122,&
               galwemge_struc(n)%w212,galwemge_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(galwemge_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci(n,:),&
               galwemge_struc(n)%n113)
       endif
    enddo

  end subroutine init_galwemge
end module galwemge_forcingMod
