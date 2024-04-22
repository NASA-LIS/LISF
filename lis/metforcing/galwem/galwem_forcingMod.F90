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
module galwem_forcingMod
!BOP
! !MODULE: galwem_forcingMod

! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of GALWEM (17km) forecast data used as forcing
!  within LIS.

! REVISION HISTORY:
! 11 Mar 2022; Yeosang Yoon; Initial Specification
! 08 Sep 2022; Yeosang Yoon, Add codes to read GALWEM 25 DEG dataset
! 11 Jan 2024; Eric Kemp, added third entries for fcsttime and metdata
!              for temporary storage.
! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_galwem      !defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: galwem_struc

  type, public ::  galwem_type_dec
     real                              :: ts
     integer                           :: nc, nr, vector_len   
     real*8                            :: fcsttime1,fcsttime2,fcsttime3
     character(len=LIS_CONST_PATH_LEN) :: odir      !GALWEM forecast forcing Directory
     character*20                      :: runmode
     integer                           :: resol     !GALWEM forecast resolution (17km or 25deg)

     integer, allocatable   :: gindex(:,:)
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
     integer                :: fcst_hour
     integer                :: init_yr, init_mo, init_da, init_hr
     real, allocatable      :: metdata1(:,:) 
     real, allocatable      :: metdata2(:,:)
     real, allocatable      :: metdata3(:,:)
     integer                :: nmodels   

  end type galwem_type_dec

  type(galwem_type_dec), allocatable :: galwem_struc(:)
!EOP
contains

!BOP
!
! !ROUTINE: init_galwem
! \label{init_galwem}
! 
! !INTERFACE:
  subroutine init_galwem(findex)
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
    
    write(LIS_logunit,*) "[INFO] Initializing the GALWEM forecast inputs "

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GALWEM forecast forcing reader'
       write(LIS_logunit,*) '[ERR] is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR] LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(galwem_struc(LIS_rc%nnest))
    call readcrd_galwem()

    do n=1, LIS_rc%nnest
       galwem_struc(n)%ts = 3600  !check
       call LIS_update_timestep(LIS_rc, n, galwem_struc(n)%ts)
    enddo

    do n=1, LIS_rc%nnest
       if(galwem_struc(n)%resol == 17) then ! galwem-17km
          galwem_struc(:)%nc = 1536  
          galwem_struc(:)%nr = 1152
       elseif(galwem_struc(n)%resol == 25) then ! galwem-25deg
          galwem_struc(:)%nc = 1440
          galwem_struc(:)%nr = 721
       else
          write(LIS_logunit,*) '[ERR] Currently the GALWEM forecast forcing reader'
          write(LIS_logunit,*) '[ERR] supports 17 km and 25 deg datasets.'
          write(LIS_logunit,*) '[ERR] LIS forecast run-time ending.'
          call LIS_endrun()
       endif
    enddo

    ! 8 - key met field
    LIS_rc%met_nf(findex) = 8  

    do n=1,LIS_rc%nnest
     
       ! Check if starting hour of LIS run matches 00/06/12/18Z:
       if((LIS_rc%shr .ne.  0) .and. (LIS_rc%shr .ne.  6) .and. &
          (LIS_rc%shr .ne. 12) .and. (LIS_rc%shr .ne. 18)) then
          write(LIS_logunit,*) "[ERR] GALWEM forecast type begins"
          write(LIS_logunit,*) "[ERR] at 00/06/12/18Z for a forecast window, so the "
          write(LIS_logunit,*) "[ERR] 'Starting hour:' should be set to 0/6/12/18 in"
          write(LIS_logunit,*) "[ERR]  your lis.config file.."
          call LIS_endrun()
       endif
       
       allocate(galwem_struc(n)%metdata1(LIS_rc%met_nf(findex),LIS_rc%ngrid(n)))
       allocate(galwem_struc(n)%metdata2(LIS_rc%met_nf(findex),LIS_rc%ngrid(n)))
       allocate(galwem_struc(n)%metdata3(LIS_rc%met_nf(findex),LIS_rc%ngrid(n)))

       ! Initialize the forecast initial date-time and grib record:
       galwem_struc(n)%init_yr = LIS_rc%syr
       galwem_struc(n)%init_mo = LIS_rc%smo
       galwem_struc(n)%init_da = LIS_rc%sda
       galwem_struc(n)%init_hr = LIS_rc%shr

       galwem_struc(n)%fcst_hour = 0
       galwem_struc(n)%metdata1 = 0
       galwem_struc(n)%metdata2 = 0
       galwem_struc(n)%metdata3 = 0
       gridDesci = 0

       if(galwem_struc(n)%resol == 17) then   !galwem-17km 
          gridDesci(n,1) = 0
          gridDesci(n,2) = galwem_struc(n)%nc !gnc
          gridDesci(n,3) = galwem_struc(n)%nr !gnr
          gridDesci(n,4) = -89.921875         !lat(1,1)
          gridDesci(n,5) = -179.882813        !lon(1,1)
          gridDesci(n,6) = 128
          gridDesci(n,7) = 89.921875          !lat(gnc,gnr)
          gridDesci(n,8) = 179.882813         !lon(gnc,gnr)
          gridDesci(n,9) = 0.234375           !dx
          gridDesci(n,10) = 0.15625           !dy
          gridDesci(n,20) = 0
       endif

       if(galwem_struc(n)%resol == 25) then   !galwem-25deg
          gridDesci(n,1) = 0
          gridDesci(n,2) = galwem_struc(n)%nc !gnc
          gridDesci(n,3) = galwem_struc(n)%nr !gnr
          gridDesci(n,4) = -90.0              !lat(1,1)
          gridDesci(n,5) = -180.0             !lon(1,1)
          gridDesci(n,6) = 128
          gridDesci(n,7) = 90.0               !lat(gnc,gnr)
          gridDesci(n,8) = 179.75             !lon(gnc,gnr)
          gridDesci(n,9) = 0.25               !dx
          gridDesci(n,10) = 0.25              !dy
          gridDesci(n,20) = 0
       endif

       galwem_struc(n)%mi = galwem_struc(n)%nc*galwem_struc(n)%nr
       galwem_struc(n)%fcsttime1 = 3000.0
       galwem_struc(n)%fcsttime2 = 0.0
       galwem_struc(n)%fcsttime3 = 0.0
    enddo

    do n=1,LIS_rc%nnest
       !Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(galwem_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               galwem_struc(n)%n111,galwem_struc(n)%n121,&
               galwem_struc(n)%n211,galwem_struc(n)%n221,&
               galwem_struc(n)%w111,galwem_struc(n)%w121,&
               galwem_struc(n)%w211,galwem_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(galwem_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwem_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               galwem_struc(n)%n111,galwem_struc(n)%n121,&
               galwem_struc(n)%n211,galwem_struc(n)%n221,&
               galwem_struc(n)%w111,galwem_struc(n)%w121,&
               galwem_struc(n)%w211,galwem_struc(n)%w221)

          allocate(galwem_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwem_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwem_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwem_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwem_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwem_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwem_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwem_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               galwem_struc(n)%n112,galwem_struc(n)%n122,&
               galwem_struc(n)%n212,galwem_struc(n)%n222,&
               galwem_struc(n)%w112,galwem_struc(n)%w122,&
               galwem_struc(n)%w212,galwem_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(galwem_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci(n,:),&
               galwem_struc(n)%n113)
       endif
    enddo

  end subroutine init_galwem
end module galwem_forcingMod
