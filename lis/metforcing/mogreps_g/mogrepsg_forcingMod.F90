!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module mogrepsg_forcingMod
!BOP
! !MODULE: mogrepsg_forcingMod

! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of MOGREPS-G 8-day forecast data used as forcing
!  within LIS.

! REVISION HISTORY:
! 26 Jan 2023; Yeosang Yoon; Initial Specification
! 01 Jan 2024; Yeosang Yoon; update codes for precpi. bias-correction

! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_mogrepsg      !defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: mogrepsg_struc

  type, public ::  mogrepsg_type_dec
     real                              :: ts
     integer                           :: nc, nr, vector_len
     real*8                            :: fcsttime1,fcsttime2
     character(len=LIS_CONST_PATH_LEN) :: odir     !MOGREPS-G forecast forcing Directory
     character*20                      :: runmode
     integer                           :: max_ens_members

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
     real, allocatable      :: metdata1(:,:,:)
     real, allocatable      :: metdata2(:,:,:)

     integer                :: nmodels

     ! only for v-wind due to difference resolution
     integer                :: nrv
     integer, allocatable   :: nv111(:)
     integer, allocatable   :: nv121(:)
     integer, allocatable   :: nv211(:)
     integer, allocatable   :: nv221(:)
     real, allocatable      :: wv111(:), wv121(:)
     real, allocatable      :: wv211(:), wv221(:)

     integer, allocatable   :: nv112(:,:)
     integer, allocatable   :: nv122(:,:)
     integer, allocatable   :: nv212(:,:)
     integer, allocatable   :: nv222(:,:)
     real, allocatable      :: wv112(:,:), wv122(:,:)
     real, allocatable      :: wv212(:,:), wv222(:,:)

     integer, allocatable   :: nv113(:)

     ! precipitation bias correction
     integer                           :: bc        !option for bias correction
     character(len=LIS_CONST_PATH_LEN) :: cdf_fname !MOGREPS-G model CDF file name
     real, allocatable                 :: pcp_bc(:,:)
     real, allocatable                 :: bc_param_a(:,:)
     real, allocatable                 :: bc_param_b(:,:)
     real, allocatable                 :: bc_mean(:,:)
     real, allocatable                 :: bc_std(:,:)

  end type mogrepsg_type_dec

  type(mogrepsg_type_dec), allocatable :: mogrepsg_struc(:)
!EOP
contains

!BOP
!
! !ROUTINE: init_mogrepsg
! \label{init_mogrepsg}
!
! !INTERFACE:
  subroutine init_mogrepsg(findex)
! !USES:
    use LIS_coreMod,    only : LIS_rc
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !USES:
    integer, intent(in)  :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for MOGREPS-G
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!EOP

    integer :: n
    real    :: gridDesci(LIS_rc%nnest,50)
    real    :: gridDesci_v(LIS_rc%nnest,50) ! v-wind

    external :: readcrd_mogrepsg
    external :: bilinear_interp_input
    external :: conserv_interp_input
    external :: neighbor_interp_input
    external :: get_cdf_params

    write(LIS_logunit,*) &
         "[INFO] Initializing the MOGREPS-G forecast inputs "

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) &
            '[ERR] Currently the MOGREPS-G forecast forcing reader'
       write(LIS_logunit,*) &
            '[ERR] is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR] LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(mogrepsg_struc(LIS_rc%nnest))
    call readcrd_mogrepsg()

    do n=1, LIS_rc%nnest
       mogrepsg_struc(n)%ts = 3600*3  !check
       call LIS_update_timestep(LIS_rc, n, mogrepsg_struc(n)%ts)
    enddo

    do n=1, LIS_rc%nnest
       mogrepsg_struc(:)%nc  = 1280  ! mogreps-g
       mogrepsg_struc(:)%nr  = 960
       mogrepsg_struc(:)%nrv = 961   ! v-wind
    enddo

    ! 8 - key met field
    LIS_rc%met_nf(findex) = 8

    do n = 1, LIS_rc%nnest

       ! Check if starting hour of LIS run matches 00/06/12/18 UTC:
       if((LIS_rc%shr .ne.  0) .and. (LIS_rc%shr .ne. 6) .and. &
          (LIS_rc%shr .ne.  12) .and. (LIS_rc%shr .ne. 18)) then
          write(LIS_logunit,*) "[ERR] MOGREPS-G forecast type begins"
          write(LIS_logunit,*) "[ERR] at 00/12Z for a forecast window, so the "
          write(LIS_logunit,*) "[ERR] 'Starting hour:' should be set to 0/12 in"
          write(LIS_logunit,*) "[ERR]  your lis.config file.."
          call LIS_endrun()
       endif

       ! Allocate and initialize MOGREPS-G metforcing data structures:
       LIS_rc%met_nensem(findex) = mogrepsg_struc(n)%max_ens_members

       allocate(mogrepsg_struc(n)%metdata1(LIS_rc%met_nf(findex),&
                mogrepsg_struc(n)%max_ens_members,LIS_rc%ngrid(n)))
       allocate(mogrepsg_struc(n)%metdata2(LIS_rc%met_nf(findex),&
                mogrepsg_struc(n)%max_ens_members,LIS_rc%ngrid(n)))

       ! Initialize the forecast initial date-time and grib record:
       mogrepsg_struc(n)%init_yr = LIS_rc%syr
       mogrepsg_struc(n)%init_mo = LIS_rc%smo
       mogrepsg_struc(n)%init_da = LIS_rc%sda
       mogrepsg_struc(n)%init_hr = LIS_rc%shr

       mogrepsg_struc(n)%fcst_hour = 0
       mogrepsg_struc(n)%metdata1 = 0
       mogrepsg_struc(n)%metdata2 = 0
       gridDesci = 0

       gridDesci(n,1)  = 0
       gridDesci(n,2)  = real(mogrepsg_struc(n)%nc) !gnc
       gridDesci(n,3)  = real(mogrepsg_struc(n)%nr) !gnr
       gridDesci(n,4)  =  -89.906250   !lat(1,1)
       !NOTE:  gfortran complains about non-significant digits in below
       !assignment.  For now we ignore the warning message from the
       !compiler.
       gridDesci(n,5)  = -179.859375   !lon(1,1)
       gridDesci(n,6)  = 128
       gridDesci(n,7)  =  89.906250    !lat(gnc,gnr)
       !NOTE:  gfortran complains about non-significant digits in below
       !assignment.  For now we ignore the warning message from the
       !compiler.
       gridDesci(n,8)  = 179.859375    !lon(gnc,gnr)
       gridDesci(n,9)  = 0.28125       !dx
       gridDesci(n,10) = 0.18750       !dy
       gridDesci(n,20) = 0             !for 0 to 360?

       mogrepsg_struc(n)%mi = mogrepsg_struc(n)%nc*mogrepsg_struc(n)%nr
       mogrepsg_struc(n)%fcsttime1 = dble(3000.0)
       mogrepsg_struc(n)%fcsttime2 = dble(0.0)

       ! v-wind
       gridDesci_v(n,1)  = 0
       gridDesci_v(n,2)  = real(mogrepsg_struc(n)%nc)  !gnc
       gridDesci_v(n,3)  = real(mogrepsg_struc(n)%nrv) !gnr
       gridDesci_v(n,4)  =  -90.000000   !lat(1,1)
       !NOTE:  gfortran complains about non-significant digits in below
       !assignment.  For now we ignore the warning message from the
       !compiler.
       gridDesci_v(n,5)  = -179.859375   !lon(1,1)
       gridDesci_v(n,6)  = 128
       gridDesci_v(n,7)  =  90.000000    !lat(gnc,gnr)
       !NOTE:  gfortran complains about non-significant digits in below
       !assignment.  For now we ignore the warning message from the
       !compiler.
       gridDesci_v(n,8)  = 179.859375    !lon(gnc,gnr)
       gridDesci_v(n,9)  = 0.28125       !dx
       gridDesci_v(n,10) = 0.18750       !dy
       gridDesci(n,20) = 0             !for 0 to 360?
    enddo

    do n=1,LIS_rc%nnest
       !Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(mogrepsg_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               mogrepsg_struc(n)%n111,mogrepsg_struc(n)%n121,&
               mogrepsg_struc(n)%n211,mogrepsg_struc(n)%n221,&
               mogrepsg_struc(n)%w111,mogrepsg_struc(n)%w121,&
               mogrepsg_struc(n)%w211,mogrepsg_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(mogrepsg_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               mogrepsg_struc(n)%n111,mogrepsg_struc(n)%n121,&
               mogrepsg_struc(n)%n211,mogrepsg_struc(n)%n221,&
               mogrepsg_struc(n)%w111,mogrepsg_struc(n)%w121,&
               mogrepsg_struc(n)%w211,mogrepsg_struc(n)%w221)

          allocate(mogrepsg_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               mogrepsg_struc(n)%n112,mogrepsg_struc(n)%n122,&
               mogrepsg_struc(n)%n212,mogrepsg_struc(n)%n222,&
               mogrepsg_struc(n)%w112,mogrepsg_struc(n)%w122,&
               mogrepsg_struc(n)%w212,mogrepsg_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(mogrepsg_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci(n,:),&
               mogrepsg_struc(n)%n113)
       endif

       ! v-wind
       !Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(mogrepsg_struc(n)%nv111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%nv121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%nv211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%nv221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%wv111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%wv121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%wv211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%wv221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci_v(n,:),&
               mogrepsg_struc(n)%nv111,mogrepsg_struc(n)%nv121,&
               mogrepsg_struc(n)%nv211,mogrepsg_struc(n)%nv221,&
               mogrepsg_struc(n)%wv111,mogrepsg_struc(n)%wv121,&
               mogrepsg_struc(n)%wv211,mogrepsg_struc(n)%wv221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(mogrepsg_struc(n)%nv111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%nv121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%nv211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%nv221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%wv111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%wv121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%wv211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(mogrepsg_struc(n)%wv221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci_v(n,:),&
               mogrepsg_struc(n)%nv111,mogrepsg_struc(n)%nv121,&
               mogrepsg_struc(n)%nv211,mogrepsg_struc(n)%nv221,&
               mogrepsg_struc(n)%wv111,mogrepsg_struc(n)%wv121,&
               mogrepsg_struc(n)%wv211,mogrepsg_struc(n)%wv221)

          allocate(mogrepsg_struc(n)%nv112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%nv122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%nv212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%nv222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%wv112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%wv122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%wv212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(mogrepsg_struc(n)%wv222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci_v(n,:),&
               mogrepsg_struc(n)%nv112,mogrepsg_struc(n)%nv122,&
               mogrepsg_struc(n)%nv212,mogrepsg_struc(n)%nv222,&
               mogrepsg_struc(n)%wv112,mogrepsg_struc(n)%wv122,&
               mogrepsg_struc(n)%wv212,mogrepsg_struc(n)%wv222)

       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(mogrepsg_struc(n)%nv113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci_v(n,:),&
               mogrepsg_struc(n)%nv113)
       endif
    enddo

    ! precipitation bias correction
    do n = 1, LIS_rc%nnest
       if (mogrepsg_struc(n)%bc == 1) then
          allocate(mogrepsg_struc(n)%pcp_bc(mogrepsg_struc(n)%max_ens_members,LIS_rc%ngrid(n)))
          allocate(mogrepsg_struc(n)%bc_param_a(LIS_rc%ngrid(n),8)) !8: lead time
          allocate(mogrepsg_struc(n)%bc_param_b(LIS_rc%ngrid(n),8))
          allocate(mogrepsg_struc(n)%bc_mean(LIS_rc%ngrid(n),8))
          allocate(mogrepsg_struc(n)%bc_std(LIS_rc%ngrid(n),8))

          mogrepsg_struc(n)%pcp_bc = 0
          mogrepsg_struc(n)%bc_param_a = 0
          mogrepsg_struc(n)%bc_param_b = 0
          mogrepsg_struc(n)%bc_mean = 0
          mogrepsg_struc(n)%bc_std = 0

          ! read cdf parameters
          call get_cdf_params(n,mogrepsg_struc(n)%cdf_fname,LIS_rc%mo, &
               mogrepsg_struc(n)%bc_param_a, &
               mogrepsg_struc(n)%bc_param_b, &
               mogrepsg_struc(n)%bc_mean, mogrepsg_struc(n)%bc_std)
       endif
    enddo

  end subroutine init_mogrepsg
end module mogrepsg_forcingMod
