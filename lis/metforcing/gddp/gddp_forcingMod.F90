!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gddp_forcingMod
!BOP
! !MODULE: gddp_forcingMod
!
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data from the NASA Earth
!  Exchange Global Daily Downscaled Projections (NEX-GDDP)  
!
!  https://www.nasa.gov/nex/gddp
!  https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6
!
! ! !REVISION HISTORY:
!
! 03 Feb 2022; Sujay Kumar; Initial Specification

! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDDP      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gddp_struc
!EOP

  type, public ::  gddp_type_dec
     real     :: ts
     integer  :: nc, nr, vector_len   
     real*8   :: gddptime1,gddptime2
     character(len=LIS_CONST_PATH_LEN) :: odir
     character(len=LIS_CONST_PATH_LEN) :: scenario

     character(len=LIS_CONST_PATH_LEN) :: ref_dclimodir
     character(len=LIS_CONST_PATH_LEN) :: ref_hclimodir

     integer, allocatable   :: gindex(:,:)

     integer :: mi
     integer                :: day_check1
     integer                :: day_check2
     
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
     real, allocatable      ::  w112(:,:),w122(:,:)
     real, allocatable      ::  w212(:,:),w222(:,:)

     integer, allocatable   :: n113(:)
     
     integer                :: findtime1, findtime2
     real, allocatable      :: metdata1(:,:,:) 
     real, allocatable      :: metdata2(:,:,:)
     real, allocatable      :: metdata(:,:,:) 

     integer                :: nmodels   

     real, allocatable      :: tair_climo(:,:)
     real, allocatable      :: qair_climo(:,:)
     real, allocatable      :: swdown_climo(:,:)
     real, allocatable      :: lwdown_climo(:,:)
     real, allocatable      :: wind_climo(:,:)
     real, allocatable      :: psurf_climo(:,:)
     real, allocatable      :: prcp_climo(:,:)
  

  end type gddp_type_dec
  
  type(gddp_type_dec), allocatable :: gddp_struc(:)
!EOP
contains
  
!BOP
!
! !ROUTINE: init_GDDP
! \label{init_GDDP}
! 
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_GDDP(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !USES: 
    integer, intent(in)  :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GDDP
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!EOP


    integer :: n 
    real  :: gridDesci(LIS_rc%nnest, 50)

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GDDP forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif
    
    allocate(gddp_struc(LIS_rc%nnest))
    call readcrd_gddp()

    do n=1, LIS_rc%nnest
       gddp_struc(n)%nmodels = 25
       gddp_struc(n)%ts = 3600
       call LIS_update_timestep(LIS_rc, n, gddp_struc(n)%ts)
    enddo

    gridDesci = 0 

    LIS_rc%met_nf(findex) = 7

    gddp_struc(:)%nc = 1440
    gddp_struc(:)%nr = 600

    do n=1,LIS_rc%nnest

       allocate(gddp_struc(n)%metdata1(gddp_struc(n)%nmodels,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gddp_struc(n)%metdata2(gddp_struc(n)%nmodels,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gddp_struc(n)%metdata(gddp_struc(n)%nmodels,&
            LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gddp_struc(n)%metdata1 = 0
       gddp_struc(n)%metdata2 = 0
       gddp_struc(n)%metdata = 0
       
       gddp_struc(n)%day_check1 = -1
       gddp_struc(n)%day_check2 = -1
       
       gridDesci(n,1) = 0
       gridDesci(n,2) = gddp_struc(n)%nc
       gridDesci(n,3) = gddp_struc(n)%nr
       gridDesci(n,4) = -59.875
       gridDesci(n,5) = -179.875
       gridDesci(n,7) = 89.875
       gridDesci(n,8) = 179.875
       gridDesci(n,6) = 128
       gridDesci(n,9) = 0.25
       gridDesci(n,10) = 0.25
       gridDesci(n,20) = 0.0
       gddp_struc(n)%mi = gddp_struc(n)%nc*gddp_struc(n)%nr
       gddp_struc(n)%gddptime1 = 3000.0
       gddp_struc(n)%gddptime2 = 0.0

!Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(gddp_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               gddp_struc(n)%n111,gddp_struc(n)%n121,&
               gddp_struc(n)%n211,gddp_struc(n)%n221,&
               gddp_struc(n)%w111,gddp_struc(n)%w121,&
               gddp_struc(n)%w211,gddp_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

          allocate(gddp_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gddp_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               gddp_struc(n)%n111,gddp_struc(n)%n121,&
               gddp_struc(n)%n211,gddp_struc(n)%n221,&
               gddp_struc(n)%w111,gddp_struc(n)%w121,&
               gddp_struc(n)%w211,gddp_struc(n)%w221)

          allocate(gddp_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gddp_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gddp_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gddp_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gddp_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gddp_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gddp_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(gddp_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               gddp_struc(n)%n112,gddp_struc(n)%n122,&
               gddp_struc(n)%n212,gddp_struc(n)%n222,&
               gddp_struc(n)%w112,gddp_struc(n)%w122,&
               gddp_struc(n)%w212,gddp_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 

          allocate(gddp_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci(n,:),&
               gddp_struc(n)%n113)
       endif
       allocate(gddp_struc(n)%tair_climo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(gddp_struc(n)%qair_climo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(gddp_struc(n)%swdown_climo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(gddp_struc(n)%lwdown_climo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(gddp_struc(n)%wind_climo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(gddp_struc(n)%psurf_climo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(gddp_struc(n)%prcp_climo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       
    enddo
    
  end subroutine init_GDDP
end module gddp_forcingMod
