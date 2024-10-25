!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module narr_forcingMod
!BOP
! !MODULE: narr_forcingMod
! 
! !DESCRIPTION: 
!  
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_NARR      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: narr_struc
!EOP

  type, public :: narr_type_dec
     
     real             :: ts
     real*8           :: narrtime1,narrtime2
     character(len=LIS_CONST_PATH_LEN)    :: narrdir
     integer          :: nc,nr
     integer          :: nlevels

     integer, allocatable :: n111(:)
     integer, allocatable :: n121(:)
     integer, allocatable :: n211(:)
     integer, allocatable :: n221(:)
     real,    allocatable :: w111(:)
     real,    allocatable :: w121(:)
     real,    allocatable :: w211(:)
     real,    allocatable :: w221(:)

     integer          :: findtime1, findtime2

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 
  end type narr_type_dec
  
  type(narr_type_dec), allocatable :: narr_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_narr
!  \label{init_narr}
!
! !REVISION HISTORY: 
! 30 APR 2009: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_narr(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
    
    integer, intent(in)  :: findex
! 
! !DESCRIPTION: 
!
!EOP
    integer   :: n
    real      :: gridDesci(50)

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the NARR forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(narr_struc(LIS_rc%nnest))

    call readcrd_narr()

    do n=1, LIS_rc%nnest
       narr_struc(n)%ts = 3*3600 
       call LIS_update_timestep(LIS_rc, n, narr_struc(n)%ts)
    enddo

!remove this warning once the met_nf variable is set.     
    print*, 'set the number of forcing variables in init_narr'
    print*, 'stopping..'
    stop

    gridDesci = 0 

    do n=1,LIS_rc%nnest

       allocate(narr_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(narr_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       
       narr_struc(n)%metdata1 = 0
       narr_struc(n)%metdata2 = 0

       gridDesci(1) = 3        ! lambert conformal grid
       gridDesci(2) = 349
       gridDesci(3) = 277
       gridDesci(4) = 1.0      ! latitude of origin
       gridDesci(5) = -145.5   ! longitude of origin
       gridDesci(6) = 128
       gridDesci(7) = 50.0     ! truelat1
       gridDesci(8) = 32.463   ! grid spacing in meters
       gridDesci(9) = 0.469
       gridDesci(10) = 50.0    ! truelat2
       gridDesci(11) = -107.00 ! standard longitude
       gridDesci(20) = 0.0       

       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 

          allocate(narr_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(narr_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(narr_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(narr_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          allocate(narr_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(narr_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(narr_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(narr_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call bilinear_interp_input(n,gridDesci,&
               narr_struc(n)%n111,&
               narr_struc(n)%n121,narr_struc(n)%n211,&
               narr_struc(n)%n221,narr_struc(n)%w111,&
               narr_struc(n)%w121,narr_struc(n)%w211,&
               narr_struc(n)%w221)

       endif
    enddo

  end subroutine init_narr

end module narr_forcingMod

