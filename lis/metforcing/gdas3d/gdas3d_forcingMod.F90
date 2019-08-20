!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdas3d_forcingMod
!BOP
! !MODULE: gdas3d_forcingMod
! 
! !DESCRIPTION: 
!  
!
! !USES: 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDAS3D      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdas3d_struc
!EOP

  type, public :: gdas3d_type_dec
     real             :: ts
     real*8           :: gdastime1,gdastime2
     character*100    :: gdasdir
     integer          :: nc,nr
     integer          :: nlayer

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
  end type gdas3d_type_dec
  
  type(gdas3d_type_dec), allocatable :: gdas3d_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_GDAS3d
!  \label{init_GDAS3d}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_GDAS3d(findex)
! !USES: 
   use LIS_coreMod,    only : LIS_rc, LIS_domain
   use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
   use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
    
    integer,  intent(in)     :: findex
! 
! !DESCRIPTION: 
!
!EOP
    integer              :: n
    real                 :: gridDesci(50)

    allocate(gdas3d_struc(LIS_rc%nnest))

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GDAS-3D forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    call readcrd_gdas3d()

    do n=1, LIS_rc%nnest
       gdas3d_struc(n)%ts = 6*3600 
       call LIS_update_timestep(LIS_rc, n, gdas3d_struc(n)%ts)
    enddo

    gridDesci = 0 

!    LIS_rc%met_nf(findex) = ??
!remove this warning once the met_nf variable is set.     
    print*, 'set the number of forcing variables in init_gdas3d'
    print*, 'stopping..'
    stop

    do n=1,LIS_rc%nnest
       
       allocate(gdas3d_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gdas3d_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gdas3d_struc(n)%metdata1 = 0
       gdas3d_struc(n)%metdata2 = 0

       gdas3d_struc(n)%findtime1 = 0 
       gdas3d_struc(n)%findtime2 = 0 

       gridDesci(1) = 4
       gridDesci(2) = 768
       gridDesci(3) = 384
       gridDesci(4) = 89.642
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) = -89.642
       gridDesci(8) = -0.469
       gridDesci(9) = 0.469
       gridDesci(10) = 192
       gridDesci(20) = 0.0       

       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
          allocate(gdas3d_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdas3d_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdas3d_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdas3d_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          allocate(gdas3d_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdas3d_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdas3d_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdas3d_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call bilinear_interp_input(n,gridDesci,&
               gdas3d_struc(n)%n111,&
               gdas3d_struc(n)%n121,gdas3d_struc(n)%n211,&
               gdas3d_struc(n)%n221,gdas3d_struc(n)%w111,&
               gdas3d_struc(n)%w121,gdas3d_struc(n)%w211,&
               gdas3d_struc(n)%w221)

       endif
    enddo

  end subroutine init_GDAS3d
end module gdas3d_forcingMod

