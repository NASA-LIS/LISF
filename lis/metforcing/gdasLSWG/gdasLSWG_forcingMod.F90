!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdasLSWG_forcingMod
!BOP
! !MODULE: gdasLSWG_forcingMod
! 
! !DESCRIPTION: 
!  
!  This module contains variables and subroutines used for the implementation 
!  of the forcing data generated from the GDAS output for the LSWG activity. 
!
! !USES: 
  use ESMF

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDASLSWG      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdasLSWG_struc
!EOP

  type, public :: gdasLSWG_type_dec

     real                     :: ts
     real*8                   :: time1, time2
     character*100            :: gdasLSWGfile
     integer                  :: mi
     logical                  :: startRead
     real, allocatable        :: rlat3(:)
     real, allocatable        :: rlon3(:)
     integer, allocatable     :: n113(:)
     type(ESMF_Time)          :: startTime
     type(ESMF_Time)          :: btime1, btime2
     real*8                   :: st_real
     type(ESMF_TimeInterval)  :: timeStep
     real,  allocatable       :: tmp(:,:,:)
     real,  allocatable       :: rh(:,:,:)
     real                     :: gridDesci(50)
     integer                  :: findtime1, findtime2
  end type gdasLSWG_type_dec
  
  type(gdasLSWG_type_dec), allocatable :: gdasLSWG_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_gdasLSWG
!  \label{init_gdasLSWG}
!
! !REVISION HISTORY: 
! 30 APR 2009: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_gdasLSWG(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_calendar, LIS_date2time, LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun, LIS_verify

    implicit none
    
    integer,  intent(in)     :: findex
! 
! !DESCRIPTION: 
!   This routine performs required initializations, including reading the 
!   configurable options, allocating memory and initializing data structures
!   and setting up interpolation weights. 
!EOP
    integer              :: n
    integer              :: doy
    real                 :: gmt
    integer              :: status

    allocate(gdasLSWG_struc(LIS_rc%nnest))

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GDAS-LSWG forcing reader'
       write(LIS_logunit,*) '[ERR]  is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR]  May be added in future releases.'
       write(LIS_logunit,*) '[ERR]  LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    do n=1, LIS_rc%nnest
       gdasLSWG_struc(n)%ts = 6*3600 
       call LIS_update_timestep(LIS_rc, n, gdasLSWG_struc(n)%ts)
    enddo

    do n=1,LIS_rc%nnest
       gdasLSWG_struc(n)%gridDesci = 0 
    enddo

!remove this warning once the met_nf variable is set.     
    print*, 'set the number of forcing variables in init_gdasLSWG'
    print*, 'stopping..'
    stop

    call readcrd_gdasLSWG()

    do n=1,LIS_rc%nnest

       gdasLSWG_struc(n)%gridDesci(1) = 0 
       gdasLSWG_struc(n)%gridDesci(2) = nint((gdasLSWG_struc(n)%gridDesci(8)-&
            gdasLSWG_struc(n)%gridDesci(5))/gdasLSWG_struc(n)%gridDesci(9))+1
       gdasLSWG_struc(n)%gridDesci(3) = nint((gdasLSWG_struc(n)%gridDesci(7)-&
            gdasLSWG_struc(n)%gridDesci(4))/gdasLSWG_struc(n)%gridDesci(10))+1
       gdasLSWG_struc(n)%gridDesci(6) = 128
       gdasLSWG_struc(n)%gridDesci(20) = 64

       gdasLSWG_struc(n)%mi = nint(gdasLSWG_struc(n)%gridDesci(2))*&
            nint(gdasLSWG_struc(n)%gridDesci(3))

       allocate(gdasLSWG_struc(n)%tmp(4504,gdasLSWG_struc(n)%mi,26))
       allocate(gdasLSWG_struc(n)%rh(4504,gdasLSWG_struc(n)%mi,21))
       
       gdasLSWG_struc(n)%tmp = LIS_rc%udef
       gdasLSWG_struc(n)%rh = LIS_rc%udef

       if(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
          allocate(gdasLSWG_struc(n)%rlat3(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasLSWG_struc(n)%rlon3(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(gdasLSWG_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          gdasLSWG_struc(n)%rlat3 = LIS_domain(n)%lat(:)
          gdasLSWG_struc(n)%rlon3 = LIS_domain(n)%lon(:)

          call neighbor_interp_input(n,gdasLSWG_struc(n)%gridDesci(:),&
               gdasLSWG_struc(n)%n113)
       else
          write(LIS_logunit,*) 'This interpolation option is not supported for GDAS LSWG'
          write(LIS_logunit,*) 'Program stopping ... '
          call LIS_endrun()
       endif

!      data interval is 6 hours
       call ESMF_TimeIntervalSet(gdasLSWG_struc(n)%timeStep, s=6*60*60,rc=status)

       !julian day 183, 2004
        call ESMF_TimeSet(gdasLSWG_struc(n)%startTime, yy=2004, &
             mm = 7, dd = 1, h=0, m = 0, s=0, calendar=LIS_calendar, &
             rc=status)
        call LIS_verify(status, 'ESMF_TimeSet : init_gdasLSWG')
       call LIS_date2time(gdasLSWG_struc(n)%st_real, doy, gmt, &
            2004, 7, 1, 0, 0, 0)
       
     enddo

  end subroutine init_gdasLSWG
end module gdasLSWG_forcingMod

