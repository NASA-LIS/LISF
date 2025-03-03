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
!BOP
!
! !ROUTINE: setup_MODIS_RT_lai
! \label{setup_MODIS_RT_lai}
!
! !REVISION HISTORY:
!  16 Jul 2008: Sujay Kumar; Initial Specification
!  17 Oct 2011: Yudong Tian; Modified from gfrac to lai 
!
! !INTERFACE:
subroutine setup_MODIS_RT_lai(n)
! !USES:
  use LIS_vegDataMod,  only : LIS_lai
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod, only: LIS_calendar
  use ESMF

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n
  real :: gridDesci(50)
  integer :: rc
!  type(ESMF_Time), intent(out)     :: time_begin

! !DESCRIPTION:
!  This subroutine sets up the variables required for reading the LAI 
!  climatology data from MODIS
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!EOP      
!  type(ESMF_Time)   :: mytime
!  type(ESMF_TimeInterval) :: deltaT
!  integer :: it, status
! MOIDS RT is 8-day data

!  call ESMF_TimeSet(mytime, yy=LIS_rc%yr, &
!       mm=LIS_rc%mo, &
!       dd=LIS_rc%da, &
!       h =0, &
!       m =0, &
!       s =0, &
!       calendar = LIS_calendar, &
!       rc = status)
!
!  it=8  ! look for file within 8 days prior as 8 day product
!  call ESMF_TimeIntervalSet(deltaT,d=it,rc=status)
!
!  time_begin = mytime-deltaT
!
!!  Grid
!  gridDesci = 0
!  gridDesci(1) = 0
!  gridDesci(2) = 1440
!  gridDesci(3) = 600
!  gridDesci(4) = -59.875
!  gridDesci(5) = -179.875
!  gridDesci(6) = 128         
!  gridDesci(7) = 89.875
!  gridDesci(8) = 179.875
!  gridDesci(9) = 0.25
!  gridDesci(10) = 0.25
!  gridDesci(20) = 255  ! need to check; guess here
!
!  allocate(LIS_lai(n)%rlat1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%rlon1(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!  allocate(LIS_lai(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
!
!  call bilinear_interp_input( n, gridDesci,
!       LIS_lai(n)%n111, LIS_lai(n)%n121, &
!       LIS_lai(n)%n211, LIS_lai(n)%n221, &
!       LIS_lai(n)%w111, LIS_lai(n)%w121, &
!       LIS_lai(n)%w211, LIS_lai(n)%w221 )

  call ESMF_ConfigFindLabel(LIS_config,&
       "MODIS RT LAI data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config, &
       LIS_lai(n)%laifile,rc=rc)
  call LIS_verify(rc,'MODIS RT LAI data directory: not defined')

end subroutine setup_MODIS_RT_lai
