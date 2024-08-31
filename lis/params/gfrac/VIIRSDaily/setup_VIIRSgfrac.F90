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
! !ROUTINE: setup_VIIRSgfrac
! \label{setup_VIIRSgfrac}
!
! !REVISION HISTORY:
!  16 Jul 2008: Sujay Kumar; Initial Specification
!   4 Nov 2014: Jonathan Case: Modified for NESDIS/VIIRS GVF, 
!                              following SPORT/MODIS GVF code.
!
! !INTERFACE:
subroutine setup_VIIRSgfrac(n)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_vegDataMod,  only : LIS_gfrac
  use LIS_timeMgrMod
  use VIIRSgfracMod

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n
  integer             :: rc
! J.Case (10/27/2014)
  integer             :: k

! !DESCRIPTION:
!  This subroutine sets up the variables required for reading the greenness
!  fraction daily composite data from VIIRS
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!EOP      

! This is a daily real-time data source. Set the intervalType = 5.
! J.Case (8/10/2010) 5 = daily interval type, as defined in core/LIS_gfracMod.F90 and 
! lis.config file.
  real :: gridDesci(50)
  integer :: status

  LIS_gfrac(n)%gfracIntervalType = "daily"
  LIS_gfrac(n)%gfracInterval = 86400 

  gridDesci = 0
  
  call ESMF_ConfigGetattribute(LIS_config,LIS_gfrac(n)%realtimemode,&
       label="VIIRS GVF use realtime mode:",rc=status)
  call LIS_verify(status,'VIIRS GVF use realtime mode: not defined')

  call ESMF_ConfigGetattribute(LIS_config,gridDesci(4),&
       label="VIIRS GVF lower left lat:",rc=status)
  call LIS_verify(status,'VIIRS GVF lower left lat: not defined')

  call ESMF_ConfigGetattribute(LIS_config,gridDesci(5),&
       label="VIIRS GVF lower left lon:",rc=status)
  call LIS_verify(status,'VIIRS GVF lower left lon: not defined')


  call ESMF_ConfigGetattribute(LIS_config,gridDesci(7),&
       label="VIIRS GVF upper right lat:",rc=status)
  call LIS_verify(status,'VIIRS GVF upper right lat: not defined')


  call ESMF_ConfigGetattribute(LIS_config,gridDesci(8),&
       label="VIIRS GVF upper right lon:",rc=status)
  call LIS_verify(status,'VIIRS GVF upper right lon: not defined')

  call ESMF_ConfigGetattribute(LIS_config,gridDesci(9),&
       label="VIIRS GVF resolution (dx):",rc=status)
  call LIS_verify(status,'VIIRS GVF resolution (dx): not defined')


  call ESMF_ConfigGetattribute(LIS_config,gridDesci(10),&
       label="VIIRS GVF resolution (dy):",rc=status)
  call LIS_verify(status,'VIIRS GVF resolution (dy): not defined')
  
  gridDesci(1) = 0
  gridDesci(2) = nint((gridDesci(8)-gridDesci(5))/gridDesci(9))+1
  gridDesci(3) = nint((gridDesci(7)-gridDesci(4))/gridDesci(10))+1
  gridDesci(6) = 128         
  gridDesci(20) = 64

  allocate(LIS_gfrac(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  allocate(LIS_gfrac(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(LIS_gfrac(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  call bilinear_interp_input( n, gridDesci,  &
       LIS_gfrac(n)%n111, LIS_gfrac(n)%n121, &
       LIS_gfrac(n)%n211, LIS_gfrac(n)%n221, &
       LIS_gfrac(n)%w111, LIS_gfrac(n)%w121, &
       LIS_gfrac(n)%w211, LIS_gfrac(n)%w221 )

  call ESMF_ConfigFindLabel(LIS_config,"VIIRS greenness data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config,LIS_gfrac(n)%gfracfile,rc=rc)
  call LIS_verify(rc,'VIIRS greenness data directory: not defined')

  call LIS_registerAlarm("LIS gfrac read alarm",LIS_rc%ts, &
       LIS_gfrac(n)%gfracInterval,intervalType=LIS_gfrac(n)%gfracIntervalType) 

! J.Case (10/24/2014)
  !do k = 1,10
  !  LIS_rc%gridDesc(n,k+30) = gridDesci(k)
  !enddo
  gvf_nc = nint(gridDesci(2))
  gvf_nr = nint(gridDesci(3))
  write(LIS_logunit,*) 'GVF gridDesci: ',gridDesci
  write(LIS_logunit,*) 'GVF nc: ',gvf_nc
  write(LIS_logunit,*) 'GVF nr: ',gvf_nr
  !write(LIS_logunit,*) 'LIS gridDescr: ',LIS_rc%gridDesc(n,:)

end subroutine setup_VIIRSgfrac

