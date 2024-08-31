!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: hyssib_coldstart
! \label{hyssib_coldstart}
!
! !REVISION HISTORY:
! 15 Dec 2003: Luis-Gustavo Goncalves, Initial version
! 29 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
! 25 Aug 2007: Chuck Alonge, Updates for LIS 5.0 compliance
!
! !INTERFACE:
subroutine hyssib_coldstart()
! !USES:
  use LIS_coreMod, only: LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use LIS_logMod, only : LIS_logunit
  use hyssib_lsmMod    ! HY-SSiB tile variables
!
! !DESCRIPTION:
!  
!  This routine initializes the hyssib state variables with some 
!  predefined values uniformly for the entire domain. These initial values
!  will be overwritten by the values read from the supplied 
!  HYSSIB model restart file. 
! 
!EOP
  implicit none
  integer :: t,n

  do n=1,LIS_rc%nnest
     if ( trim(LIS_rc%startcode).eq."coldstart" ) then
        write(LIS_logunit,*)'MSG: hyssib_coldstart -- cold-starting hyssib'
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           hyssib_struc(n)%hyssib(t)%tc          =  hyssib_struc(n)%initTemp
           hyssib_struc(n)%hyssib(t)%tg          =  hyssib_struc(n)%initTemp
           hyssib_struc(n)%hyssib(t)%tsn         =  hyssib_struc(n)%initTemp
           hyssib_struc(n)%hyssib(t)%td          =  hyssib_struc(n)%initTemp
           hyssib_struc(n)%hyssib(t)%www(1)      =  hyssib_struc(n)%initSm
           hyssib_struc(n)%hyssib(t)%www(2)      =  hyssib_struc(n)%initSm
           hyssib_struc(n)%hyssib(t)%www(3)      =  hyssib_struc(n)%initSm
           hyssib_struc(n)%hyssib(t)%capac(1)    =  0.0
           hyssib_struc(n)%hyssib(t)%capac(2)    =  0.0
           hyssib_struc(n)%hyssib(t)%snow(1)     =  0.0
           hyssib_struc(n)%hyssib(t)%snow(2)     =  0.0
           hyssib_struc(n)%hyssib(t)%sgfg        =  0.0
           hyssib_struc(n)%hyssib(t)%sdens       =  0.0
        enddo
     endif
     LIS_rc%yr=LIS_rc%syr
     LIS_rc%mo=LIS_rc%smo 
     LIS_rc%da=LIS_rc%sda
     LIS_rc%hr=LIS_rc%shr
     LIS_rc%mn=LIS_rc%smn
     LIS_rc%ss=LIS_rc%sss
     
     call LIS_date2time(LIS_rc%time,LIS_rc%doy,LIS_rc%gmt,LIS_rc%yr,&
          LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn,LIS_rc%ss) 
     write(LIS_logunit,*) 'MSG: hyssib_coldstart --', & 
          'Using the specified start time ',LIS_rc%time
  enddo
end subroutine hyssib_coldstart

