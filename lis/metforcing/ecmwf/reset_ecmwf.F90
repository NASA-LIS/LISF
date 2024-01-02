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
! !MODULE: reset_ecmwf
!  \label{reset_ecmwf}
!
! !REVISION HISTORY: 
! 2Oct2015; Hiroko Beaudoing, Initial Code
!
! !INTERFACE:
subroutine reset_ecmwf
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use ecmwf_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for ECMWF forcing.
!
!EOP
  implicit none

  real    :: gridDesci(LIS_rc%nnest,50)
  integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
  real    :: upgmt
  integer :: n

   gridDesci = 0

   do n=1,LIS_rc%nnest

     gridDesci(n,1) = 0
     gridDesci(n,2) = ecmwf_struc(n)%ncold
     gridDesci(n,3) = ecmwf_struc(n)%nrold
     gridDesci(n,4) = 90.000
     gridDesci(n,5) = -180.000
     gridDesci(n,6) = 128
     gridDesci(n,7) = -60.000
     gridDesci(n,8) = 179.750
     gridDesci(n,9) = 0.250
     gridDesci(n,10) = 0.250
     gridDesci(n,20) = 64
     ecmwf_struc(n)%mi = ecmwf_struc(n)%ncold*ecmwf_struc(n)%nrold

     yr1 = 2003     !grid update time to IFS25R1 after data gap
     mo1 = 01
     da1 = 08
     hr1 = 00
     mn1 = 0; ss1 = 0
     call LIS_date2time(ecmwf_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

     yr1 = 2006     !grid update time to IFS30R1
     mo1 = 02
     da1 = 01
     hr1 = 06
     mn1 = 0; ss1 = 0
     call LIS_date2time(ecmwf_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)

     yr1 = 2008     !grid update time to IFS33R1
     mo1 = 06
     da1 = 03
     hr1 = 06
     mn1 = 0; ss1 = 0
     call LIS_date2time(ecmwf_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)

     yr1 = 2009     !grid update time to IFS35R2
     mo1 = 03
     da1 = 10
     hr1 = 06
     mn1 = 0; ss1 = 0
     call LIS_date2time(ecmwf_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)

      yr1 = 2009     !grid update time to IFS35R3
      mo1 = 09
      da1 = 08
      hr1 = 06
      mn1 = 0; ss1 = 0
      call LIS_date2time(ecmwf_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)

      yr1 = 2010     !grid update time to IFS36R1
      mo1 = 01
      da1 = 26
      hr1 = 06
      mn1 = 0; ss1 = 0
      call LIS_date2time(ecmwf_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)

      yr1 = 2011     !grid update time to IFS37R2
      mo1 = 05
      da1 = 18
      hr1 = 06
      mn1 = 0; ss1 = 0
      call LIS_date2time(ecmwf_struc(n)%griduptime7,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1)

      ecmwf_struc(n)%gridchange1 = .true.
      ecmwf_struc(n)%gridchange2 = .true.
      ecmwf_struc(n)%gridchange3 = .true.
      ecmwf_struc(n)%gridchange4 = .true.
      ecmwf_struc(n)%gridchange5 = .true.
      ecmwf_struc(n)%gridchange6 = .true.
      ecmwf_struc(n)%gridchange7 = .true.

    enddo

end subroutine reset_ecmwf
