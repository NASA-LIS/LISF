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
! !MODULE: reset_gdasT1534
! \label{reset_gdasT1534}
! 
! !REVISION HISTORY: 
!  20 June 2014: Sujay Kumar; initial implementation
! 
! !INTERFACE:
subroutine reset_gdasT1534
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use gdasT1534_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GDAST1534 forcing. 
!
!EOP  
  implicit none

  real :: gridDesci(50)
  integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
  real :: upgmt
  integer :: n 
  
  do n=1,LIS_rc%nnest

     gridDesci = 0 
     gridDesci(1) = 4
     gridDesci(2) = 3072
     gridDesci(3) = 1536
     gridDesci(4) = 89.91153
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) =  -89.91153
     gridDesci(8) = -0.117187
     gridDesci(9) = 0.117187
     gridDesci(10) = 768
     gridDesci(20) = 0
     gdasT1534_struc(n)%mi = gdasT1534_struc(n)%ncold*gdasT1534_struc(n)%nrold
       
  enddo
end subroutine reset_gdasT1534
