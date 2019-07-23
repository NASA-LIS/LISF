!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: reset_gdas
! \label{reset_gdas}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_gdas
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use gdas_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GDAS forcing. 
!
!EOP  
  implicit none

  real :: gridDesci(20)
  integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
  real :: upgmt
  integer :: n 
  
  do n=1,LDT_rc%nnest
     
       ! This grid is good for some time in the 1980's.
       ! Look up the exact dates.
     gridDesci = 0 
     gridDesci(1) = 4
     gridDesci(2) = 192
     gridDesci(3) = 94
     gridDesci(4) = 88.542
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) =  -88.542
     gridDesci(8) = -1.875
     gridDesci(9) = 1.875
     gridDesci(10) = 47
     gridDesci(20) = 0
     gdas_struc(n)%mi = gdas_struc(n)%nc*gdas_struc(n)%nr
     gdas_struc(n)%reset_flag = .true.
       
       ! This grid is good for some time in the 1990's.
       ! Look up the exact dates.
     yr1 = 1991
     mo1 = 01
     da1 = 01
     hr1 = 12
     mn1 = 0; ss1 = 0
     call LDT_date2time( gdas_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
     
     yr1 = 2000
     mo1 = 01
     da1 = 24
     hr1 = 12
     mn1 = 0; ss1 = 0
     call LDT_date2time( gdas_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
     
     yr1 = 2002     !grid update time ~ 0.469
     mo1 = 10
     da1 = 29
     hr1 = 12
     mn1 = 0; ss1 = 0
     call LDT_date2time(gdas_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )
     
     yr1 = 2005     !grid update time ~ 0.313
     mo1 = 05
     da1 = 31
     hr1 = 12
     mn1 = 0; ss1 = 0
     call LDT_date2time(gdas_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

     yr1 = 2010     !grid update time ~ 0.205
     mo1 = 07
     da1 = 28
     hr1 = 12
     mn1 = 0; ss1 = 0
     call LDT_date2time(gdas_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

     yr1 = 2015     !grid update time ~ 0.117
     mo1 = 01
     da1 = 14
     hr1 = 6
     mn1 = 0; ss1 = 0
     call LDT_date2time(gdas_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,da1,hr1,mn1,ss1 )

     
     gdas_struc(n)%gridchange1 = .true.
     gdas_struc(n)%gridchange2 = .true.
     gdas_struc(n)%gridchange3 = .true.
     gdas_struc(n)%gridchange4 = .true.
     gdas_struc(n)%gridchange5 = .true.
     gdas_struc(n)%gridchange6 = .true.
  enddo

end subroutine reset_gdas
