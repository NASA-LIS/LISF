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
! !MODULE: reset_plumber2
! \label{reset_plumber2}
! 
! !REVISION HISTORY: 
! October 2021: MCB@ASRC - Initial version
! (adapted from PUMET by Eric Kemp & David Mocko)
!
! !INTERFACE:
subroutine reset_plumber2
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use plumber2_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for PLUMBER2 forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LIS_rc%nnest
     plumber2_struc(n)%plumber2time1 = 3000.0
     plumber2_struc(n)%plumber2time2 = 0.0
     plumber2_struc(n)%ringtime = 0.0
     plumber2_struc(n)%reset_flag = .true.
  enddo

end subroutine reset_plumber2
