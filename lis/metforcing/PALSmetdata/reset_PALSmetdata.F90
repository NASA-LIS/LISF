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
! !MODULE: reset_PALSmetdata
!  \label{reset_PALSmetdata}
!
! !REVISION HISTORY: 
! 7 Mar 2013: Sujay Kumar, initial specification
! 
! !INTERFACE:
subroutine reset_PALSmetdata()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use PALSmetdata_forcingMod, only : PALSmetdata_struc
!
! !DESCRIPTION:
!  Routine to reset PALS forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     PALSmetdata_struc(n)%fcsttime1 = 3000.0
     PALSmetdata_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine reset_PALSmetdata
