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
! !MODULE: reset_imerg
!  \label{reset_imerg}
!
! !REVISION HISTORY: 
! Feb 15 2021  Wanshu Nie: Added to support Parameter Estimation 
! 
! !INTERFACE:
subroutine reset_imerg
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use imerg_forcingMod, only : imerg_struc
!
! !DESCRIPTION:
!  Routine to reset nldas2 forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n

  do n=1,LIS_rc%nnest
     imerg_struc(n)%imergtime = 0.0
  enddo

end subroutine reset_imerg
