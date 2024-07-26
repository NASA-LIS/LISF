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
! !MODULE: reset_AWRAL
!  \label{reset_AWRAL}
!
! !REVISION HISTORY: 
! 30 Jan 2017; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_AWRAL()

! !USES:
  use LIS_coreMod, only : LIS_rc
  use AWRAL_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup Stage IV forcing related memory allocations.   
! 
!EOP
  implicit none
  integer   :: n

  do n=1,LIS_rc%nnest  
     AWRAL_struc(n)%AWRALtime = 0.0
  enddo

end subroutine reset_AWRAL
