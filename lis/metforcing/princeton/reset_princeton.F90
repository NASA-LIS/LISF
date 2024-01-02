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
! !MODULE: reset_princeton
!  \label{reset_princeton}
!
! !REVISION HISTORY: 
! 21/08/2014 Bailing Li;initial creation 
! 
! !INTERFACE:
subroutine reset_princeton()
! !USES:
  use LIS_coreMod,          only : LIS_rc
  use princeton_forcingMod, only : princeton_struc
!
! !DESCRIPTION:
!  Routine to reset nldas2 forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     princeton_struc(n)%princetontime1 = 3000.0
     princeton_struc(n)%princetontime2 = 0.0
  enddo

end subroutine reset_princeton
