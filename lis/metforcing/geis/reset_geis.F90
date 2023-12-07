!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: reset_geis
!  \label{reset_geis}
!
! !REVISION HISTORY: 
! 15 Nov 2023; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_geis()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use geis_forcingMod, only : geis_struc
!
! !DESCRIPTION:
!  Routine to reset geis forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     geis_struc(n)%geistime1 = 3000.0
     geis_struc(n)%geistime2 = 0.0
  enddo

end subroutine reset_geis
