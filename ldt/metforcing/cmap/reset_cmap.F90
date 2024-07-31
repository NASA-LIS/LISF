!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !MODULE: reset_cmap
! \label{reset_cmap}
! 
! !REVISION HISTORY: 
! 25Oct2005; Sujay Kumar, Initial Code
! 
! !INTERFACE:
subroutine reset_cmap
! !USES:
  use LDT_coreMod,    only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use cmap_forcingMod
!
! !DESCRIPTION:
!  Routine to reset gridchange fields for CMAP forcing. 
!
!EOP  
  implicit none

  integer :: n 
  
  do n=1,LDT_rc%nnest

     cmap_struc(n)%gridchange1 = .true.
     cmap_struc(n)%gridchange2 = .true.
     cmap_struc(n)%gridchange3 = .true.
     cmap_struc(n)%gridchange4 = .true.
     cmap_struc(n)%gridchange5 = .true.

  enddo

end subroutine reset_cmap
